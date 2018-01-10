# encoding: utf-8
'''Tools for creating and working with Arakawa C-Grids'''
__docformat__ = "restructuredtext en"

import os
import sys
import ctypes
import pickle
from warnings import warn
from copy import deepcopy

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.artist import Artist
from matplotlib.patches import Polygon, CirclePolygon
from matplotlib.lines import Line2D
#from matplotlib.numerix.mlab import amin
from matplotlib.mlab import dist_point_to_segment
#from matplotlib.nxutils import points_inside_poly      #decrepeted in version 1.3.0. Use maplotlib.path instead.


from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import pyproj

try:
    import scipy.spatial.cKDTree as KDTree
except:
    #  no scipy
    from pyroms.extern import KDTree

import pyroms
from pyroms.vgrid import *
from pyroms.extern import GreatCircle

class BoundaryInteractor(object):
    """
    Interactive grid creation
        
    bry = BoundaryClick(x=[], y=[], beta=None, ax=gca(), **gridgen_options)
    
    The initial boundary polygon points (x and y) are
    counterclockwise, starting in the upper left corner of the
    boundary. 
    
    Key commands:
        
        t : toggle visibility of verticies
        d : delete a vertex
        i : insert a vertex at a point on the polygon line
        
        p : set vertex as beta=1 (a Positive turn, marked with green triangle)
        m : set vertex as beta=1 (a Negative turn, marked with red triangle)
        z : set vertex as beta=0 (no corner, marked with a black dot)
        
        G : generate grid from the current boundary using gridgen
        T : toggle visability of the current grid
    
    Methods:
    
        bry.dump(bry_file)
            Write the current boundary informtion (bry.x, bry.y, bry.beta) to
            a cPickle file bry_file.
        
        bry.load(bry_file)
            Read in boundary informtion (x, y, beta) from the cPickle file
            bry_file.
        
        bry.remove_grid()  
            Remove gridlines from axes.
    
    Attributes:
        bry.x : the X boundary points
        bry.y : the Y boundary points
        bry.verts : the verticies of the grid
        bry.grd : the CGrid object
        
    """
    
    _showverts = True
    _showbetas = True
    _showgrid = True
    _epsilon = 5  # max pixel distance to count as a vertex hit
    
    def _update_beta_lines(self):
        """Update m/pline by finding the points where self.beta== -/+ 1"""
        x, y = list(zip(*self._poly.xy))
        num_points = len(x)-1  # the first and last point are repeated
        
        xp = [x[n] for n in range(num_points) if self.beta[n]==1]
        yp = [y[n] for n in range(num_points) if self.beta[n]==1]
        self._pline.set_data(xp, yp)
        
        xm = [x[n] for n in range(num_points) if self.beta[n]==-1]
        ym = [y[n] for n in range(num_points) if self.beta[n]==-1]
        self._mline.set_data(xm, ym)
        
        xz = [x[n] for n in range(num_points) if self.beta[n]==0]
        yz = [y[n] for n in range(num_points) if self.beta[n]==0]
        self._zline.set_data(xz, yz)
        
        if len(x)-1 < self.gridgen_options['ul_idx']:
            self.gridgen_options['ul_idx'] = len(x)-1
        xs = x[self.gridgen_options['ul_idx']]
        ys = y[self.gridgen_options['ul_idx']]
        self._sline.set_data(xs, ys)
    
    def remove_grid(self):
        """Remove a generated grid from the BoundaryClick figure"""
        if hasattr(self, '_gridlines'):
            for line in self._gridlines:
                self._ax.lines.remove(line)
            delattr(self, '_gridlines')
    
    def _draw_callback(self, event):
        self._background = self._canvas.copy_from_bbox(self._ax.bbox)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._pline)
        self._ax.draw_artist(self._mline)
        self._ax.draw_artist(self._zline)
        self._ax.draw_artist(self._sline)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    def _poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self._line.get_visible()
        Artist.update_from(self._line, poly)
        self._line.set_visible(vis)  # don't use the poly visibility state
    
    def _get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
        try:
            x, y = list(zip(*self._poly.xy))
            
            # display coords
            xt, yt = self._poly.get_transform().numerix_x_y(x, y)
            d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
            indseq = np.nonzero(np.equal(d, np.amin(d)))
            ind = indseq[0]
        
            if d[ind]>=self._epsilon:
                ind = None
        
            return ind
        except:
            # display coords
            xy = np.asarray(self._poly.xy)
            xyt = self._poly.get_transform().transform(xy)
            xt, yt = xyt[:, 0], xyt[:, 1]
            d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
            indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
            ind = indseq[0]
            
            if d[ind]>=self._epsilon:
                ind = None
            
            return ind
    
    def _button_press_callback(self, event):
        'whenever a mouse button is pressed'
        # if not self._showverts: return
        if event.inaxes==None: return
        if event.button != 1: return
        self._ind = self._get_ind_under_point(event)
    
    def _button_release_callback(self, event):
        'whenever a mouse button is released'
        # if not self._showverts: return
        if event.button != 1: return
        self._ind = None
    
    def _key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes: return
        if event.key=='shift': return
        
        if event.key=='t':
            self._showbetas = not self._showbetas
            self._line.set_visible(self._showbetas)
            self._pline.set_visible(self._showbetas)
            self._mline.set_visible(self._showbetas)
            self._zline.set_visible(self._showbetas)
            self._sline.set_visible(self._showbetas)
        elif event.key=='d':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self._poly.xy = [tup for i,tup in enumerate(self._poly.xy) \
                                 if i!=ind]
                self._line.set_data(list(zip(*self._poly.xy)))
                self.beta = [beta for i,beta in enumerate(self.beta) \
                             if i!=ind]
        elif event.key=='p':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = 1.0
        elif event.key=='m':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = -1.0
        elif event.key=='z':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.beta[ind] = 0.0
        elif event.key=='s':
            ind = self._get_ind_under_point(event)
            if ind is not None:
                self.gridgen_options['ul_idx'] = ind
        elif event.key=='i':
            xys = self._poly.get_transform().transform(self._poly.xy)
            p = event.x, event.y # display coords
            for i in range(len(xys)-1):
                s0 = xys[i]
                s1 = xys[i+1]
                d = dist_point_to_segment(p, s0, s1)
                if d<=self._epsilon:
                    self._poly.xy = np.array(
                        list(self._poly.xy[:i+1]) +
                        [(event.xdata, event.ydata)] +
                        list(self._poly.xy[i+1:]))
                    self._line.set_data(list(zip(*self._poly.xy)))
                    self.beta.insert(i+1, 0)
                    break
            s0 = xys[-1]
            s1 = xys[0]
            d = dist_point_to_segment(p, s0, s1)
            if d<=self._epsilon:
                self._poly.xy = np.array(
                    list(self._poly.xy) +
                    [(event.xdata, event.ydata)])
                self._line.set_data(list(zip(*self._poly.xy)))
                self.beta.append(0)
        elif event.key=='G' or event.key == '1':
            options = deepcopy(self.gridgen_options)
            shp = options.pop('shp')
            if self.proj is None:
                x = self.x
                y = self.y
                self.grd = Gridgen(x, y, self.beta, shp,
                                   proj=self.proj, **options)
            else:
                lon, lat = self.proj(self.x, self.y, inverse=True)
                self.grd = Gridgen(lon, lat, self.beta, shp, 
                                   proj=self.proj, **options)
            self.remove_grid()
            self._showgrid = True
            gridlineprops = {'linestyle':'-', 'color':'k', 'lw':0.2}
            self._gridlines = []
            for line in self._ax._get_lines(*(self.grd.x, self.grd.y),
                                            **gridlineprops):
                self._ax.add_line(line)
                self._gridlines.append(line)
            for line in self._ax._get_lines(*(self.grd.x.T, self.grd.y.T),
                                            **gridlineprops):
                self._ax.add_line(line)
                self._gridlines.append(line)
        elif event.key=='T' or event.key == '2':
            self._showgrid = not self._showgrid
            if hasattr(self, '_gridlines'):
                for line in self._gridlines:
                    line.set_visible(self._showgrid)
        
        self._update_beta_lines()
        self._draw_callback(event)
        self._canvas.draw()
    
    def _motion_notify_callback(self, event):
        'on mouse movement'
        # if not self._showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        x,y = event.xdata, event.ydata
        self._poly.xy[self._ind] = x, y
        if self._ind == 0:
            self._poly.xy[-1] = x, y
        
        x, y = list(zip(*self._poly.xy))
        self._line.set_data(x[:-1], y[:-1])
        self._update_beta_lines()
        
        self._canvas.restore_region(self._background)
        self._ax.draw_artist(self._poly)
        self._ax.draw_artist(self._pline)
        self._ax.draw_artist(self._mline)
        self._ax.draw_artist(self._zline)
        self._ax.draw_artist(self._sline)
        self._ax.draw_artist(self._line)
        self._canvas.blit(self._ax.bbox)
    
    
    def __init__(self, x, y=None, beta=None, ax=None, proj=None, 
                 **gridgen_options):
        
        if isinstance(x, str):
            bry_dict = np.load(x)
            x = bry_dict['x']
            y = bry_dict['y']
            beta = bry_dict['beta']
        
        assert len(x) >= 4, 'Boundary must have at least four points.'
        
        if ax is None: 
            ax = plt.gca()
        
        self._ax = ax
        
        self.proj = proj
        
        # Set default gridgen option, and copy over specified options.
        self.gridgen_options = {'ul_idx': 0, 'shp': (32, 32)}
        
        for key, value in gridgen_options.items():
            self.gridgen_options[key] = gridgen_options[key]
        
        x = list(x); y = list(y)
        assert len(x)==len(y), 'arrays must be equal length'
        
        if beta is None:
            self.beta = [0 for xi in x]
        else:
            assert len(x)==len(beta), 'beta must have same length as x and y'
            self.beta = list(beta)
        
        self._line = Line2D(x, y, animated=True, 
                            ls='-', color='k', alpha=0.5, lw=1)
        self._ax.add_line(self._line)
        
        self._canvas = self._line.figure.canvas        
        
        self._poly = Polygon(self.verts, alpha=0.1, fc='k', animated=True)
        self._ax.add_patch(self._poly)
        
        # Link in the lines that will show the beta values
        # pline for positive turns, mline for negative (minus) turns
        # otherwize zline (zero) for straight sections
        self._pline = Line2D([], [], marker='^', ms=12, mfc='g',\
                             animated=True, lw=0)
        self._mline = Line2D([], [], marker='v', ms=12, mfc='r',\
                             animated=True, lw=0)
        self._zline = Line2D([], [], marker='o', mfc='k', animated=True, lw=0)
        self._sline = Line2D([], [], marker='s', mfc='k', animated=True, lw=0)
        
        self._update_beta_lines()
        self._ax.add_artist(self._pline)
        self._ax.add_artist(self._mline)
        self._ax.add_artist(self._zline)
        self._ax.add_artist(self._sline)
        
        # get the canvas and connect the callback events
        cid = self._poly.add_callback(self._poly_changed)
        self._ind = None # the active vert
        
        self._canvas.mpl_connect('draw_event', self._draw_callback)
        self._canvas.mpl_connect('button_press_event',\
                                 self._button_press_callback)
        self._canvas.mpl_connect('key_press_event', self._key_press_callback)
        self._canvas.mpl_connect('button_release_event',\
                                 self._button_release_callback)
        self._canvas.mpl_connect('motion_notify_event',\
                                 self._motion_notify_callback)
    
    def save_bry(self, bry_file='bry.pickle'):
        f = open(bry_file, 'wb')
        bry_dict = {'x': self.x, 'y': self.y, 'beta': self.beta}
        pickle.dump(bry_dict, f, protocol=-1)
        f.close()
    
    def load_bry(self, bry_file='bry.pickle'):
        bry_dict = np.load(bry_file)
        x = bry_dict['x']
        y = bry_dict['y']
        self._line.set_data(x, y)
        self.beta = bry_dict['beta']
        if hasattr(self, '_poly'):
            self._poly.xy = list(zip(x, y))
            self._update_beta_lines()
            self._draw_callback(None)
            self._canvas.draw()
    
    def save_grid(self, grid_file='grid.pickle'):
        f = open(grid_file, 'wb')
        pickle.dump(self.grd, f, protocol=-1)
        f.close()
    
    def _get_verts(self): return list(zip(self.x, self.y))
    verts = property(_get_verts)    
    def get_xdata(self): return self._line.get_xdata()
    x = property(get_xdata)
    def get_ydata(self): return self._line.get_ydata()
    y = property(get_ydata)
    


def _approximate_erf(x):
    '''
    Return approximate solution to error function
    see http://en.wikipedia.org/wiki/Error_function
    '''
    a = -(8*(np.pi-3.0)/(3.0*np.pi*(np.pi-4.0)))
    return np.sign(x) * \
           np.sqrt(1.0 - np.exp( -x**2*(4.0/np.pi+a*x*x)/(1.0+a*x*x) ))
    

class _Focus_x(object):
    """
    Return a transformed, uniform grid, focused in the x-direction
    
    This class may be called with a uniform grid, with limits from [0, 1], to
    create a focused grid in the x-directions centered about xo. The output
    grid is also uniform from [0, 1] in both x and y.
    
    Parameters
    ----------
    xo : float
        Location about which to focus grid
    factor : float
        amount to focus grid. Creates cell sizes that are factor smaller in
        the focused
        region.
    Rx : float
        Lateral extent of focused region, similar to a lateral spatial scale
        for the focusing region.
    
    Returns
    -------
    foc : class
        The class may be called with arguments of a grid. The returned
        transformed grid (x, y) will be focused as per the parameters above.
    """
    
    def __init__(self, xo, factor=2.0, Rx=0.1):
        self.xo = xo
        self.factor = factor
        self.Rx = Rx
    
    def __call__(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        assert not np.any(x>1.0) or not np.any(x<0.0)  \
            or not np.any(y>1.0) or not np.any(x<0.0), \
                'x and y must both be within the range [0, 1].'
        
        alpha = 1.0 - 1.0/self.factor
        def xf(x):
            return x - 0.5*( np.sqrt(np.pi)*self.Rx*alpha
                            *_approximate_erf((x-self.xo)/self.Rx) )
        
        xf0 = xf(0.0); xf1 = xf(1.0)
        
        return (xf(x)-xf0)/(xf1-xf0), y

class _Focus_y(object):
    """
    Return a transformed, uniform grid, focused in the y-direction
    
    This class may be called with a uniform grid, with limits from [0, 1], 
    to create a focused grid in the y-directions centered about yo. 
    The output grid is also uniform from [0, 1] in both x and y.
    
    Parameters
    ----------
    yo : float
        Location about which to focus grid
    factor : float
        amount to focus grid. Creates cell sizes that are factor 
        smaller in the focused region.
    Ry : float
        Lateral extent of focused region, similar to a lateral 
        spatial scale for the focusing region.
    
    Returns
    -------
    foc : class
        The class may be called with arguments of a grid. The returned 
        transformed grid (x, y) will be focused as per the parameters above.
    """
    
    def __init__(self, yo, factor=2.0, Ry=0.1):
        self.yo = yo
        self.factor = factor
        self.Ry = Ry
    
    def __call__(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        assert not np.any(x>1.0) or not np.any(x<0.0)  \
            or not np.any(y>1.0) or not np.any(x<0.0), \
                'x and y must both be within the range [0, 1].'
        
        alpha = 1.0 - 1.0/self.factor
        
        def yf(y):
            return y - 0.5*( np.sqrt(np.pi)*self.Ry*alpha 
                            *_approximate_erf((y-self.yo)/self.Ry) )
        
        yf0 = yf(0.0); yf1 = yf(1.0)
        
        return x, (yf(y)-yf0)/(yf1-yf0)

class Focus(object):
    """
    Return a container for a sequence of Focus objects
    
    foc = Focus()
    
    The sequence is populated by using the 'add_focus_x' and 'add_focus_y'
    methods. These methods define a point ('xo' or 'yo'), around witch to
    focus, a focusing factor of 'focus', and x and y extent of focusing given
    by Rx or Ry. The region of focusing will be approximately Gausian, and the
    resolution will be increased by approximately the value of factor.
    
    Methods
    -------
    foc.add_focus_x(xo, factor=2.0, Rx=0.1)
    foc.add_focus_y(yo, factor=2.0, Ry=0.1)
    
    Calls to the object return transformed coordinates:
        xf, yf = foc(x, y)
    where x and y must be within [0, 1], and are typically a uniform,
    normalized grid. The focused grid will be the result of applying each of
    the focus elements in the sequence they are added to the series.
    
    
    EXAMPLES
    --------
    
    >>> foc = pyroms.grid.Focus()
    >>> foc.add_focus_x(0.2, factor=3.0, Rx=0.2)
    >>> foc.add_focus_y(0.6, factor=5.0, Ry=0.35)
    
    >>> x, y = np.mgrid[0:1:3j,0:1:3j]
    >>> xf, yf = foc(x, y)
    
    >>> print xf
    [[ 0.          0.          0.        ]
     [ 0.36594617  0.36594617  0.36594617]
     [ 1.          1.          1.        ]]
    >>> print yf
    [[ 0.          0.62479833  1.        ]
     [ 0.          0.62479833  1.        ]
     [ 0.          0.62479833  1.        ]]
    """
    def __init__(self):
        self._focuspoints = []
    
    def add_focus_x(self, xo, factor=2.0, Rx=0.1):
        """docstring for add_point"""
        self._focuspoints.append(_Focus_x(xo, factor, Rx))
    
    def add_focus_y(self, yo, factor=2.0, Ry=0.1):
        """docstring for add_point"""
        self._focuspoints.append(_Focus_y(yo, factor, Ry))
    
    def __call__(self, x, y):
        """docstring for __call__"""
        for focuspoint in self._focuspoints:
            x, y = focuspoint(x, y)
        return x, y



class CGrid(object):
    """
    Curvilinear Arakawa C-Grid
     
    The basis for the CGrid class are two arrays defining the verticies of the
    grid in Cartesian (for geographic coordinates, see CGrid_geo). An optional
    mask may be defined on the cell centers. Other Arakawa C-grid properties,
    such as the locations of the cell centers (rho-points), cell edges (u and
    v velocity points), cell widths (dx and dy) and other metrics (angle,
    dmde, and dndx) are all calculated internally from the vertex points.
     
    Input vertex arrays may be either type np.array or np.ma.MaskedArray. If
    masked arrays are used, the mask will be a combination of the specified
    mask (if given) and the masked locations.
     
    EXAMPLES:
    --------
     
    >>> x, y = mgrid[0.0:7.0, 0.0:8.0]
    >>> x = np.ma.masked_where( (x<3) & (y<3), x)
    >>> y = np.ma.MaskedArray(y, x.mask)
    >>> grd = pyroms.grid.CGrid(x, y)
    >>> print grd.x_rho
    [[-- -- -- 0.5 0.5 0.5 0.5]
     [-- -- -- 1.5 1.5 1.5 1.5]
     [-- -- -- 2.5 2.5 2.5 2.5]
     [3.5 3.5 3.5 3.5 3.5 3.5 3.5]
     [4.5 4.5 4.5 4.5 4.5 4.5 4.5]
     [5.5 5.5 5.5 5.5 5.5 5.5 5.5]]
    >>> print grd.mask
    [[ 0.  0.  0.  1.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]
     [ 1.  1.  1.  1.  1.  1.  1.]]
    """
    
    def __init__(self, x_vert, y_vert, x_rho=None, y_rho=None, x_u=None, y_u=None, x_v=None, y_v=None, \
                    x_psi=None, y_psi=None, dx=None, dy=None, dndx=None, dmde=None, angle_rho=None):
                
        assert np.ndim(x_vert)==2 and np.ndim(y_vert)==2 and np.shape(x_vert)==np.shape(y_vert), \
            'x and y must be 2D arrays of the same size.'
        
        if np.any(np.isnan(x_vert)) or np.any(np.isnan(y_vert)):
            x_vert = np.ma.masked_where( (isnan(x_vert)) | (isnan(y_vert)) , x_vert)
            y_vert = np.ma.masked_where( (isnan(x_vert)) | (isnan(y_vert)) , y_vert)
            
        self.x_vert = x_vert
        self.y_vert = y_vert

        self.f = None
        self.spherical = 'F'
        
        mask_shape = tuple([n-1 for n in self.x_vert.shape])
        self.mask_rho = np.ones(mask_shape, dtype='d')
        
        # If maskedarray is given for verticies, modify the mask such that 
        # non-existant grid points are masked.  A cell requires all four
        # verticies to be defined as a water point.
        if isinstance(self.x_vert, np.ma.MaskedArray):
            mask = (self.x_vert.mask[:-1,:-1] | self.x_vert.mask[1:,:-1] | \
                    self.x_vert.mask[:-1,1:] | self.x_vert.mask[1:,1:])
            self.mask_rho = np.asarray(~(~np.bool_(self.mask_rho) | mask), dtype='d')
        
        if isinstance(self.y_vert, np.ma.MaskedArray):
            mask = (self.y_vert.mask[:-1,:-1] | self.y_vert.mask[1:,:-1] | \
                    self.y_vert.mask[:-1,1:] | self.y_vert.mask[1:,1:])
            self.mask_rho = np.asarray(~(~np.bool_(self.mask_rho) | mask), dtype='d')

        if x_rho is None or y_rho is None or x_u is None or y_u is None or \
            x_v is None or y_v is None or x_psi is None or y_psi is None:
            self._calculate_subgrids()
        else:
            self.x_rho = x_rho
            self.y_rho = y_rho
            self.x_u = x_u
            self.y_u = y_u
            self.x_v = x_v
            self.y_v = y_v
            self.x_psi = x_psi
            self.y_psi = y_psi

        if dx is None or dy is None:
            self._calculate_metrics()
        else:
            self.dx = dx
            self.dy = dy

        self.xl = np.maximum(self.dx[0,:].sum(), self.dx[-1,:].sum())
        self.el = np.maximum(self.dy[:,0].sum(), self.dy[:,-1].sum())

        if dndx is None or dmde is None:
            self._calculate_derivative_metrics()
        else:
            self.dndx = dndx
            self.dmde = dmde

        if angle_rho is None:
            self._calculate_angle_rho()
        else:
            self.angle_rho = angle_rho

        self._calculate_angle()
   
    
    def _calculate_subgrids(self):
        self.x_rho = 0.25*(self.x_vert[1:,1:]+self.x_vert[1:,:-1]+ \
                           self.x_vert[:-1,1:]+self.x_vert[:-1,:-1])
        self.y_rho = 0.25*(self.y_vert[1:,1:]+self.y_vert[1:,:-1]+ \
                           self.y_vert[:-1,1:]+self.y_vert[:-1,:-1])
        self.x_u = 0.5*(self.x_vert[:-1,1:-1] + self.x_vert[1:,1:-1])
        self.y_u = 0.5*(self.y_vert[:-1,1:-1] + self.y_vert[1:,1:-1])
        self.x_v = 0.5*(self.x_vert[1:-1,:-1] + self.x_vert[1:-1,1:])
        self.y_v = 0.5*(self.y_vert[1:-1,:-1] + self.y_vert[1:-1,1:])
        self.x_psi = self.x_vert[1:-1,1:-1]
        self.y_psi = self.y_vert[1:-1,1:-1]
    
    def _calculate_metrics(self):
        'Calculates pm, pn, dndx, dmde from x_vert and y_vert'
        x_temp = 0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])
        y_temp = 0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])
        self.dx = np.sqrt(np.diff(x_temp, axis=1)**2 + np.diff(y_temp, axis=1)**2)
        x_temp = 0.5*(self.x_vert[:,1:]+self.x_vert[:,:-1])
        y_temp = 0.5*(self.y_vert[:,1:]+self.y_vert[:,:-1])
        self.dy = np.sqrt(np.diff(x_temp, axis=0)**2 + np.diff(y_temp, axis=0)**2)

    def _calculate_derivative_metrics(self):     
        if isinstance(self.dy, np.ma.MaskedArray):
            self.dndx = np.ma.zeros(self.x_rho.shape, dtype='d')
        else:
            self.dndx = np.zeros(self.x_rho.shape, dtype='d')
        
        if isinstance(self.dx, np.ma.MaskedArray):
            self.dmde = np.ma.zeros(self.x_rho.shape, dtype='d')
        else:
            self.dmde = np.zeros(self.x_rho.shape, dtype='d')
        
        self.dndx[1:-1,1:-1] = 0.5*(self.dy[1:-1,2:] - self.dy[1:-1,:-2])
        self.dmde[1:-1,1:-1] = 0.5*(self.dx[2:,1:-1] - self.dx[:-2,1:-1])

    def _calculate_angle(self):     
        if isinstance(self.x_vert, np.ma.MaskedArray) or \
           isinstance(self.y_vert, np.ma.MaskedArray):
            self.angle = np.ma.zeros(self.x_vert.shape, dtype='d')
        else:
            self.angle = np.zeros(self.x_vert.shape, dtype='d')
        
        angle_ud = np.arctan2(np.diff(self.y_vert, axis=1), np.diff(self.x_vert, axis=1))
        angle_lr = np.arctan2(np.diff(self.y_vert, axis=0), np.diff(self.x_vert, axis=0)) - np.pi/2.0        
        # domain center
        self.angle[1:-1,1:-1] = 0.25*(angle_ud[1:-1,1:]+angle_ud[1:-1,:-1]\
                                     +angle_lr[1:,1:-1]+angle_lr[:-1,1:-1])
        # edges
        self.angle[0,1:-1] = (1.0/3.0)*(angle_lr[0,1:-1]+angle_ud[0,1:]+angle_ud[0,:-1])
        self.angle[-1,1:-1] = (1.0/3.0)*(angle_lr[-1,1:-1]+angle_ud[-1,1:]+angle_ud[-1,:-1])
        self.angle[1:-1,0] = (1.0/3.0)*(angle_ud[1:-1,0]+angle_lr[1:,0]+angle_lr[:-1,0])
        self.angle[1:-1,-1] = (1.0/3.0)*(angle_ud[1:-1,-1]+angle_lr[1:,-1]+angle_lr[:-1,-1])
        #conrers
        self.angle[0,0] = 0.5*(angle_lr[0,0]+angle_ud[0,0])
        self.angle[0,-1] = 0.5*(angle_lr[0,-1]+angle_ud[0,-1])
        self.angle[-1,0] = 0.5*(angle_lr[-1,0]+angle_ud[-1,0])
        self.angle[-1,-1] = 0.5*(angle_lr[-1,-1]+angle_ud[-1,-1])
        
    def _calculate_angle_rho(self): 
        self.angle_rho = np.arctan2(np.diff(0.5*(self.y_vert[1:,:]+self.y_vert[:-1,:])), \
                                    np.diff(0.5*(self.x_vert[1:,:]+self.x_vert[:-1,:])))        
    
    def calculate_orthogonality(self):
        '''
        Calculate orthogonality error in radians
        '''
        z = self.x_vert + 1j*self.y_vert
        du = np.diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang1 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,:-1]
        ang2 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[:-1,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang3 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        du = np.diff(z, axis=1); du = (du/abs(du))[1:,:]
        dv = np.diff(z, axis=0); dv = (dv/abs(dv))[:,1:]
        ang4 = np.arccos(du.real*dv.real + du.imag*dv.imag)
        ang = np.mean([abs(ang1), abs(ang2), abs(ang3), abs(ang4)], axis=0)
        ang = (ang-np.pi/2.0)
        return ang
    
    def mask_polygon(self, polyverts, mask_value=0.0):
        """
        Mask Cartesian points contained within the polygon defined by polyverts
        
        A cell is masked if the cell center (x_rho, y_rho) is within the
        polygon. Other sub-masks (mask_u, mask_v, and mask_psi) are updated
        automatically.
        
        mask_value [=0.0] may be specified to alter the value of the mask set
        within the polygon.  E.g., mask_value=1 for water points.
        """
        
        polyverts = np.asarray(polyverts)
        assert polyverts.ndim == 2, \
            'polyverts must be a 2D array, or a similar sequence'
        assert polyverts.shape[1] == 2, \
            'polyverts must be two columns of points'
        assert polyverts.shape[0] > 2, \
            'polyverts must contain at least 3 points'
        
        mask = self.mask_rho
        #inside = points_inside_poly(
        #    np.vstack( (self.x_rho.flatten(), self.y_rho.flatten()) ).T,
        #    polyverts)
        path = mpl.path.Path(polyverts)
        inside = path.contains_points(np.vstack( (self.x_rho.flatten(), self.y_rho.flatten()) ).T)
        if np.any(inside):
            self.mask_rho.flat[inside] = mask_value
    
    def _get_mask_u(self):
        return self.mask_rho[:,1:]*self.mask_rho[:,:-1]
    
    def _get_mask_v(self):
        return self.mask_rho[1:,:]*self.mask_rho[:-1,:]
    
    def _get_mask_psi(self):
        return self.mask_rho[1:,1:]*self.mask_rho[:-1,1:]* \
               self.mask_rho[1:,:-1]*self.mask_rho[:-1,:-1]

    def _set_mask_rho(self, mask_rho):
        self.mask_rho = mask_rho
 
    x = property(lambda self: self.x_vert, None, None, 'Return x_vert')
    y = property(lambda self: self.y_vert, None, None, 'Return x_vert')
    mask = property(lambda self: self.mask_rho, _set_mask_rho, None, 'Return mask_rho')
    mask_u   = property(_get_mask_u, None, None, 'Return mask_u')
    mask_v   = property(_get_mask_v, None, None, 'Return mask_v')
    mask_psi = property(_get_mask_psi, None, None, 'Return mask_psi')


class CGrid_geo(CGrid):
    """
    Curvilinear Arakawa C-grid defined in geographic coordinates
    
    For a geographic grid, a projection may be specified, or The default
    projection for will be defined by the matplotlib.toolkits.Basemap
    projection:
    
    proj = Basemap(projection='merc', resolution=None, lat_ts=0.0)
    
    For a geographic grid, the cell widths are determined by the great
    circle distances. Angles, however, are defined using the projected
    coordinates, so a projection that conserves angles must be used. This
    means typically either Mercator (projection='merc') or Lambert
    Conformal Conic (projection='lcc').
    """
    def _calculate_metrics(self):
        # calculate metrics based on x and y grid
        super(CGrid_geo, self)._calculate_metrics()

        # optionally calculate dx and dy based on great circle distances
        # for more accurate cell sizes.
        if self.use_gcdist:
            geod = pyproj.Geod(ellps=self.ellipse)
            az_forward, az_back, dx = geod.inv(self.lon[:,1:],  self.lat[:,1:], \
                                               self.lon[:,:-1], self.lat[:,:-1])
            self.dx = 0.5*(dx[1:,:]+dx[:-1,:])
            self.pm = 1.0/self.dx
            az_forward, az_back, dy = geod.inv(self.lon[1:,:],  self.lat[1:,:], \
                                               self.lon[:-1,:], self.lat[:-1,:])
            self.dy = 0.5*(dy[:,1:]+dy[:,:-1])
            self.pn = 1.0/self.dy


    def _calculate_derivative_metrics(self):
        if isinstance(self.dy, np.ma.MaskedArray):
            self.dndx = np.ma.zeros(self.dy.shape, dtype='d')
        else:
            self.dndx = np.zeros(self.dy.shape, dtype='d')
        
        if isinstance(self.dx, np.ma.MaskedArray):
            self.dmde = np.ma.zeros(self.dx.shape, dtype='d')
        else:
            self.dmde = np.zeros(self.dx.shape, dtype='d')
        
        self.dndx[1:-1,1:-1] = 0.5*(self.dy[1:-1,2:] - self.dy[1:-1,:-2])
        self.dmde[1:-1,1:-1] = 0.5*(self.dx[2:,1:-1] - self.dx[:-2,1:-1])
        
    def _calculate_angle_rho(self):
        if isinstance(self.lon, np.ma.MaskedArray) or \
           isinstance(self.lat, np.ma.MaskedArray):
            self.angle_rho = np.ma.zeros(self.lon.shape, dtype='d')
        else:
            self.angle_rho = np.zeros(self.lon.shape, dtype='d')

        # calculate metrics based on x and y grid
        super(CGrid_geo, self)._calculate_angle_rho()

        # optionally calculate dx and dy based on great circle distances
        # for more accurate cell sizes.
        if self.use_gcdist:
            geod = pyproj.Geod(ellps=self.ellipse)
            az_forward, az_back, dx = geod.inv(self.lon[:,:-1], self.lat[:,:-1], \
                                               self.lon[:,1:], self.lat[:,1:])

            angle = 0.5 * (az_forward[1:,:] + az_forward[:-1,:])
            self.angle_rho = (90 - angle) * np.pi/180.


    def __init__(self, lon_vert, lat_vert, proj, use_gcdist=True, ellipse='WGS84', \
                    lon_rho=None, lat_rho=None, lon_u=None, lat_u=None, \
                    lon_v=None, lat_v=None, lon_psi=None, lat_psi=None, dx=None, dy=None, \
                    dndx=None, dmde=None, angle_rho=None):

        x, y = proj(lon_vert, lat_vert)
        self.lon_vert = lon_vert
        self.lat_vert = lat_vert
        self.proj = proj

        self.use_gcdist = use_gcdist
        self.ellipse = ellipse

        if lon_rho is None or lat_rho is None or lon_u is None or lat_u is None or \
             lon_v is None or lat_v is None or lon_psi is None or lat_psi is None:

            super(CGrid_geo, self).__init__(x, y)

            self.lon_rho, self.lat_rho = self.proj(self.x_rho, self.y_rho,
                                                   inverse=True)
            self.lon_u, self.lat_u = self.proj(self.x_u, self.y_u, inverse=True)
            self.lon_v, self.lat_v = self.proj(self.x_v, self.y_v, inverse=True)
            self.lon_psi, self.lat_psi = self.proj(self.x_psi, self.y_psi,
                                                   inverse=True)
        else:
            self.lon_rho = lon_rho
            self.lat_rho = lat_rho
            self.lon_u = lon_u
            self.lat_u = lat_u
            self.lon_v = lon_v
            self.lat_v = lat_v
            self.lon_psi = lon_psi
            self.lat_psi = lat_psi
            #calculate cartesian position
            self.x_vert, self.y_vert = proj(lon_vert, lat_vert)
            self.x_rho, self.y_rho = proj(lon_rho, lat_rho)
            self.x_u, self.y_u = proj(lon_u, lat_u)
            self.x_v, self.y_v = proj(lon_v, lat_v)
            self.x_psi, self.y_psi = proj(lon_psi, lat_psi)

        if dx is None or dy is None:
            self._calculate_metrics()
        else:
            self.dx = dx
            self.dy = dy

        self.xl = np.maximum(self.dx[0,:].sum(), self.dx[-1,:].sum())
        self.el = np.maximum(self.dy[:,0].sum(), self.dy[:,-1].sum())

        if dndx is None or dmde is None:
            self._calculate_derivative_metrics()
        else:
            self.dndx = dndx
            self.dmde = dmde

        if angle_rho is None:
            self._calculate_angle_rho()
        else:
            self.angle_rho = angle_rho

        self.f = 2.0 * 7.29e-5 * np.sin(self.lat_rho * np.pi / 180.0)
        self.spherical = 'T'


    def mask_polygon_geo(lonlat_verts, mask_value=0.0):
        lon, lat = list(zip(*lonlat_verts))
        x, y = proj(lon, lat, inverse=True)
        self.mask_polygon(list(zip(x, y)), mask_value)
    
    lon = property(lambda self: self.lon_vert, None, None, 'Shorthand for lon_vert')
    lat = property(lambda self: self.lat_vert, None, None, 'Shorthand for lat_vert')
        


class Gridgen(CGrid):
    """
    docstring for Gridgen
    """

    
    def generate_grid(self):
        
        if self._gn is not None:
            self._libgridgen.gridnodes_destroy(self._gn)
        
        nbry = len(self.xbry)
        
        nsigmas = ctypes.c_int(0)
        sigmas = ctypes.c_void_p(0)
        nrect = ctypes.c_int(0)
        xrect =  ctypes.c_void_p(0)
        yrect = ctypes.c_void_p(0)
        
        if self.focus is None:
            ngrid = ctypes.c_int(0)
            xgrid = ctypes.POINTER(ctypes.c_double)()
            ygrid = ctypes.POINTER(ctypes.c_double)()
        else:
            y, x =  np.mgrid[0:1:self.ny*1j, 0:1:self.nx*1j]
            xgrid, ygrid = self.focus(x, y)
            ngrid = ctypes.c_int(xgrid.size)
            xgrid = (ctypes.c_double * xgrid.size)(*xgrid.flatten())
            ygrid = (ctypes.c_double * ygrid.size)(*ygrid.flatten())
        
        self._gn = self._libgridgen.gridgen_generategrid2(
             ctypes.c_int(nbry), 
             (ctypes.c_double * nbry)(*self.xbry), 
             (ctypes.c_double * nbry)(*self.ybry), 
             (ctypes.c_double * nbry)(*self.beta),
             ctypes.c_int(self.ul_idx), 
             ctypes.c_int(self.nx), 
             ctypes.c_int(self.ny), 
             ngrid, 
             xgrid, 
             ygrid,
             ctypes.c_int(self.nnodes), 
             ctypes.c_int(self.newton), 
             ctypes.c_double(self.precision),
             ctypes.c_int(self.checksimplepoly), 
             ctypes.c_int(self.thin), 
             ctypes.c_int(self.nppe),
             ctypes.c_int(self.verbose),
             ctypes.byref(nsigmas), 
             ctypes.byref(sigmas), 
             ctypes.byref(nrect),
             ctypes.byref(xrect), 
             ctypes.byref(yrect) )
        
        x = self._libgridgen.gridnodes_getx(self._gn)
        x = np.asarray([x[0][i] for i in range(self.ny*self.nx)])
        # x = np.asarray([x[j][i] for j in range(self.ny) for i in range(self.nx)])
        x.shape = (self.ny, self.nx)
        
        y = self._libgridgen.gridnodes_gety(self._gn)    
        y = np.asarray([y[0][i] for i in range(self.ny*self.nx)])            
        # y = np.asarray([y[j][i] for j in range(self.ny) for i in range(self.nx)])
        y.shape = (self.ny, self.nx)
        
        if np.any(np.isnan(x)) or np.any(np.isnan(y)):
            x = np.ma.masked_where(np.isnan(x), x)
            y = np.ma.masked_where(np.isnan(y), y)
        
        # if self.proj is not None:
        #     lon, lat = self.proj(x, y, inverse=True)
        #     super(Gridgen, self).__init__(lon, lat, proj=self.proj)
        # else:
        super(Gridgen, self).__init__(x, y)



    def __init__(self, xbry, ybry, beta, shape, ul_idx=0, \
                 focus=None, proj=None, \
                 nnodes=14, precision=1.0e-12, nppe=3, \
                 newton=True, thin=True, checksimplepoly=True, verbose=False):

        #self._libgridgen = np.ctypeslib.load_library('libgridgen',__file__)
        self._libgridgen = np.ctypeslib.load_library('libgridgen', pyroms.__path__[0])

        # In MacOSX, use of c_void_p does not return proper structure.
        # (An integer address is returned and subsequent use results in a 
        #   Segmentation Fault)
        # All structures have to be declared.  

        # NODETYPE is enum(erated)
        (NT_NONE, NT_DD, NT_CEN, NT_COR) = (0, 1, 2, 3)
        
        class GRIDSTATS(ctypes.Structure):
          _fields_ = [
            ("mdo", ctypes.c_double),
            ("imdo", ctypes.c_int),
            ("jmdo", ctypes.c_int),
            ("ado", ctypes.c_double),
            ("mar", ctypes.c_double),
            ("aar", ctypes.c_double)
          ]
        
        class GRIDNODES(ctypes.Structure):
          _fields_ = [
            ("nx", ctypes.c_int),
            ("ny", ctypes.c_int),
            ("gx", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
            ("gy", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
            ("type", ctypes.c_int),
            ("validated", ctypes.c_int),
            ("stats", ctypes.POINTER(GRIDSTATS)),
            ("nextpoint", ctypes.c_int)
          ]
        
        class EXTENT(ctypes.Structure):
          _fields_ = [
            ("xmin", ctypes.c_double),
            ("xmax", ctypes.c_double),
            ("ymin", ctypes.c_double),
            ("ymax", ctypes.c_double)
          ]
        
        class POLY(ctypes.Structure):
          _fields_ = [
            ("n", ctypes.c_int),
            ("nallocated", ctypes.c_int),
            ("e", EXTENT),
            ("x", ctypes.POINTER(ctypes.c_double)),
            ("y", ctypes.POINTER(ctypes.c_double))
          ]
        
        # SUBGRID
        # A forward declaration of this structure is used
        # (1) Defined it first with pass
        # (2) Define the fields next
        
        # SUBGRID (1)
        class SUBGRID(ctypes.Structure):
          pass
        
        class GRIDMAP(ctypes.Structure):
          _fields_ = [
            ("bound", ctypes.POINTER(POLY)),
            ("trunk", ctypes.POINTER(SUBGRID)),
            ("nleaves", ctypes.c_int),
            ("nce1", ctypes.c_int),
            ("nce2", ctypes.c_int),
            ("gx", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
            ("gy", ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),
            ("sign", ctypes.c_int)
          ]
        
        # SUBGRID (2)
        SUBGRID._fields_ = [
            ("gmap", ctypes.POINTER(GRIDMAP)),
            ("bound", ctypes.POINTER(POLY)),
            ("mini", ctypes.c_int),
            ("maxi", ctypes.c_int),
            ("minj", ctypes.c_int),
            ("maxj", ctypes.c_int),
            ("half1", ctypes.POINTER(SUBGRID)),
            ("half2", ctypes.POINTER(SUBGRID))
        ]
                
        #self._libgridgen.gridgen_generategrid2.restype = ctypes.c_void_p
        self._libgridgen.gridgen_generategrid2.restype = ctypes.POINTER(GRIDNODES)
        self._libgridgen.gridnodes_getx.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
        self._libgridgen.gridnodes_gety.restype = ctypes.POINTER(ctypes.POINTER(ctypes.c_double))
        self._libgridgen.gridnodes_getnce1.restype = ctypes.c_int
        self._libgridgen.gridnodes_getnce2.restype = ctypes.c_int
        #self._libgridgen.gridmap_build.restype = ctypes.c_void_p
        self._libgridgen.gridmap_build.restype = ctypes.POINTER(GRIDMAP)
        
        self.xbry = np.asarray(xbry, dtype='d')
        self.ybry = np.asarray(ybry, dtype='d')
        self.beta = np.asarray(beta, dtype='d')
        assert self.beta.sum() == 4.0, 'sum of beta must be 4.0'
        self.shape = shape
        self.ny = shape[0]
        self.nx = shape[1]
        self.ul_idx = ul_idx
        self.focus = focus
        self.nnodes = nnodes
        self.precision = precision
        self.nppe = nppe
        self.newton = newton
        self.thin = thin
        self.checksimplepoly = checksimplepoly
        self.verbose = verbose
        
        self.proj = proj
        if self.proj is not None:
            self.xbry, self.ybry = proj(self.xbry, self.ybry)
        
        self._gn = None
        self.generate_grid()
    
    def __del__(self):
        """delete gridnode object upon deletion"""
        self._libgridgen.gridnodes_destroy(self._gn)


def rho_to_vert(xr, yr, pm, pn, ang):
    Mp, Lp = xr.shape
    x = np.empty((Mp+1, Lp+1), dtype='d')
    y = np.empty((Mp+1, Lp+1), dtype='d')
    x[1:-1, 1:-1] = 0.25*(xr[1:,1:]+xr[1:,:-1]+xr[:-1,1:]+xr[:-1,:-1])
    y[1:-1, 1:-1] = 0.25*(yr[1:,1:]+yr[1:,:-1]+yr[:-1,1:]+yr[:-1,:-1])
    
    # east side
    theta = 0.5*(ang[:-1,-1]+ang[1:,-1])
    dx = 0.5*(1.0/pm[:-1,-1]+1.0/pm[1:,-1])
    dy = 0.5*(1.0/pn[:-1,-1]+1.0/pn[1:,-1])
    x[1:-1,-1] = x[1:-1,-2] + dx*np.cos(theta)
    y[1:-1,-1] = y[1:-1,-2] + dx*np.sin(theta)
    
    # west side
    theta = 0.5*(ang[:-1,0]+ang[1:,0])
    dx = 0.5*(1.0/pm[:-1,0]+1.0/pm[1:,0])
    dy = 0.5*(1.0/pn[:-1,0]+1.0/pn[1:,0])
    x[1:-1,0] = x[1:-1,1] - dx*np.cos(theta)
    y[1:-1,0] = y[1:-1,1] - dx*np.sin(theta)
    
    # north side
    theta = 0.5*(ang[-1,:-1]+ang[-1,1:])
    dx = 0.5*(1.0/pm[-1,:-1]+1.0/pm[-1,1:])
    dy = 0.5*(1.0/pn[-1,:-1]+1.0/pn[-1,1:])
    x[-1,1:-1] = x[-2,1:-1] - dy*np.sin(theta)
    y[-1,1:-1] = y[-2,1:-1] + dy*np.cos(theta)
    
    # here we are now going to the south side..
    theta = 0.5*(ang[0,:-1]+ang[0,1:])
    dx = 0.5*(1.0/pm[0,:-1]+1.0/pm[0,1:])
    dy = 0.5*(1.0/pn[0,:-1]+1.0/pn[0,1:])
    x[0,1:-1] = x[1,1:-1] + dy*np.sin(theta)
    y[0,1:-1] = y[1,1:-1] - dy*np.cos(theta)
    
    #Corners
    x[0,0] = 4.0*xr[0,0]-x[1,0]-x[0,1]-x[1,1]
    x[-1,0] = 4.0*xr[-1,0]-x[-2,0]-x[-1,1]-x[-2,1]
    x[0,-1] = 4.0*xr[0,-1]-x[0,-2]-x[1,-1]-x[1,-2]
    x[-1,-1] = 4.0*xr[-1,-1]-x[-2,-2]-x[-2,-1]-x[-1,-2]
    
    y[0,0] = 4.0*yr[0,0]-y[1,0]-y[0,1]-y[1,1]
    y[-1,0] = 4.0*yr[-1,0]-y[-2,0]-y[-1,1]-y[-2,1]
    y[0,-1] = 4.0*yr[0,-1]-y[0,-2]-y[1,-1]-y[1,-2]
    y[-1,-1] = 4.0*yr[-1,-1]-y[-2,-2]-y[-2,-1]-y[-1,-2]
    
    return x, y


def rho_to_vert_geo(lonr, latr, lonp, latp):
    Mm, Lm = lonr.shape
    lon = np.zeros((Mm+1,Lm+1))
    lat = np.zeros((Mm+1,Lm+1))

    lon[1:-1, 1:-1] = lonp[:,:]
    lat[1:-1, 1:-1] = latp[:,:]

    #North edge
    lon[Mm,0:-2] = lonr[Mm-1,0:-1] - ( lonp[Mm-2,:] - lonr[Mm-1,0:-1] )
    lon[Mm,-2:] = lonr[Mm-1,-2:] - ( lonp[Mm-2,-2:] - lonr[Mm-1,-2:] )
    lat[Mm,0:-2] = latr[Mm-1,0:-1] - ( latp[Mm-2,:] - latr[Mm-1,0:-1] )
    lat[Mm,-2:] = latr[Mm-1,-2:] - ( latp[Mm-2,-2:] - latr[Mm-1,-2:] )
        
    #South edge
    lon[0,0:-2] = lonr[0,0:-1] - ( lonp[0,:] - lonr[0,0:-1] )
    lon[0,-2:] = lonr[0,-2:] - ( lonp[0,-2:] - lonr[0,-2:] )
    lat[0,0:-2] = latr[0,0:-1] - ( latp[0,:] - latr[0,0:-1] )
    lat[0,-2:] = latr[0,-2:] - ( latp[0,-2:] - latr[0,-2:] )

    #East edge
    lon[0:-2,Lm] = lonr[0:-1,Lm-1] - ( lonp[:,Lm-2] - lonr[0:-1,Lm-1] )
    lon[-2:,Lm] = lonr[-2:,Lm-1] - ( lonp[-2:,Lm-2] - lonr[-2:,Lm-1] )
    lat[0:-2,Lm] = latr[0:-1,Lm-1] - ( latp[:,Lm-2] - latr[0:-1,Lm-1] )
    lat[-2:,Lm] = latr[-2:,Lm-1] - ( latp[-2:,Lm-2] - latr[-2:,Lm-1] )

    #West edge
    lon[0:-2,0] = lonr[0:-1,0] - ( lonp[:,0] - lonr[0:-1,0] )
    lon[-2:,0] = lonr[-2:,0] - ( lonp[-2:,0] - lonr[-2:,0] )
    lat[0:-2,0] = latr[0:-1,0] - ( latp[:,0] - latr[0:-1,0] )
    lat[-2:,0] = latr[-2:,0] - ( latp[-2:,0] - latr[-2:,0] )

    return lon, lat


class edit_mask_mesh(object):
    """
    Interactive mask editor

    edit_mask_mesh(grd, proj)

    Edit grd mask. Mask/Unsmask cell by a simple click on the cell.
    Mask modification are store in mask_change.txt for further use.

    Key commands:
        e : toggle between Editing/Viewing mode
    """
    
    def _on_key(self, event):
        if event.key == 'e':
            self._clicking = not self._clicking
            plt.title('Editing %s -- click "e" to toggle' % self._clicking)
            plt.draw()
    
    def _on_click(self, event):
        x, y = event.xdata, event.ydata
        if event.button==1 and event.inaxes is not None and self._clicking == True:
            d = (x-self._xc)**2 + (y-self._yc)**2
            if isinstance(self.xv, np.ma.MaskedArray):
                idx = np.argwhere(d[~self._xc.mask] == d.min())
            else:
                idx = np.argwhere(d.flatten() == d.min())
            self._mask[idx] = float(not self._mask[idx])
            j, i = np.argwhere(d == d.min())[0]
            self.mask[j, i] = float(not self.mask[j, i])
            #open output file
            f = open('mask_change.txt','a')
            #value = (i, j, self.mask[i, j])
            #s = str(value)
            s = '%d %d %f' %(i, j, self.mask[j, i])
            f.write(s + '\n')
            #close file
            f.close()
            self._pc.set_array(self._mask)
            self._pc.changed()
            plt.draw()
    
    def __init__(self, grd, proj=None, **kwargs):

        if type(grd).__name__ == 'ROMS_Grid':
            try:
                xv = grd.hgrid.lon_vert
                yv = grd.hgrid.lat_vert
                mask = grd.hgrid.mask_rho
            except:
                xv = grd.hgrid.x_vert
                yv = grd.hgrid.y_vert
                mask = grd.hgrid.mask_rho

        if type(grd).__name__ == 'CGrid_geo':
            try:
                xv = grd.lon_vert
                yv = grd.lat_vert
                mask = grd.mask_rho
            except:
                xv = grd.x_vert
                yv = grd.y_vert
                mask = grd.mask_rho

        assert xv.shape == yv.shape, 'xv and yv must have the same shape'
        for dx, dq in zip(xv.shape, mask.shape):
             assert dx==dq+1, \
             '''xv and yv must be cell verticies
             (i.e., one cell bigger in each dimension)'''
        
        self.xv = xv
        self.yv = yv
        
        self.mask = mask

        self.proj = proj
        
        land_color = kwargs.pop('land_color', (0.6, 1.0, 0.6))
        sea_color = kwargs.pop('sea_color', (0.6, 0.6, 1.0))

        cm = plt.matplotlib.colors.ListedColormap([land_color, sea_color], 
                                                 name='land/sea')

        if self.proj is None:
            self._pc = plt.pcolor(xv, yv, mask, cmap=cm, vmin=0, vmax=1, edgecolor='k', **kwargs)
        else:
            xv, yv = self.proj(xv, yv)
            self._pc = Basemap.pcolor(self.proj, xv, yv, mask, cmap=cm, vmin=0, vmax=1, edgecolor='k', **kwargs)
            self.proj.drawcoastlines()

        self._xc = 0.25*(xv[1:,1:]+xv[1:,:-1]+xv[:-1,1:]+xv[:-1,:-1])
        self._yc = 0.25*(yv[1:,1:]+yv[1:,:-1]+yv[:-1,1:]+yv[:-1,:-1])
        
        if isinstance(self.xv, np.ma.MaskedArray):
            self._mask = mask[~self._xc.mask]
        else:
            self._mask = mask.flatten()
        
        plt.connect('button_press_event', self._on_click)
        plt.connect('key_press_event', self._on_key)
        self._clicking = False
        plt.title('Editing %s -- click "e" to toggle' % self._clicking)
        plt.draw()



class edit_mask_mesh_ij(object):
    """
    Interactive mask editor

    edit_mask_mesh_ij(grd)

    Edit grd mask. Mask/Unsmask cell by a simple click on the cell.
    Mask modification are store in mask_change.txt for further use.

    Key commands:
        e : toggle between Editing/Viewing mode
    """
    
    def _on_key(self, event):
        if event.key == 'e':
            self._clicking = not self._clicking
            plt.title('Editing %s -- click "e" to toggle' % self._clicking)
            plt.draw()
    
    def _on_click(self, event):
        x, y = event.xdata, event.ydata
        if event.button==1 and event.inaxes is not None and self._clicking == True:
            d = (x-self._xc)**2 + (y-self._yc)**2
            j, i = np.argwhere(d == d.min())[0]
            self.mask[j, i] = float(not self.mask[j, i])
            #open output file
            f = open('mask_change.txt','a')
            #value = (i, j, self.mask[i, j])
            #s = str(value)
            s = '%d %d %f' %(i, j, self.mask[j, i])
            f.write(s + '\n')
            #close file
            f.close()
            self._pc.set_array(self.mask)
            self._pc.changed()
            plt.draw()
    
    def __init__(self, grd, coast=None, **kwargs):

        if type(grd).__name__ == 'ROMS_Grid':
            try:
                x = list(range(grd.hgrid.lon_vert.shape[1]))
                y = list(range(grd.hgrid.lat_vert.shape[0]))
                xv, yv = np.meshgrid(x,y)
                mask = grd.hgrid.mask_rho
            except:
                x = list(range(grd.hgrid.x_vert.shape[1]))
                y = list(range(grd.hgrid.y_vert.shape[0]))
                xv, yv = np.meshgrid(x,y)
                mask = grd.hgrid.mask_rho

        if type(grd).__name__ == 'CGrid_geo':
            try:
                x = list(range(grd.lon_vert.shape[1]))
                y = list(range(grd.lat_vert.shape[0]))
                xv, yv = np.meshgrid(x,y)
                mask = grd.mask_rho
            except:
                x = list(range(grd.x_vert.shape[1]))
                y = list(range(grd.y_vert.shape[0]))
                xv, yv = np.meshgrid(x,y)
                mask = grd.mask_rho

        assert xv.shape == yv.shape, 'xv and yv must have the same shape'
        for dx, dq in zip(xv.shape, mask.shape):
             assert dx==dq+1, \
             '''xv and yv must be cell verticies
             (i.e., one cell bigger in each dimension)'''
        
        self.xv = xv
        self.yv = yv
        
        self.mask = mask

        if coast is None:
            #get the coast from Basemap (GEOS Library)
            map = pyroms.utility.get_grid_proj(grd, resolution='h', area_thresh=5)
            coast = pyroms.utility.get_coast_from_map(map)

        self.coast = coast

        ijcoast = pyroms.utility.ijcoast(self.coast, grd)

        self.ijcoast = ijcoast + 0.5 

        land_color = kwargs.pop('land_color', (0.6, 1.0, 0.6))
        sea_color = kwargs.pop('sea_color', (0.6, 0.6, 1.0))

        cm = plt.matplotlib.colors.ListedColormap([land_color, sea_color], 
                                                 name='land/sea')

        self._pc = plt.imshow(mask, cmap=cm, origin='lower', \
          interpolation='nearest', extent=(0,xv.shape[1]-1,0,xv.shape[0]-1), \
          **kwargs)
        plt.plot(xv, yv, color='k', ls='-', lw=.5)
        plt.plot(xv.T, yv.T, color='k', ls='-', lw=.5)
        plt.plot(self.ijcoast[:,0], self.ijcoast[:,1], color='k', ls='-', lw=1.25)

        self._xc = 0.25*(xv[1:,1:]+xv[1:,:-1]+xv[:-1,1:]+xv[:-1,:-1])
        self._yc = 0.25*(yv[1:,1:]+yv[1:,:-1]+yv[:-1,1:]+yv[:-1,:-1])
        
        plt.connect('button_press_event', self._on_click)
        plt.connect('key_press_event', self._on_key)
        self._clicking = False
        plt.title('Editing %s -- click "e" to toggle' % self._clicking)
        plt.draw()



def uvp_masks(rmask):
    '''
    return u-, v-, and psi-masks based on input rho-mask
    
    Parameters
    ----------
    
    rmask : ndarray
        mask at CGrid rho-points
    
    Returns
    -------
    (umask, vmask, pmask) : ndarrays
        masks at u-, v-, and psi-points
    '''
    rmask = np.asarray(rmask)
    assert rmask.ndim == 2, 'rmask must be a 2D array'
    assert np.all((rmask==0)|(rmask==1)), 'rmask array must contain only ones and zeros.'

    umask = rmask[:, :-1] * rmask[:, 1:]
    vmask = rmask[:-1, :] * rmask[1:, :]
    pmask = rmask[:-1, :-1] * rmask[:-1, 1:] * rmask[1:, :-1] * rmask[1:, 1:]

    return umask, vmask, pmask



class get_position_from_map(object):
    """
    Get cell index position Interactively

    get_position_from_map(grd, proj)

    Get index i, j as well as lon, lat coordinates for one cell
    simply by clicking on the cell.

    Key commands:
        i : toggle between Interactive/Viewing mode
    """
    def _on_key(self, event):
        if event.key == 'i':
            self._clicking = not self._clicking
            plt.title('Interactive %s -- click "i" to toggle' % self._clicking)
            plt.draw()

    def _on_click(self, event):
        x, y = event.xdata, event.ydata
        if event.button==1 and event.inaxes is not None and self._clicking == True:
            d = (x-self._xc)**2 + (y-self._yc)**2
            if isinstance(self.xv, np.ma.MaskedArray):
                idx = np.argwhere(d[~self._xc.mask] == d.min())
            else:
                idx = np.argwhere(d.flatten() == d.min())
            j, i = np.argwhere(d == d.min())[0]
            print('Position on the grid (rho point): i =', i, ', j =', j)
            if self.proj is not None:
                lon, lat = self.proj(self._xc[j,i], self._yc[j,i], inverse=True)
                print('corresponding geographical position : lon = ', lon, ', lat =', lat)
            else:
                print('corresponding cartesian position : x = ', self._xc[j,i], ', y =', self._yc[j,i])
    
    def __init__(self, grd, proj=None, **kwargs):

        try:
            xv = grd.hgrid.lon_vert
            yv = grd.hgrid.lat_vert
            mask = grd.hgrid.mask_rho
        except:
            xv = grd.hgrid.x_vert
            yv = grd.hgrid.y_vert
            mask = grd.hgrid.mask_rho
          
        assert xv.shape == yv.shape, 'xv and yv must have the same shape'
        for dx, dq in zip(xv.shape, mask.shape):
             assert dx==dq+1, \
             '''xv and yv must be cell verticies
             (i.e., one cell bigger in each dimension)'''
        
        self.xv = xv
        self.yv = yv
        
        self.mask = mask

        self.proj = proj
        
        land_color = kwargs.pop('land_color', (0.6, 1.0, 0.6))
        sea_color = kwargs.pop('sea_color', (0.6, 0.6, 1.0))

        cm = plt.matplotlib.colors.ListedColormap([land_color, sea_color], 
                                                 name='land/sea')

        if self.proj is None:
            self._pc = plt.pcolor(xv, yv, mask, cmap=cm, vmin=0, vmax=1, edgecolor='k', **kwargs)

        else:
            xv, yv = self.proj(xv, yv)
            self._pc = Basemap.pcolor(self.proj, xv, yv, mask, cmap=cm, vmin=0, vmax=1, edgecolor='k', **kwargs)
            self.proj.drawcoastlines()

        self._xc = 0.25*(xv[1:,1:]+xv[1:,:-1]+xv[:-1,1:]+xv[:-1,:-1])
        self._yc = 0.25*(yv[1:,1:]+yv[1:,:-1]+yv[:-1,1:]+yv[:-1,:-1])
        
        plt.connect('button_press_event', self._on_click)
        plt.connect('key_press_event', self._on_key)
        self._clicking = False
        plt.title('Interactive %s -- click "i" to toggle' % self._clicking)
        plt.draw()
    

if __name__ == '__main__':
    geographic = False
    if geographic:
        from mpl_toolkits.basemap import Basemap
        proj = Basemap(projection='lcc',
                       resolution='i',
                       llcrnrlon=-72.0,
                       llcrnrlat= 40.0,
                       urcrnrlon=-63.0,
                       urcrnrlat=47.0,
                       lat_0=43.0,
                       lon_0=-62.5)
        
        lon = (-71.977385177601761, -70.19173825913137,
               -63.045075098584945,-64.70104074097425)
        lat = (42.88215610827428, 41.056141745853786,
               44.456701607935841, 46.271758064353897)
        beta = [1.0, 1.0, 1.0, 1.0]

        grd = Gridgen(lon, lat, beta, (32, 32), proj=proj)
        
        for seg in proj.coastsegs:
            grd.mask_polygon(seg)
        
        plt.pcolor(grd.x, grd.y, grd.mask)
        plt.show()
    else:
        x = [0.2, 0.85, 0.9, 0.82, 0.23]
        y = [0.2, 0.25, 0.5, 0.82, .83]
        beta = [1.0, 1.0, 0.0, 1.0, 1.0]
        
        grd = Gridgen(x, y, beta, (32, 32))
        
        ax = plt.subplot(111)
        BoundaryInteractor(x, y, beta)
        plt.show()
