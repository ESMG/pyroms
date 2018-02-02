# encoding: utf-8

import sys
import os
import numpy as np
from mpl_toolkits.basemap import Basemap
from datetime import datetime
try:
  import netCDF4 as netCDF
except:
  import netCDF3 as netCDF

import pyroms
from pyroms.hgrid import *
from pyroms.vgrid import *
from pyroms.grid import *
from pyroms import io

#define a dictionary that will remember gridid's that are defined from
#a history and grid file. Because this is defined in this model's name
#space, it will remain persistent.  The keys are the gridid, and the
#values are ROMS_gridinfo objects.
gridid_dictionary={}

class ROMS_Grid(object):
    """
    grd = ROMS_Grid(hgrid, vgrid)

    ROMS Grid object combining horizontal and vertical grid
    """

    def __init__(self, name, hgrid=CGrid, vgrid=s_coordinate):
        self.name = name
        self.hgrid = hgrid
        self.vgrid = vgrid


class ROMS_gridinfo(object):
    '''
    gridinfo = ROMS_gridinfo(gridid,grid_file=None,hist_file=None)

    Return an object with grid information for gridid.

    There are two ways to define the grid information.  If grid_file
    and hist_file are not passed to the object when it is created, the
    information is retrieved from gridid.txt.
    To add new grid please edit your gridid.txt. You need to define
    an environment variable PYROMS_GRIDID_FILE pointing to your
    gridid.txt file. Just copy an existing grid and modify the
    definition accordingly to your case (Be carefull with
    space and blank line).

    If grid_file is the path to a ROMS grid file, and hist_file is the
    path to a ROMS history file, then the grid information will be
    read from those files.  Gridid can then be used to refer to this
    grid information so that the grid and history files do not be
    included in subsequent calls.
    '''

    def __init__(self, gridid,grid_file=None,hist_file=None):
      #first determine if the information for the gridid has already been obtained.
      if gridid in gridid_dictionary:
        #print 'CJMP> gridid found in gridid_dictionary, grid retrieved from dictionary'
        saved_self=gridid_dictionary[gridid]
        for attrib in list(saved_self.__dict__.keys()):
          setattr(self,attrib,getattr(saved_self,attrib))
      else:
        #nope, we need to get the information from gridid.txt or from
        #the grid and history files from the model
        self.id = gridid
        self._get_grid_info(grid_file,hist_file)

        #now save the data in the dictionary, so we don't need to get it again
        gridid_dictionary[gridid]=self

    def _get_grid_info(self,grid_file,hist_file):

      #check if the grid_file and hist_files are both null; if so get data from gridid.txt
      if (type(grid_file)==type(None))&(type(hist_file)==type(None)):
        #print 'CJMP> gridid not in dictionary, data will be retrieved from gridid.txt'
        gridid_file =  os.getenv("PYROMS_GRIDID_FILE")
        data = open(gridid_file,'r')
        lines = data.readlines()
        data.close()

        line_nb = 0
        info = []
        for line in lines:
            s = line.split()
            if s[0] == 'id':
                if s[2] == self.id:
                    for l in range(line_nb, line_nb+5):
                        s = lines[l].split()
                        info.append(s[2])
                        line_nb = line_nb + 1
                    if info[4] == 'roms':
                        for l in range(line_nb, line_nb+4):
                            s = lines[l].split()
                            info.append(s[2])
                    if info[4] == 'z':
                        s = lines[line_nb].split()
                        info.append(s[3:-1])
                        while s[-1:] == ['\\']:
                            line_nb = line_nb + 1
                            s = lines[line_nb].split()
                            info.append(s[:-1])
            line_nb = line_nb + 1

        if info == []:
            raise ValueError('Unknown gridid. Please check your gridid.txt file')

        if info[4] == 'roms':
            self.name     =          info[1]
            self.grdfile  =          info[2]
            self.N        =   np.int(info[3])
            self.grdtype  =          info[4]
            self.Vtrans   =   np.int(info[5])
            self.theta_s  = np.float(info[6])
            self.theta_b  = np.float(info[7])
            self.Tcline   = np.float(info[8])

        elif info[4] == 'z':
            nline = len(info)
            dep = info[5]
            for line in range(6,nline):
                dep = dep + info[line]
            dep = np.array(dep, dtype=np.float)

            self.name    =        info[1]
            self.grdfile =        info[2]
            self.N       = np.int(info[3])
            self.grdtype =        info[4]
            self.depth   = dep

        else:
            raise ValueError('Unknown grid type. Please check your gridid.txt file')

      else: #lets get the grid information from the history and grid files
        #print 'CJMP> getting grid info from ROMS history and grid files'
        assert type(grid_file)!=type(None), 'if specify history file you must specify grid file'
        assert type(hist_file)!=type(None), 'if specify grid file you must specify history file'

        #open history file and get necessary grid information from it.
        hist=netCDF.Dataset(hist_file,'r')

        #put data into ROMS_gridinfo object
        self.name=self.id
        self.grdfile=grid_file
        self.N=len(hist.dimensions['s_rho'])
        self.grdtype='roms'

        #now write this to deal with both ROMS 3 and 2
        try:
          self.Vtrans=np.float(hist.Vstretching)
          self.theta_s=np.float(hist.theta_s)
          self.theta_b=np.float(hist.theta_b)
          self.Tcline=np.float(hist.Tcline)
        except AttributeError:
          try:
            self.Vtrans=np.float(hist.variables['Vstretching'][:])
          except:
            print('variable Vtransform not found in history file. Defaulting to Vtrans=1')
            self.Vtrans=1
          self.theta_s=np.float(hist.variables['theta_s'][:])
          self.theta_b=np.float(hist.variables['theta_b'][:])
          self.Tcline=np.float(hist.variables['Tcline'][:])


def print_ROMS_gridinfo(gridid):
    """
    print_ROMS_gridinfo(gridid)

    return the grid information for gridid
    """

    gridinfo = ROMS_gridinfo(gridid)

    print(' ')
    print('grid information for gridid ', gridinfo.id, ':')
    print(' ')
    print('grid name : ', gridinfo.name)
    print('grid file path : ', gridinfo.grdfile)
    print('number of vertical level : ', gridinfo.N)
    print('grid type : ', gridinfo.grdtype)
    if gridinfo.grdtype == 'roms':
        print('theta_s = ', gridinfo.theta_s)
        print('theta_b = ', gridinfo.theta_b)
        print('Tcline  = ', gridinfo.Tcline)
        #print 'hc      = ', gridinfo.hc
    elif gridinfo.grdtype == 'z':
        print('depth = ', gridinfo.depth)


def list_ROMS_gridid():
    """
    list_ROMS_gridid()

    return the list of the defined gridid
    """

    gridid_file = os.getenv("PYROMS_GRIDID_FILE")
    data = open(gridid_file,'r')
    lines = data.readlines()
    data.close()

    gridid_list = []
    for line in lines:
        s = line.split()
        if s[0] == 'id':
            gridid_list.append(s[2])

    print('List of defined gridid : ', gridid_list)


def get_ROMS_hgrid(gridid):
    """
    hgrid = get_ROMS_hgrid(gridid)

    Load ROMS horizontal grid object
    """

    gridinfo = ROMS_gridinfo(gridid)
    grdfile = gridinfo.grdfile

    nc = io.Dataset(grdfile)

    #Check for cartesian or geographical grid
    spherical = nc.variables['spherical'][0]

    #if it is type byte, then convert to string
    try:
      spherical=spherical.decode('utf8')
    except:
      print('Assuming spherical is integer',spherical, type(spherical))

    #Get horizontal grid
    if ((spherical == 0) or (spherical == 'F')):
        #cartesian grid
        print('Load cartesian grid from file')
        if 'x_vert' in list(nc.variables.keys()) and 'y_vert' in list(nc.variables.keys()):
            x_vert = nc.variables['x_vert'][:]
            y_vert = nc.variables['y_vert'][:]
        elif 'x_rho' in list(nc.variables.keys()) and 'y_rho' in list(nc.variables.keys()) \
                 and 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            x_rho = nc.variables['x_rho'][:]
            y_rho = nc.variables['y_rho'][:]
            pm = nc.variables['pm'][:]
            pn = nc.variables['pn'][:]
            try: angle = nc.variables['angle'][:]
            except: angle = np.zeros(x_rho.shape)
            #compute verts from rho point, pm, pn, angle
            x_vert, y_vert = rho_to_vert(x_rho, y_rho, pm, pn, angle)
        else:
            raise ValueError('NetCDF file must contain x_vert and y_vert \
                     or x_rho, y_rho, pm, pn and angle for a cartesian grid')

        if 'x_rho' in list(nc.variables.keys()) and 'y_rho' in list(nc.variables.keys()) and \
             'x_u' in list(nc.variables.keys()) and 'y_u' in list(nc.variables.keys()) and \
             'x_v' in list(nc.variables.keys()) and 'y_v' in list(nc.variables.keys()) and \
             'x_psi' in list(nc.variables.keys()) and 'y_psi' in list(nc.variables.keys()):
            x_rho = nc.variables['x_rho'][:]
            y_rho = nc.variables['y_rho'][:]
            x_u = nc.variables['x_u'][:]
            y_u = nc.variables['y_u'][:]
            x_v = nc.variables['x_v'][:]
            y_v = nc.variables['y_v'][:]
            x_psi = nc.variables['x_psi'][:]
            y_psi = nc.variables['y_psi'][:]
        else:
            x_rho = None
            y_rho = None
            x_u = None
            y_u = None
            x_v = None
            y_v = None
            x_psi = None
            y_psi = None

        if 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            pm = nc.variables['pm'][:]
            dx = 1. / pm
            pn = nc.variables['pn'][:]
            dy = 1. / pn
        else:
            dx = None
            dy = None

        if 'dndx' in list(nc.variables.keys()) and 'dmde' in list(nc.variables.keys()):
            dndx = nc.variables['dndx'][:]
            dmde = nc.variables['dmde'][:]
        else:
            dndx = None
            dmde = None

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get cartesian grid
        hgrd = CGrid(x_vert, y_vert, x_rho=x_rho, y_rho=y_rho, \
                     x_u=x_u, y_u=y_u, x_v=x_v, y_v=y_v, \
                     x_psi=x_psi, y_psi=y_psi, dx=dx, dy=dy, \
                     dndx=dndx, dmde=dmde, angle_rho=angle)

        #load the mask
        try:
            hgrd.mask_rho = np.array(nc.variables['mask_rho'][:])
        except:
            hgrd.mask_rho = np.ones(hgrd.x_rho.shape)

    else:
        #geographical grid
        print('Load geographical grid from file')
        proj = Basemap(projection='merc', resolution=None, lat_0=0, lon_0=0)
        if 'lon_vert' in list(nc.variables.keys()) and 'lat_vert' in list(nc.variables.keys()):
            lon_vert = nc.variables['lon_vert'][:]
            lat_vert = nc.variables['lat_vert'][:]
        elif 'lon_rho' in list(nc.variables.keys()) and 'lat_rho' in list(nc.variables.keys()) \
                and 'lon_psi' in list(nc.variables.keys()) and 'lat_psi' in list(nc.variables.keys()):
            lon_rho = nc.variables['lon_rho'][:]
            lat_rho = nc.variables['lat_rho'][:]
            lon_psi = nc.variables['lon_psi'][:]
            lat_psi = nc.variables['lat_psi'][:]
            #compute verts from rho and psi point
            lon_vert, lat_vert = rho_to_vert_geo(lon_rho, lat_rho, lon_psi, lat_psi)
        else:
            raise ValueError('NetCDF file must contain lon_vert and lat_vert \
                  or lon_rho, lat_rho, lon_psi, lat_psi for a geographical grid')

        if 'lon_rho' in list(nc.variables.keys()) and 'lat_rho' in list(nc.variables.keys()) and \
              'lon_u' in list(nc.variables.keys()) and 'lat_u' in list(nc.variables.keys()) and \
              'lon_v' in list(nc.variables.keys()) and 'lat_v' in list(nc.variables.keys()) and \
              'lon_psi' in list(nc.variables.keys()) and 'lat_psi' in list(nc.variables.keys()):
            lon_rho = nc.variables['lon_rho'][:]
            lat_rho = nc.variables['lat_rho'][:]
            lon_u = nc.variables['lon_u'][:]
            lat_u = nc.variables['lat_u'][:]
            lon_v = nc.variables['lon_v'][:]
            lat_v = nc.variables['lat_v'][:]
            lon_psi = nc.variables['lon_psi'][:]
            lat_psi = nc.variables['lat_psi'][:]
        else:
            lon_rho = None
            lat_rho = None
            lon_u = None
            lat_u = None
            lon_v = None
            lat_v = None
            lon_psi = None
            lat_psi = None

        if 'pm' in list(nc.variables.keys()) and 'pn' in list(nc.variables.keys()):
            pm = nc.variables['pm'][:]
            dx = 1. / pm
            pn = nc.variables['pn'][:]
            dy = 1. / pn
        else:
            dx = None
            dy = None

        if 'dndx' in list(nc.variables.keys()) and 'dmde' in list(nc.variables.keys()):
            dndx = nc.variables['dndx'][:]
            dmde = nc.variables['dmde'][:]
        else:
            dndx = None
            dmde = None

        if 'angle' in list(nc.variables.keys()):
            angle = nc.variables['angle'][:]
        else:
            angle = None

        #Get geographical grid
        hgrd = CGrid_geo(lon_vert, lat_vert, proj, \
                         lon_rho=lon_rho, lat_rho=lat_rho, \
                         lon_u=lon_u, lat_u=lat_u, lon_v=lon_v, lat_v=lat_v, \
                         lon_psi=lon_psi, lat_psi=lat_psi, dx=dx, dy=dy, \
                         dndx=dndx, dmde=dmde, angle_rho=angle)

        #load the mask
        try:
            hgrd.mask_rho = np.array(nc.variables['mask_rho'][:])
        except:
            hgrd.mask_rho = np.ones(hgrd.lat_rho.shape)

    return hgrd


def get_ROMS_vgrid(gridid, zeta=None):
    """
    vgrid = get_ROMS_vgrid(gridid)

    Load ROMS vertical grid object. vgrid is a s_coordinate or
    a z_coordinate object, depending on gridid.grdtype.
    vgrid.z_r and vgrid.z_w (vgrid.z for a z_coordinate object)
    can be indexed in order to retreive the actual depths. The
    free surface time serie zeta can be provided as an optional
    argument. Note that the values of zeta are not calculated
    until z is indexed, so a netCDF variable for zeta may be passed,
    even if the file is large, as only the values that are required
    will be retrieved from the file.
    """

    gridinfo = ROMS_gridinfo(gridid)
    grdfile = gridinfo.grdfile

    nc = io.Dataset(grdfile)

    #Get vertical grid
    try:
        h = nc.variables['h'][:]
    except:
        raise ValueError('NetCDF file must contain the bathymetry h')

    try:
        hraw = nc.variables['hraw'][:]
    except:
        hraw = None

    if gridinfo.grdtype == 'roms':
        Vtrans = gridinfo.Vtrans
        theta_b = gridinfo.theta_b
        theta_s = gridinfo.theta_s
        Tcline = gridinfo.Tcline
        N = gridinfo.N
        if Vtrans == 1:
            vgrid = s_coordinate(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 2:
            vgrid = s_coordinate_2(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 4:
            vgrid = s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        elif Vtrans == 5:
            vgrid = s_coordinate_5(h, theta_b, theta_s, Tcline, N, hraw=hraw, zeta=zeta)
        else:
            raise Warning('Unknown vertical transformation Vtrans')

    elif  gridinfo.grdtype == 'z':
        N = gridinfo.N
        depth = gridinfo.depth
        vgrid = z_coordinate(h, depth, N)

    else:
        raise ValueError('Unknown grid type')

    return vgrid


def get_ROMS_grid(gridid, zeta=None, hist_file=None,grid_file=None):
    """
    grd = get_ROMS_grid(gridid,hist_file=None,grid_file=None)

    Load ROMS grid object.

    gridid is a string with the name of the grid in it.  If hist_file
       and grid_file are not passed into the function, or are set to
       None, then gridid is used to get the grid data from the
       gridid.txt file.

       if hist_file and grid_file are given, and they are the file
       paths to a ROMS history file and grid file respectively, the
       grid information will be extracted from those files, and gridid
       will be used to name that grid for the rest of the python
       session.

    grd.vgrid is a s_coordinate or
    a z_coordinate object, depending on gridid.grdtype.
    grd.vgrid.z_r and grd.vgrid.z_w (grd.vgrid.z for a
    z_coordinate object) can be indexed in order to retreive the
    actual depths. The free surface time serie zeta can be provided
    as an optional argument. Note that the values of zeta are not
    calculated until z is indexed, so a netCDF variable for zeta may
    be passed, even if the file is large, as only the values that
    are required will be retrieved from the file.
    """

    #in this first call to ROMS_gridinfo, we pass in the history file
    #and gridfile info.  If hist_file and grid_file are defined, the
    #grid info will be extracted from those files and will able to be
    #accessed later by gridid
    gridinfo = ROMS_gridinfo(gridid,hist_file=hist_file,grid_file=grid_file)
    name = gridinfo.name

    #we need not pass in hist_file and grid_file here, because the
    #gridinfo file will already have been initialized by the call to
    #ROMS_gridinfo above.
    hgrd = get_ROMS_hgrid(gridid)
    vgrid = get_ROMS_vgrid(gridid, zeta=zeta)

    #Get ROMS grid
    return ROMS_Grid(name, hgrd, vgrid)


def write_ROMS_grid(grd, filename='roms_grd.nc'):
    """
    write_ROMS_grid(grd, filename)

    Write ROMS_CGrid class on a NetCDF file.
    """

    Mm, Lm = grd.hgrid.x_rho.shape


    # Write ROMS grid to file
    nc = netCDF.Dataset(filename, 'w', format='NETCDF3_64BIT')
    nc.Description = 'ROMS grid'
    nc.Author = 'pyroms.grid.write_grd'
    nc.Created = datetime.now().isoformat()
    nc.type = 'ROMS grid file'

    nc.createDimension('xi_rho', Lm)
    nc.createDimension('xi_u', Lm-1)
    nc.createDimension('xi_v', Lm)
    nc.createDimension('xi_psi', Lm-1)

    nc.createDimension('eta_rho', Mm)
    nc.createDimension('eta_u', Mm)
    nc.createDimension('eta_v', Mm-1)
    nc.createDimension('eta_psi', Mm-1)

    nc.createDimension('xi_vert', Lm+1)
    nc.createDimension('eta_vert', Mm+1)

    nc.createDimension('bath', None)

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        N, = grd.vgrid.s_rho.shape
        nc.createDimension('s_rho', N)
        nc.createDimension('s_w', N+1)

    def write_nc_var(var, name, dimensions, long_name=None, units=None):
        nc.createVariable(name, 'f8', dimensions)
        if long_name is not None:
            nc.variables[name].long_name = long_name
        if units is not None:
            nc.variables[name].units = units
        nc.variables[name][:] = var
        print(' ... wrote ', name)

    if hasattr(grd.vgrid, 's_rho') is True and grd.vgrid.s_rho is not None:
        write_nc_var(grd.vgrid.theta_s, 'theta_s', (), 'S-coordinate surface control parameter')
        write_nc_var(grd.vgrid.theta_b, 'theta_b', (), 'S-coordinate bottom control parameter')
        write_nc_var(grd.vgrid.Tcline, 'Tcline', (), 'S-coordinate surface/bottom layer width', 'meter')
        write_nc_var(grd.vgrid.hc, 'hc', (), 'S-coordinate parameter, critical depth', 'meter')
        write_nc_var(grd.vgrid.s_rho, 's_rho', ('s_rho'), 'S-coordinate at RHO-points')
        write_nc_var(grd.vgrid.s_w, 's_w', ('s_w'), 'S-coordinate at W-points')
        write_nc_var(grd.vgrid.Cs_r, 'Cs_r', ('s_rho'), 'S-coordinate stretching curves at RHO-points')
        write_nc_var(grd.vgrid.Cs_w, 'Cs_w', ('s_w'), 'S-coordinate stretching curves at W-points')

    write_nc_var(grd.vgrid.h, 'h', ('eta_rho', 'xi_rho'), 'bathymetry at RHO-points', 'meter')
    #ensure that we have a bath dependancy for hraw
    if len(grd.vgrid.hraw.shape) == 2:
        hraw = np.zeros((1, grd.vgrid.hraw.shape[0], grd.vgrid.hraw.shape[1]))
        hraw[0,:] = grd.vgrid.hraw
    else:
        hraw = grd.vgrid.hraw
    write_nc_var(hraw, 'hraw', ('bath', 'eta_rho', 'xi_rho'), 'raw bathymetry at RHO-points', 'meter')
    write_nc_var(grd.hgrid.f, 'f', ('eta_rho', 'xi_rho'), 'Coriolis parameter at RHO-points', 'second-1')
    write_nc_var(1./grd.hgrid.dx, 'pm', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in XI', 'meter-1')
    write_nc_var(1./grd.hgrid.dy, 'pn', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in ETA', 'meter-1')
    write_nc_var(grd.hgrid.dmde, 'dmde', ('eta_rho', 'xi_rho'), 'XI derivative of inverse metric factor pn', 'meter')
    write_nc_var(grd.hgrid.dndx, 'dndx', ('eta_rho', 'xi_rho'), 'ETA derivative of inverse metric factor pm', 'meter')
    write_nc_var(grd.hgrid.xl, 'xl', (), 'domain length in the XI-direction', 'meter')
    write_nc_var(grd.hgrid.el, 'el', (), 'domain length in the ETA-direction', 'meter')

    write_nc_var(grd.hgrid.x_rho, 'x_rho', ('eta_rho', 'xi_rho'), 'x location of RHO-points', 'meter')
    write_nc_var(grd.hgrid.y_rho, 'y_rho', ('eta_rho', 'xi_rho'), 'y location of RHO-points', 'meter')
    write_nc_var(grd.hgrid.x_u, 'x_u', ('eta_u', 'xi_u'), 'x location of U-points', 'meter')
    write_nc_var(grd.hgrid.y_u, 'y_u', ('eta_u', 'xi_u'), 'y location of U-points', 'meter')
    write_nc_var(grd.hgrid.x_v, 'x_v', ('eta_v', 'xi_v'), 'x location of V-points', 'meter')
    write_nc_var(grd.hgrid.y_v, 'y_v', ('eta_v', 'xi_v'), 'y location of V-points', 'meter')
    write_nc_var(grd.hgrid.x_psi, 'x_psi', ('eta_psi', 'xi_psi'), 'x location of PSI-points', 'meter')
    write_nc_var(grd.hgrid.y_psi, 'y_psi', ('eta_psi', 'xi_psi'), 'y location of PSI-points', 'meter')
    write_nc_var(grd.hgrid.x_vert, 'x_vert', ('eta_vert', 'xi_vert'), 'x location of cell verticies', 'meter')
    write_nc_var(grd.hgrid.y_vert, 'y_vert', ('eta_vert', 'xi_vert'), 'y location of cell verticies', 'meter')

    if hasattr(grd.hgrid, 'lon_rho'):
        write_nc_var(grd.hgrid.lon_rho, 'lon_rho', ('eta_rho', 'xi_rho'), 'longitude of RHO-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_rho, 'lat_rho', ('eta_rho', 'xi_rho'), 'latitude of RHO-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_u, 'lon_u', ('eta_u', 'xi_u'), 'longitude of U-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_u, 'lat_u', ('eta_u', 'xi_u'), 'latitude of U-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_v, 'lon_v', ('eta_v', 'xi_v'), 'longitude of V-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_v, 'lat_v', ('eta_v', 'xi_v'), 'latitude of V-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_psi, 'lon_psi', ('eta_psi', 'xi_psi'), 'longitude of PSI-points', 'degree_east')
        write_nc_var(grd.hgrid.lat_psi, 'lat_psi', ('eta_psi', 'xi_psi'), 'latitude of PSI-points', 'degree_north')
        write_nc_var(grd.hgrid.lon_vert, 'lon_vert', ('eta_vert', 'xi_vert'), 'longitude of cell verticies', 'degree_east')
        write_nc_var(grd.hgrid.lat_vert, 'lat_vert', ('eta_vert', 'xi_vert'), 'latitude of cell verticies', 'degree_north')

    nc.createVariable('spherical', 'c')
    nc.variables['spherical'].long_name = 'Grid type logical switch'
    nc.variables['spherical'][:] = grd.hgrid.spherical
    print(' ... wrote ', 'spherical')

    write_nc_var(grd.hgrid.angle_rho, 'angle', ('eta_rho', 'xi_rho'), 'angle between XI-axis and EAST', 'radians')

    write_nc_var(grd.hgrid.mask_rho, 'mask_rho', ('eta_rho', 'xi_rho'), 'mask on RHO-points')
    write_nc_var(grd.hgrid.mask_u, 'mask_u', ('eta_u', 'xi_u'), 'mask on U-points')
    write_nc_var(grd.hgrid.mask_v, 'mask_v', ('eta_v', 'xi_v'), 'mask on V-points')
    write_nc_var(grd.hgrid.mask_psi, 'mask_psi', ('eta_psi', 'xi_psi'), 'mask on psi-points')

    nc.close()
