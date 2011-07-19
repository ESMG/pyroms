#!/usr/bin/env python
# encoding: utf-8
"""
geo_anim.py

Created by Rob Hetland on 2008-01-14.
Copyright (c) 2008 Texas A&M Univsersity. All rights reserved.
"""

import matplotlib
matplotlib.use('Agg')

from numpy import *
from matplotlib.pyplot import *
import pylab
import zipfile
import octant
import os

kml_preamble = '''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
<Document>
  <name>Time Animation Test</name>
  <open>1</open>
  <description>
    by Rob Hetland
  </description>
  <Folder>
    <name>Frames</name>
'''

kml_frame = '''  <GroundOverlay>
      <TimeSpan>
        <begin>__TIMEBEGIN__</begin>
        <end>__TIMEEND__</end>
      </TimeSpan>
      <color>__COLOR__</color>
      <Icon>
        <href>__FRAME__</href>
      </Icon>
      <LatLonBox>
        <north>__NORTH__</north>
        <south>__SOUTH__</south>
        <east> __EAST__</east>
        <west> __WEST__</west>
      </LatLonBox>
  </GroundOverlay>
'''

kml_legend = '''<ScreenOverlay>
    <name>Legend</name>
    <Icon>
        <href>legend.png</href>
    </Icon>
    <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
    <screenXY x="0.015" y="0.075" xunits="fraction" yunits="fraction"/>
    <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
    <size x="0" y="0" xunits="pixels" yunits="pixels"/>
</ScreenOverlay>
'''

kml_closing = '''  </Folder>
</Document>
</kml>
'''


def kmz_anim(lon, lat, time, prop, **kwargs):
        lon = asarray(lon)
        lat = asarray(lat)
        
        jd = pylab.date2num(time)
        jd_edges = hstack((1.5*jd[0]-0.5*jd[1], 
                           0.5*(jd[1:]+jd[:-1]), 
                           1.5*jd[-1]-0.5*jd[-2]))
        time_edges = pylab.num2date(jd_edges)
        time_starts = time_edges[:-1]
        time_stops = time_edges[1:]
        
        name = kwargs.pop('name', 'overlay')
        color = kwargs.pop('color', '9effffff')
        visibility = str( kwargs.pop('visibility', 1) )
        kmzfile = kwargs.pop('kmzfile', 'overlay.kmz')
        pixels = kwargs.pop('pixels', 300)  # pixels of the max. dimension
        units = kwargs.pop('units', '')
        vmax = kwargs.pop('vmax', prop.max())
        kwargs['vmax'] = vmax
        vmin = kwargs.pop('vmin', prop.min())
        kwargs['vmin'] = vmin
        
        geo_aspect = cos(lat.mean()*pi/180.0)
        xsize = lon.ptp()*geo_aspect
        ysize = lat.ptp()
        
        aspect = ysize/xsize
        if aspect > 1.0:
            figsize = (10.0/aspect, 10.0)
        else:
            figsize = (10.0, 10.0*aspect)
        
        kml_text = kml_preamble
        
        ioff()
        fig = figure(figsize=figsize, dpi=pixels//10, facecolor=None, frameon=False)
        ax = fig.add_axes([0, 0, 1, 1])
        
        f = zipfile.ZipFile(kmzfile, 'w')
        
        for frame in range(prop.shape[0]):
            tstart = time_starts[frame]
            tstop = time_stops[frame]
            print 'Writing frame ', frame, tstart.isoformat(), tstop.isoformat()
            ax.cla()
            pc = ax.pcolor(lon, lat, prop[frame], **kwargs)
            ax.set_xlim(lon.min(), lon.max())
            ax.set_ylim(lat.min(), lat.max())
            ax.set_axis_off()
            icon = 'overlay_%d.png' % frame
            savefig(icon)
            kml_text += kml_frame.replace('__NAME__', name)\
                                 .replace('__COLOR__', color)\
                                 .replace('__VISIBILITY__', visibility)\
                                 .replace('__SOUTH__', str(lat.min()))\
                                 .replace('__NORTH__', str(lat.max()))\
                                 .replace('__EAST__', str(lon.max()))\
                                 .replace('__WEST__', str(lon.min()))\
                                 .replace('__FRAME__', icon)\
                                 .replace('__TIMEBEGIN__', tstart.isoformat())\
                                 .replace('__TIMEEND__', tstop.isoformat())
            
            f.write(icon)
            os.remove(icon)
        
        # legend
        fig = figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
        cax = fig.add_axes([0.0, 0.05, 0.2, 0.90])
        cb = colorbar(pc, cax=cax)
        cb.set_label(units, color='0.9')
        for lab in cb.ax.get_yticklabels():
            setp(lab, 'color', '0.9')
        
        savefig('legend.png')
        f.write('legend.png')
        os.remove('legend.png')
        
        kml_text += kml_legend
        
        kml_text += kml_closing
        f.writestr('overlay.kml', kml_text)
        f.close()


if __name__ == '__main__':
    ncll = octant.io.Dataset('/Users/rob/Archive/GWB/bodden/latlon.nc')
    nc = octant.io.Dataset('/Users/rob/Archive/GWB/bodden/bsh_elev_2001-10.nc')
    
    lat = ncll.variables['lat'][:]
    lon = ncll.variables['lon'][:]
    
    lon, lat = meshgrid(lon, lat)
    
    time = octant.ocean_time(nc, name='time')[:200:4]
    
    propname = 'elev'
    
    prop = nc.variables[propname][:200:4]
    mask = prop == nc.variables[propname].missing_value
    prop = ma.masked_where(mask, prop)
    
    kmz_anim(lon, lat, time.dates, prop, kmzfile='bsh_anim.kmz', 
             name='BSH model -- sea surface height', units='sea surface height [m]')





kml_groundoverlay = '''<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://earth.google.com/kml/2.0">
<Document>
<GroundOverlay>
  <name>__NAME__</name>
  <color>__COLOR__</color>
  <visibility>__VISIBILITY__</visibility>
  <Icon>
    <href>overlay.png</href>
  </Icon>
  <LatLonBox>
    <south>__SOUTH__</south>
    <north>__NORTH__</north>
    <west>__WEST__</west>
    <east>__EAST__</east>
  </LatLonBox>
</GroundOverlay>
<ScreenOverlay>
    <name>Legend</name>
    <Icon>
        <href>legend.png</href>
    </Icon>
    <overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/>
    <screenXY x="0.015" y="0.075" xunits="fraction" yunits="fraction"/>
    <rotationXY x="0.5" y="0.5" xunits="fraction" yunits="fraction"/>
    <size x="0" y="0" xunits="pixels" yunits="pixels"/>
</ScreenOverlay>
</Document>
</kml>
'''

def geo_pcolor(lon, lat, prop, **kwargs):
    """docstring for geo_pcolor"""
    
    name = kwargs.pop('name', 'overlay')
    color = kwargs.pop('color', '9effffff')
    visibility = str( kwargs.pop('visibility', 1) )
    kmzfile = kwargs.pop('kmzfile', 'overlay.kmz')
    pixels = kwargs.pop('pixels', 1024)  # pixels of the max. dimension
    units = kwargs.pop('units', '')
    vmax = kwargs.pop('vmax', prop.max())
    kwargs['vmax'] = vmax
    vmin = kwargs.pop('vmin', prop.min())
    kwargs['vmin'] = vmin
    
    geo_aspect = cos(lat.mean()*pi/180.0)
    xsize = lon.ptp()*geo_aspect
    ysize = lat.ptp()
    
    aspect = ysize/xsize
    if aspect > 1.0:
        figsize = (10.0/aspect, 10.0)
    else:
        figsize = (10.0, 10.0*aspect)
    
    ioff()
    fig = figure(figsize=figsize, facecolor=None, frameon=False, dpi=pixels//10)
    ax = fig.add_axes([0, 0, 1, 1])
    pc = ax.pcolor(lon, lat, prop, **kwargs)
    ax.set_xlim(lon.min(), lon.max())
    ax.set_ylim(lat.min(), lat.max())
    ax.set_axis_off()
    savefig('overlay.png')
    
    f = zipfile.ZipFile(kmzfile, 'w')
    f.writestr('overlay.kml', kml_groundoverlay.replace('__NAME__', name)\
                                               .replace('__COLOR__', color)\
                                               .replace('__VISIBILITY__', visibility)\
                                               .replace('__SOUTH__', str(lat.min()))\
                                               .replace('__NORTH__', str(lat.max()))\
                                               .replace('__EAST__', str(lon.max()))\
                                               .replace('__WEST__', str(lon.min())))
    f.write('overlay.png')
    os.remove('overlay.png')
    
    fig = figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    cb = colorbar(pc, cax=ax)
    cb.set_label(units, color='0.9')
    for lab in cb.ax.get_yticklabels():
        setp(lab, 'color', '0.9')
    
    savefig('legend.png')
    f.write('legend.png')
    os.remove('legend.png')    
    f.close()


if __name__ == '__main__':
    ncll = pyroms.Dataset('/Users/rob/Archive/GWB/bodden/latlon.nc')
    nc = pyroms.Dataset('/Users/rob/Archive/GWB/bodden/bsh_elev_2001-10.nc')
    
    lat = ncll.variables['lat'][:]
    lon = ncll.variables['lon'][:]
    
    lon, lat = meshgrid(lon, lat)
    
    propname = 'elev'
    
    prop = nc.variables[propname][-1]
    mask = prop == nc.variables[propname].missing_value
    prop = ma.masked_where(mask, prop)
    
    geo_pcolor(lon, lat, prop, kmzfile='bsh.kmz', \
               name='BSH model -- sea surface height',\
               units='sea surface height [m]')

