# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 :

# Code to compute distances and angles on a sphere.
# Based on code written by Alistair Adcroft and Matthew Harrison of GFDL

import numpy

def angle_through_center(p1, p2):
    """Angle at center of sphere between two points on the surface of the sphere.
    Positions are given as (latitude,longitude) tuples measured in degrees."""
    phi1 = numpy.deg2rad( p1[0] )
    phi2 = numpy.deg2rad( p2[0] )
    dphi_2 = 0.5 * ( phi2 - phi1 )
    dlambda_2 = 0.5 * numpy.deg2rad( p2[1] - p1[1] )
    a = numpy.sin( dphi_2 )**2 + numpy.cos( phi1 ) * numpy.cos( phi2 ) * ( numpy.sin( dlambda_2 )**2 )
    c = 2. * numpy.arctan2( numpy.sqrt(a), numpy.sqrt( 1. - a ) )
    return c

def angle_between(v1, v2, v3):
    """Returns angle v2-v1-v3 i.e betweeen v1-v2 and v1-v3."""
    # vector product between v1 and v2
    px = v1[1] * v2[2] - v1[2] * v2[1]
    py = v1[2] * v2[0] - v1[0] * v2[2]
    pz = v1[0] * v2[1] - v1[1] * v2[0]
    # vector product between v1 and v3
    qx = v1[1] * v3[2] - v1[2] * v3[1]
    qy = v1[2] * v3[0] - v1[0] * v3[2]
    qz = v1[0] * v3[1] - v1[1] * v3[0]

    ddd = (px * px + py * py + pz * pz) * (qx * qx + qy * qy + qz * qz)
    ddd = (px * qx + py * qy + pz * qz) / numpy.sqrt(ddd)
    angle = numpy.arccos( ddd )
    return angle

def quad_area(lat, lon):
    """Returns area of spherical quad (bounded by great arcs)."""
    # x,y,z are 3D coordinates
    d2r = numpy.deg2rad(1.)
    x = numpy.cos(d2r * lat) * numpy.cos(d2r * lon)
    y = numpy.cos(d2r * lat) * numpy.sin(d2r * lon)
    z = numpy.sin(d2r * lat)
    c0 = (x[ :-1, :-1], y[ :-1, :-1], z[ :-1, :-1])
    c1 = (x[ :-1,1:  ], y[ :-1,1:  ], z[ :-1,1:  ])
    c2 = (x[1:  ,1:  ], y[1:  ,1:  ], z[1:  ,1:  ])
    c3 = (x[1:  , :-1], y[1:  , :-1], z[1:  , :-1])
    a0 = angle_between(c1, c0, c2)
    a1 = angle_between(c2, c1, c3)
    a2 = angle_between(c3, c2, c0)
    a3 = angle_between(c0, c3, c1)
    return a0 + a1 + a2 + a3 - 2. * numpy.pi
