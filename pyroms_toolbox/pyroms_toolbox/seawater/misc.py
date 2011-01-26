# --- encoding: iso-8859-1 ---

"""Miscellaneous sea water functions

freezept(S[, P])    -- Freezing point
soundvel(S, T[, P]) -- Sound velocity
depth(P, lat)       -- Depth from pressure

Bjørn Ådlandsvik, <bjorn@imr.no>  07 November 2004

"""

# ----------------------------------------------------------

def freezept(S, P=0):

    """Compute freezing temperature of sea water

    Usage: freezept(S, [P])
 
    Input:               
        S = Salinity,      [psu]
        P = Pressure,      [dbar]
    P is optional, with a default value = 0

    Output:
        T = Freezing point,   [°C]

    Algorithm: UNESCO 1983 

    """
    
    a0 = -0.0575
    a1 =  1.710523e-3
    a2 = -2.154996e-4
    b  = -7.53e-4

    Tf = a0*S + a1*S**1.5 + a2*S**2 + b*P
    return Tf

# ----------------------------------------------------------------

def soundvel(S, T, P=0):
    """Compute velocity of sound

    Usage: soundvel(S, T, [P])

    Input:
        S = Salinity,     [PSS-78]
        T = Temperature,  [°C]
        P = Pressure,     [dbar]
    P is optional, with a default value = zero

    Output:
        Sound velocity,  [m/s]

    Algorithm: UNESCO 1983 

    """

    P = 0.1*P  # Conversion to bar

    c00 = 1402.388
    c01 =  5.03711
    c02 = -5.80852e-2
    c03 =  3.3420e-4
    c04 = -1.47800e-6
    c05 =  3.1464e-9

    c10 =  0.153563
    c11 =  6.8982e-4
    c12 = -8.1788e-6
    c13 =  1.3621e-7
    c14 = -6.1185e-10

    c20 =  3.1260e-5
    c21 = -1.7107e-6
    c22 =  2.5974e-8
    c23 = -2.5335e-10
    c24 =  1.0405e-12

    c30 = -9.7729e-9
    c31 =  3.8504e-10
    c32 = -2.3643e-12
 
    P2 = P*P
    P3 = P2*P
    Cw =  c00 + (c01 + (c02 + (c03 + (c04 + c05*T)*T)*T)*T)*T   \
       + (c10 + (c11 + (c12 + (c13 + c14*T)*T)*T)*T)*P          \
       + (c20 + (c21 + (c22 + (c23 + c24*T)*T)*T)*T)*P2         \
       + (c30 + (c31 + c32*T)*T)*P3
 
    a00 =  1.389
    a01 = -1.262e-2
    a02 =  7.164e-5
    a03 =  2.006e-6
    a04 = -3.21e-8

    a10 =  9.4742e-5
    a11 = -1.2580e-5
    a12 = -6.4885e-8
    a13 =  1.0507e-8
    a14 = -2.0122e-10

    a20 = -3.9064e-7
    a21 =  9.1041e-9
    a22 = -1.6002e-10
    a23 =  7.988e-12

    a30 =  1.100e-10
    a31 =  6.649e-12
    a32 = -3.389e-13

    A =  a00 + (a01 + (a02 + (a03 + a04*T)*T)*T)*T      \
      + (a10 + (a11 + (a12 + (a13 + a14*T)*T)*T)*T)*P   \
      + (a20 + (a21 + (a22 + a23*T)*T)*T)*P2            \
      + (a30 + (a31 + a32*T)*T)*P3

    b00 = -1.922e-2
    b01 = -4.42e-5
    b10 =  7.3637e-5
    b11 =  1.7945e-7

    B = b00 + b01*T + (b10 + b11*T)*P

    d00 =  1.727e-3
    d10 = -7.9836e-6

    D = d00 + d10*P

    return Cw + A*S + B*S**1.5 + D*S**2

# ----------------------------------------------------------------

def depth(P, lat):
    """Compute depth from pressure and latitude

    Usage: depth(P, lat)

    Input:
        P = Pressure,     [dbar]
        lat = Latitude    [deg]

    Output:
        Depth             [m]

    Algorithm: UNESCO 1983 

    """
    
    # Use Numeric for trigonometry if present
    try:
        from Numeric import sin, pi
    except:
        from math import sin, pi
    
    a1 =  9.72659
    a2 = -2.2512e-5
    a3 =  2.279e-10
    a4 = -1.82e-15

    b  =  1.092e-6

    g0 =  9.780318
    g1 =  5.2788e-3
    g2 =  2.36e-5

    rad = pi / 180.

    X = sin(lat*rad)
    X = X*X
    grav = g0 * (1.0 + (g1 + g2*X)*X) + b*P
    nom = (a1 + (a2 + (a3 + a4*P)*P)*P)*P

    return nom / grav
                                      



