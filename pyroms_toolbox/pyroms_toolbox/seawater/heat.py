# --- encoding: iso-8859-1 ---

"""Seawater heat module

heatcap(S, T[, P])       -- Heat capacity
adtgrad(S, T[, P])       -- Adiabatic temperature gradiente
temppot(S, T, P[, Pref]) -- Potential temperature
temppot0(S, T, P)        -- Potential temperature, relative to surface

Bjørn Ådlandsvik, <bjorn@imr.no>, 07 November 2004

"""

# -------------------------------------------------

def heatcap(S, T, P=0):
    """Compute heat capacity

    Usage: heatcap(S, T, [P])

    Input:
        S = Salinity,     [PSS-78]
        T = Temperature,  [°C]
        P = Pressure,     [dbar]
    P is optional, with a default value = zero

    Output:
        Heat capacity  [J/(kg*K)]

    Algorithm: UNESCO 1983

    """
    
    P = 0.1*P  # Conversion to bar

    # - Temperatur dependence
    c0 =  4217.4
    c1 = -3.720283
    c2 =  0.1412855
    c3 = -2.654387e-3
    c4 =  2.093236e-5

    a0 = -7.64357
    a1 =  0.1072763
    a2 = -1.38385e-3

    b0 =  0.1770383
    b1 = -4.07718e-3
    b2 =  5.148e-5

    CP0 =  c0 + c1*T + c2*T**2 + c3*T**3 + c4*T**4  \
          + (a0 + a1*T + a2*T**2)*S \
          + (b0 + b1*T + b2*T**2)*S**1.5

    # - Pressure dependence    
    a0 = -4.9592e-1
    a1 =  1.45747e-2
    a2 = -3.13885e-4
    a3 =  2.0357e-6
    a4 =  1.7168e-8

    b0 =  2.4931e-4
    b1 = -1.08645e-5
    b2 =  2.87533e-7
    b3 = -4.0027e-9
    b4 =  2.2956e-11

    c0 = -5.422e-8
    c1 =  2.6380e-9
    c2 = -6.5637e-11
    c3 =  6.136e-13

    CP1 = (a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4)*P  \
              + (b0 + b1*T + b2*T**2 + b3*T**3 + b4*T**4)*P**2 \
              + (c0 + c1*T + c2*T**2 + c3*T**3)*P**3

    # - Salinity dependence    
    d0 =  4.9247e-3
    d1 = -1.28315e-4
    d2 =  9.802e-7
    d3 =  2.5941e-8
    d4 = -2.9179e-10

    e0 = -1.2331e-4
    e1 = -1.517e-6
    e2 =  3.122e-8

    f0 = -2.9558e-6
    f1 =  1.17054e-7
    f2 = -2.3905e-9
    f3 =  1.8448e-11

    g0 =  9.971e-8
 
    h0 =  5.540e-10
    h1 = -1.7682e-11
    h2 =  3.513e-13

    j1 = -1.4300e-12
    S3_2  = S**1.5

    CP2 = ((d0 + d1*T + d2*T**2 + d3*T**3 + d4*T**4)*S \
           + (e0 + e1*T + e2*T**2)*S3_2)*P  \
	   + ((f0 + f1*T + f2*T**2 + f3*T**3)*S  \
	   +   g0*S3_2)*P**2 \
	   + ((h0 + h1*T + h2*T**2)*S + j1*T*S3_2)*P**3
     

    return CP0 + CP1 + CP2

# --------------------------------------------------------------

def adtgrad(S, T, P=0):
    """Compute adiabatic temperature gradient

    Usage: adtgrad(S, T, [P])

    Input:
        S = Salinity,     [PSS-78]
        T = Temperature,  [°C]
        P = Pressure,     [dbar]
    P is optional, with a default value = zero

    Output:
        Adiabatic temperature gradient,  [K/dbar]

    Algorithm: UNESCO 1983

    """
    
    a0 =  3.5803e-5
    a1 = +8.5258e-6
    a2 = -6.836e-8
    a3 =  6.6228e-10

    b0 = +1.8932e-6
    b1 = -4.2393e-8

    c0 = +1.8741e-8
    c1 = -6.7795e-10
    c2 = +8.733e-12
    c3 = -5.4481e-14

    d0 = -1.1351e-10
    d1 =  2.7759e-12

    e0 = -4.6206e-13
    e1 = +1.8676e-14
    e2 = -2.1687e-16

    return  a0 + (a1 + (a2 + a3*T)*T)*T  \
         + (b0 + b1*T)*(S-35)  \
	 + ( (c0 + (c1 + (c2 + c3*T)*T)*T) \
         +   (d0 + d1*T)*(S-35) )*P \
         + (e0 + (e1 + e2*T)*T )*P*P

# ---------------------------------------------------------------

def temppot(S, T, P, Pref=0):
    """Compute potential temperature

    Usage: temppot(S, T, P, [Pref])

    Input:
        S = Salinity,                [PSS-78]
        T = Temperature,             [°C]
        P = Pressure,                [dbar]
        Pref = Reference pressure,   [dbar]
    Pref is optional, with a default value = zero

    Output:
        Potential temperature,  [°C]

    Algorithm: UNESCO 1983

    """
    
    H = Pref-P
    XK = H*adtgrad(S,T,P)

    T = T + 0.5*XK
    Q = XK
    P = P + 0.5*H
    XK = H*adtgrad(S,T,P)

    T = T + 0.29289322*(XK-Q)
    Q = 0.58578644*XK + 0.121320344*Q
    XK = H*adtgrad(S,T,P)
    
    T = T + 1.707106781*(XK-Q)
    Q = 3.414213562*XK - 4.121320344*Q
    P = P + 0.5*H
    XK = H*adtgrad(S,T,P)
    
    return T + (XK-2.0*Q)/6.0    
    

# ------------------------------------------------------

def temppot0(S,T,P):
    """Compute potential temperature relative to surface

    Usage: temppot0(S, T, P)

    Input:
        S = Salinity,                [PSS-78]
        T = Temperature,             [°C]
        P = Pressure,                [dbar]
 
    Output:
        Potential temperature,       [°C]

    Algorithm: Bryden 1973

    Note: Due to different algorithms,
        temppot0(S, T, P) != tempot(S, T, P, Pref=0)
  
    """

    P = P/10  # Conversion from dbar

    a0 =  3.6504e-4
    a1 =  8.3198e-5
    a2 = -5.4065e-7
    a3 =  4.0274e-9

    b0 =  1.7439e-5
    b1 = -2.9778e-7

    c0 =  8.9309e-7
    c1 = -3.1628e-8
    c2 =  2.1987e-10

    d0 =  4.1057e-9

    e0 = -1.6056e-10
    e1 =  5.0484e-12

    S0 = S - 35.0

    return  T - (a0 + (a1 + (a2 + a3*T)*T)*T)*P  \
              - (b0 + b1*T)*P*S0                 \
              - (c0 + (c1 + c2*T)*T)*P*P         \
              + d0*S0*P*P                        \
              - (e0 + e1*T)*P*P*P                


