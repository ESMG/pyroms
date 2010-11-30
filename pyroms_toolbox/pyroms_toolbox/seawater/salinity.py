# --- encoding: iso-8859-1 ---

"""Seawater salinity module, providing salt and cond functions

S = salt(R, T[, P]) -- Salinity    
R = cond(S, T[, P]) -- Conductivity ratio

Arguments:
    R = Conductivity ratio
    S = Salinity
    T = Temperature         [°C]
    P = Pressure,           [dbar = 10**4 Pa]

Bjørn Ådlandsvik <bjorn@imr.no>, 07 November 2004

"""

# -------------------------------------------------------

def _sal(XR,XT):
    
    a0 =  0.0080
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081
    
    b0 =  0.0005
    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144
    
    k  =  0.0162

    DS = (XT / (1+k*XT) ) *        \
         (b0 + (b1 + (b2 + (b3 + (b4 + b5*XR)*XR)*XR)*XR)*XR)

    return a0 + (a1 + (a2 + (a3 + (a4 + a5*XR)*XR)*XR)*XR)*XR + DS

# ---------------------------------------------------

def _dsal(XR,XT):
    
    a1 = -0.1692
    a2 = 25.3851
    a3 = 14.0941
    a4 = -7.0261
    a5 =  2.7081

    b1 = -0.0056
    b2 = -0.0066
    b3 = -0.0375
    b4 =  0.0636
    b5 = -0.0144

    k  =  0.0162
    
    dDS = (XT / (1+k*XT) ) *      \
          (b1 + (b2*2 + (b3*3 + (b4*4 + b5*5*XR)*XR)*XR)*XR)
    
    return a1 + (a2*2 + (a3*3 + (a4*4 + a5*5*XR)*XR)*XR)*XR + dDS

# ---------------------------------------------

def _rt(T):
    
    c0 =  0.6766097
    c1 =  2.00564e-2
    c2 =  1.104259e-4
    c3 = -6.9698e-7
    c4 =  1.0031e-9
    
    return c0 + (c1 + (c2 + (c3 + c4*T)*T)*T)*T

# ---------------------------------------------------

def _c(P):
    
    e1 =  2.070e-5
    e2 = -6.370e-10
    e3 =  3.989e-15
    
    return (e1 + (e2 + e3*P)*P)*P

# ---------------------------------------------------

def _b(T):
    
    d1 =  3.426e-2
    d2 =  4.464e-4
    
    return 1.0 + (d1 + d2*T)*T

# ---------------------------------------------------

def _a(T):
    
    d3 =  4.215e-1
    d4 = -3.107e-3
    
    return d3 + d4*T

# --------------------------------------------------

def salt(R, T, P):
    """Compute salinity from conductivity, temperature, and pressure

    Usage: salt(R, T, [P])

    Input:
        R = Conductivity ratio
        T = Temperature        [°C]
        P = Pressure,          [dbar = 10**4 Pa]
    P is optional, with default value zero
      
    Output:
        S = Salinity           [PSS-78]

    """
    
    DT = T - 15.0
    RT = R/(_rt(T)*(1.0 + _c(P)/(_b(T) + _a(T)*R)))
    RT = abs(RT)**0.5
    
    return _sal(RT,DT)

# -------------------------------------------------

def cond(S, T, P):
    """Compute conductivity ratio from salinity, temperature, and pressure

    Usage: cond(S, T, [P])

    Input:
        S = Salinity      [PSS-78]
        T = Temperature   [°C]
        P = Pressure,     [dbar = 10**4 Pa]
    P is optional, with default value zero
      
    Output:
        R = Conductivity ratio

"""
    
    DT = T-15.0
    RT = (S/35.0)**0.5
    SI = _sal(RT,DT)
    # Iteration
    for n in xrange(100):
        RT = RT + (S-SI)/_dsal(RT,DT)
        SI = _sal(RT,DT)
        try:
            DELS = max(abs(SI-S))
        except TypeError: # Not sequence, i.e. scalar S
            DELS = abs(SI-S)
        if (DELS < 1.0E-4):
            break

    RTT = _rt(T)*RT*RT
    AT = _a(T)
    BT = _b(T)
    CP = _c(P)
    CP = RTT*(CP + BT)
    BT = BT - RTT*AT

    R = abs(BT*BT + 4.0*AT*CP)**0.5 - BT
    
    return 0.5*R/AT


