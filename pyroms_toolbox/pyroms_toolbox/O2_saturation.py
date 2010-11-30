# encoding: utf-8

def O2_saturation(T,S):
    """
    Return saturation value of oxygen.
    
    Parameters
    ----------
    T : array_like
        Temperature (˚C)
    S : array_like
        Salinity (PSU)
    
    Returns
    -------
    O2_sat : array_like
        concentrations of O2 [ millimole O2 / m3 ] for a given temperature and
        salinity (at STP)
    """

    A1 = -173.4292
    A2 =  249.6339
    A3 =  143.3483
    A4 =  -21.8492
    B1 =   -0.033096
    B2 =    0.014259
    B3 =   -0.0017000
    # Convert T to deg. C to deg. K
    T = T + 273.15
    # O2 Concentration in mg/l
    # [from Millero and Sohn, Chemical Oceanography, CRC Press, 1992]
    O = np.exp(A1 + A2*(100.0/T) + A3*np.log(T/100.0) + A4*(T/100.0) + \
            S*(B1 + B2*(T/100.0) + B3*((T/100.0)**2)) )
    # Convert to mmol/m3
    #  mmol/m3 = 44.66 ml/l
    #  mg/l = ml/l * 1.42903 mg/ml
    return O*(44.66*1.42903)
