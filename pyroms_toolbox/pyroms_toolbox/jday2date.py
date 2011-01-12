import numpy as np
from decimal import *
from datetime import datetime

def jday2date(jday):
    """
    return the date from the number od days since 1900/01/01.

    date = jday2date(jday) where jday is an integer or a float
                             and date is a datetime object
    """

    jd0 = 2415021 #days since 1900/01/01


    j = int(np.floor(jday)) + 32044 + jd0 
    g = j / 146097 
    dg = j % 146097 
    c = (dg / 36524 + 1) * 3 / 4 
    dc = dg - c * 36524 
    b = dc / 1461 
    db = dc % 1461 
    a = (db / 365 + 1) * 3 / 4 
    da = db - a * 365 
    y = g * 400 + c * 100 + b * 4 + a 
    m = (da * 5 + 308) / 153 - 2 
    d = da - (m + 4) * 153 / 5 + 122 
    Y = y - 4800 + (m + 2) / 12 
    M = (m + 2) % 12 + 1 
    D = d + 1 

    hr = Decimal(str(jday)) - Decimal(np.floor(jday))
    hr = hr * 24
    h = int(hr)
    m = int((hr-h)*60)
    s = int(hr*3600 - h*3600 - m*60)

    return datetime(Y,M,D,h,m,s)
