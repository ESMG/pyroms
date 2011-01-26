import numpy as np

def date2jday(date):
    """
    return days since 1900-01-01

    jday = date2jday(date) where date is a datetime object
    """

    a = np.floor((14 - date.month)/12.)
    y = date.year + 4800 - a
    m = date.month + 12*a - 3

    # For a date in the Gregorian calendar:
    jday = date.day + np.floor((153*m + 2)/5) \
             + y*365 + np.floor(y/4.) - np.floor(y/100.) \
             + np.floor(y/400.) - 32045 \
             + ( date.second + 60*date.minute + 3600*(date.hour - 12) )/86400. \
             - 2415020.5

    return jday
