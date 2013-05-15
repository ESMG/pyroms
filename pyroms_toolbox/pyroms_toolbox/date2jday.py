import numpy as np

def date2jday(date):
    """
    return days since 1900-01-01

    jday = date2jday(date) where date is a datetime object
    """

    if type(date).__name__ == 'ndarray':
        nt = np.size(date)
    else:
        nt = 1
        date = np.array([date])

    jday = []

    for t in range(nt):
        a = np.floor((14 - date[t].month)/12.)
        y = date[t].year + 4800 - a
        m = date[t].month + 12*a - 3

        # For a date in the Gregorian calendar:
        jday.append(date[t].day + np.floor((153*m + 2)/5) \
                    + y*365 + np.floor(y/4.) - np.floor(y/100.) \
                    + np.floor(y/400.) - 32045 \
                    + ( date[t].second + 60*date[t].minute + 3600*(date[t].hour - 12) )/86400. \
                    - 2415020.5)

    if np.size(jday) == 1:
        jday = jday[0]
    else:
        jday = np.array(jday)

    return jday
