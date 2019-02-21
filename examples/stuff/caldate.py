#/usr/bin/env python

# This script is used to change a Julian date into a dictionary containing
# the Gregorian equivalent.

from math import *

def leapyear(year):
    """
    Returns 1 if the provided year is a leap year, 0 if the provided
    year is not a leap year.
    """
    if year % 4 == 0:
        if year % 100 == 0:
            if year % 400 == 0:
                return 1
            else:
                return 0
        else:
            return 1
    else:
        return 0
    
def caldate_1900(Julian):
    """
    This is nearly a direct translation of a Matlab script to a Python script
    for changing a Julian date into a Gregorian date. The algorithm was not
    created by me.
    """
    JulDay = Julian
    if (JulDay < 2440000):
        JulDay = JulDay+2415020+1

    # This is meant to prevent round-off
    JulDay = JulDay+5.0e-9

    # Conversion to a Gregorian date
    j = floor(JulDay)-1721119
    
    jin = 4*j-1
    
    y = floor(jin/146097)

    j = jin-146097*y

    jin = floor(j/4)

    jin = 4*jin+3

    j = floor(jin/1461)

    d = floor(((jin-1461*j)+4)/4)
    jin = 5*d-3

    m = floor(jin/153)
    d = floor(((jin-153*m)+5)/5)
    y = y*100+j
    
    if m < 10:
        mo = m+3
        yr = y
    else:
        mo = m-9
        yr = y+1

    ivd = [1,32,60,91,121,152,182,213,244,274,305,335,366]
    ivdl = [1,32,61,92,122,153,183,214,245,275,306,336,367]

    if (leapyear(yr) == 1):
        yday = ivdl[int(mo-1)]+d-1
    else:
        yday = ivd[int(mo-1)]+d-1


    secs = (JulDay%1)*24*3600
    sec = round(secs)
    hour = floor(sec/3600)
    min = floor((sec%3600)/60)
    sec = round(sec%60)

    print("Year: "+str(yr))
    print("Year Day: "+str(yday))
    print("Month: "+str(mo))
    print("Day: "+str(d))
    print("Hour: "+str(hour))
    print("Min: "+str(min))
    print("Sec: "+str(sec))

    cal = {'year':yr,'yearday':yday,'month':mo,'day':d,\
           'hour':hour,'minute':min,'second':sec}

    return cal

if __name__ == "__main__":

    mycal = caldate_1900(36761.5)
    # You can get any of the dictionary values by typing mycal['key'] where
    # key is replaced by the key to the hash table such as 'hour'
    # print mycal['year']
