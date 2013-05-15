import datetime

def gregorian2iso(date):
    "ISO year, week and day for a given Gregorian calendar date"
    return date.isocalendar()
