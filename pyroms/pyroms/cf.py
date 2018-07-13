# encoding: utf-8
"""
cf.py - classes around CF compliant files

The cf module is made for reading CF-compliant datasets,
knowing data, its structure, units and conversions
between units afterwards.

Dependencies:
=============
numpy
netcdftime (packaged in netcdf4-python)
"""
__docformat__ = "restructuredtext en"

import numpy as np

try:
  import cftime
except:
  import netcdftime

import pyroms.io

class time (np.ndarray):
    """Return time object from netCDF file

    Parameters
    ----------
    nc : netCDF3/4 object or filename
        Time information will be read from this netCDF3/4 file.
    name : string, optional
        The name of the the variable.
    units : string, optional
        The name of the variable units.
    calendar : string, optional
        A string representing the calandar to use. See netcdftime
        documentation for possible values.

    Returns
    -------
    nctime : ndarray
        A subclass of numpy.ndarray with values equal to the time variable in
        the netCDF file referenced with nc.

    """

    _unit2sec={'seconds' : 1.0,
               'minutes' : 60.0,
               'hours'   : 3600.0,
               'days'    : 3600.0*24.0,
               'weeks'   : 3600.0*24.0*7.0,
               'years'   : 3600.0*24.0*365.242198781} #ref to udunits

    _sec2unit={'seconds' : 1.0,
               'minutes' : 1.0/60.0,
               'hours'   : 1.0/3600.0,
               'days'    : 1.0/(24.0*3600.0)}

    def __new__(self, ncfile, name='time', units=None, calendar='standard'):
        self._nc = pyroms.io.Dataset(ncfile)
        data = self._nc.variables[name][:]
        data = data.view(time)
        if units == None:
            units = self._nc.variables[name].units
        data.utime = netcdftime.utime(units, calendar=calendar)
        return data

    def __array_finalize__(self, obj):
        self.utime = getattr(obj, 'utime', {})

    def arg_nearest_date(self, dateo):
        """Return index of date nearest to query date.

        Prameters
        ---------
        dateo : datetime object
            The query date

        Returns
        -------
        idx : integer
            The index of the date closest to dateo. If two dates are
            equidistant, the smaller is returned.

        """
        to = self.utime.date2num(dateo)
        return np.min(np.where(np.abs(self-to) == \
                np.min(np.abs(self-to)))[0])

    def nearest_date(self, dateo):
        """Return the nearest date to query date.

        Prameters
        ---------
        dateo : datetime object
            The query date

        Returns
        -------
        nearest_date : datetime object
            A datetime object of the date closest to dateo. If two dates are
            equidistant, the smaller is returned.

        """
        idx = np.where(np.abs(self.dates-dateo) == \
                np.min(np.abs(self.dates-dateo)))[0]
        idx = np.min(idx)
        return self.dates[idx]

    def arg_nearest(self, to, units=None):
        """Return index of time nearest to query time.

        Prameters
        ---------
        to : float
            The query time.
        units : string, optional
            The units of the reference time. Defaults to the reference time
            string 'units' in the netcdf oject.

        Returns
        -------
        idx : integer
            The index of the date closest to to. If two times are equidistant,
            the smaller is returned.

        """
        if units is not None:
            to *= self._unit2sec[units] * self._sec2unit[self.utime.units]
        return np.min(np.where(np.abs(self-to) == np.min(np.abs(self-to)))[0])

    def nearest(self, to, units=None):
        """Return time nearest to time query.

        Prameters
        ---------
        to : float
            The query time.
        units : string, optional
            The units of the reference time. Defaults to the reference time
            string 'units' in the netcdf oject.

        Returns
        -------
        idx : integer
            The index of the date closest to to. If two times are equidistant,
            the smaller is returned.

        """
        if units is not None:
            to *= self._unit2sec[units] * self._sec2unit[self.utime.units]
        idx = np.where(np.abs(self-to) == np.min(np.abs(self-to)))[0]
        idx = np.min(idx)
        return self[idx]

    def get_seconds(self):
        fac = self._unit2sec[self.utime.units] * self._sec2unit['seconds']
        return self*fac

    def get_minutes(self):
        fac = self._unit2sec[self.utime.units] * self._sec2unit['minutes']
        return self*fac

    def get_hours(self):
        fac = self._unit2sec[self.utime.units] * self._sec2unit['hours']
        return self*fac

    def get_days(self):
        fac = self._unit2sec[self.utime.units] * self._sec2unit['days']
        return np.asarray(self,dtype='float64')*fac

    def get_jd(self):
        utime = netcdftime.utime('days since 0001-01-01 00:00:00', \
                                 calendar='proleptic_gregorian')
        return utime.date2num(self.dates)

    def get_dates(self):
        return np.array([self.utime.num2date(tval) for tval in self])

    jd = property(get_jd, None, doc="Julian day, for plotting in pylab")
    seconds = property(get_seconds, None, doc="seconds")
    minutes = property(get_minutes, None, doc="minutes")
    hours = property(get_hours, None, doc="hours")
    days = property(get_days, None, doc="days")
    dates = property(get_dates, None, doc="datetime objects")
