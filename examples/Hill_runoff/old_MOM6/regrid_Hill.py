import xarray as xr
import xesmf

def regrid_Hill(fld, coords, method='nearest_d2s', irange=None, jrange=None):

    gsource = xr.open_dataset('/import/c1/AKWATERS/kshedstrom/Runoff/discharge_1992_1993.nc')
#   p_points = xr.open_dataset('p_points.nc')
#   gsource = gsource.merge(p_points)
#   gsource = gsource.rename({'longitude': 'lon', 'latitude': 'lat'})

    if irange is not None:
        gsource = gsource.isel(lon=slice(irange[0],irange[1]))

    if jrange is not None:
        gsource = gsource.isel(lat=slice(jrange[0],jrange[1]))

    regrid = xesmf.Regridder(
        gsource,
        coords,
        method=method,
        periodic=False,
        filename='regrid_t.nc',
        reuse_weights=True
    )
    tdest = regrid(fld)
    return tdest
