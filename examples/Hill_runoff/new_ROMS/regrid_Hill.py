import xarray as xr
import xesmf
import numpy as np

def regrid_Hill(fld, method='nearest_d2s', irange=None, jrange=None):
    dst_grd = xr.open_dataset('/import/AKWATERS/kshedstrom/gridpak/Cook_inlet/Cook_Inlet_grid_1.nc')
    dst_grd = dst_grd.rename({'lon_rho': 'lon', 'lat_rho': 'lat'})
#   p_points = xr.open_dataset('p_points.nc')
#   dst_grd = dst_grd.merge(p_points)
#   dst_grd = dst_grd.rename({'lon_p': 'lon_b', 'lat_p': 'lat_b'})

    gsource = xr.open_dataset('/import/AKWATERS/kshedstrom/hydroflow/new_2019/goa_dischargex_09012014_08312015.nc')
    print(gsource["timeSeries"].shape)
#   area = np.ones((gsource["timeSeries"].shape)) * 1.e6
#   newVar = xr.Dataset({'area': (("timeSeries"), area),}, coords={"timeSeries": gsource["timeSeries"]})
#   gsource = gsource.merge(newVar)
#   gsource = gsource.rename({'longitude': 'lon', 'latitude': 'lat'})

    if irange is not None:
        gsource = gsource.isel(lon=slice(irange[0],irange[1]))

    if jrange is not None:
        gsource = gsource.isel(lat=slice(jrange[0],jrange[1]))

    regrid = xesmf.Regridder(
        gsource,
        dst_grd,
        method=method,
        periodic=False,
        filename='regrid_t.nc',
        reuse_weights=True
    )
    tdest = regrid(fld)
    return tdest
