import xarray as xr
from regrid_GLBy import regrid_GLBy

gsource = xr.open_dataset('/import/AKWATERS/kshedstrom/HYCOM/Svalbard/data/HYCOM_GLBy0.08_2018_345.nc')
myssh = regrid_GLBy(gsource.ssh, method='bilinear')
myssh.to_netcdf('myssh.nc')
