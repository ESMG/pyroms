import xarray as xr
from regrid_GLORYS import regrid_GLORYS

gsource = xr.open_dataset('/import/AKWATERS/kshedstrom/glorys/GLORYS_REANALYSIS_2018-01-01.nc')
myssh = regrid_GLORYS(gsource.zos, method='bilinear')
myssh.to_netcdf('myssh.nc')
