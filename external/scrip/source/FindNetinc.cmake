# This module for CMake is used to find the netcdf.inc file from a user-specified 
# environment variable.

set(Netcdfinc_FOUND FALSE)

find_path(Netcdf_INC netcdf.inc "$ENV{NETCDFINC}")

if (Netcdf_INC)
    set(Netcdfinc_FOUND TRUE)
endif(Netcdf_INC)
