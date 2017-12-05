#!/bin/bash

fyear=1980
lyear=1980

for year in $(seq $fyear $lyear ) ; do

    python make_clim_file.py $year

    cd clim
    ncrcat CCS_clim_${year}_??_SODA3.3.1.nc -o CCS_clim_SODA3.3.1_y${year}.nc
    rm CCS_clim_${year}_??_SODA3.3.1.nc
    cd ..

done

