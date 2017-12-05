#!/bin/bash

fyear=1980
lyear=1980

for year in $(seq $fyear $lyear ) ; do

    python make_bdry_file.py $year

    cd bdry
    ncrcat CCS_bdry_${year}_??_??_SODA3.3.1.nc -o CCS_bdry_SODA3.3.1_y${year}.nc
    rm CCS_bdry_${year}_??_??_SODA3.3.1.nc
    cd ..

done

