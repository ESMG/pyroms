#!/bin/bash

echo -n "Enter a year to process: ";read year; echo ""


python get_MERRA_Tair_from_nasa_opendap_3hours.py $year
python get_MERRA_Pair_from_nasa_opendap_3hours.py $year
python get_MERRA_Qair_from_nasa_opendap_3hours.py $year
python get_MERRA_Uwind_from_nasa_opendap_3hours.py $year
python get_MERRA_Vwind_from_nasa_opendap_3hours.py $year
python get_MERRA_lwrad_down_from_nasa_opendap_3hours.py $year
python get_MERRA_swrad_from_nasa_opendap_3hours.py $year
python get_MERRA_rain_from_nasa_opendap_3hours.py $year
python get_MERRA_snow_from_nasa_opendap_3hours.py $year
python get_MERRA_cloud_from_nasa_opendap_3hours.py $year
python get_MERRA_albedo_from_nasa_opendap_daily.py $year

