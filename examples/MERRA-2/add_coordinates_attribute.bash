lst=`ls Forcings/MERRA_Pair_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,Pair,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_Qair_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,Qair,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_Tair_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,Tair,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_Uwind_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,Uwind,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_Vwind_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,Vwind,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_cloud_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,cloud,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_lwrad_down_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,lwrad_down,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_rain_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,rain,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_snow_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,snow,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_swrad_3hours_*`
for file in $lst ; do
ncatted -O -a coordinates,swrad,o,c,'lon lat' $file
done

lst=`ls Forcings/MERRA_albedo_daily_*`
for file in $lst ; do
ncatted -O -a coordinates,albedo,o,c,'lon lat' $file
done

