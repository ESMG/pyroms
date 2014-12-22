#!/usr/bin/ksh 

print -n "Enter year to proceed: "; read year; print ""

leap=`echo $(($year % 4))`

if [ $leap == 0 ] ; then
    nday=366
else
    nday=365
fi
#nday=230

#set -A days {258..258}
set -A days {1..$nday}

for day in ${days[@]} ; do
    day=`echo $day | awk '{printf "%03d", $1}'`
    ncks -a -O HYCOM_GLBa0.08_ssh_${year}_${day}.nc HYCOM_GLBa0.08_${year}_${day}.nc
    ncks -a -A HYCOM_GLBa0.08_temp_${year}_${day}.nc HYCOM_GLBa0.08_${year}_${day}.nc
    ncks -a -A HYCOM_GLBa0.08_salt_${year}_${day}.nc HYCOM_GLBa0.08_${year}_${day}.nc
    ncks -a -A HYCOM_GLBa0.08_u_${year}_${day}.nc HYCOM_GLBa0.08_${year}_${day}.nc
    ncks -a -A HYCOM_GLBa0.08_v_${year}_${day}.nc HYCOM_GLBa0.08_${year}_${day}.nc
    rm -f HYCOM_GLBa0.08_ssh_${year}_${day}.nc HYCOM_GLBa0.08_temp_${year}_${day}.nc HYCOM_GLBa0.08_salt_${year}_${day}.nc HYCOM_GLBa0.08_u_${year}_${day}.nc HYCOM_GLBa0.08_v_${year}_${day}.nc
done
