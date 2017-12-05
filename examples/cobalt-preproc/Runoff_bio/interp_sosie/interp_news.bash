#!/bin/bash

export PATH=$PATH:/opt/esm-soft/sosie/bin/

for varname in NO3 LDON SLDON SRDON NDET PO4 LDOP SLDOP SRDOP ; do

    cat nml.sosie.skel | sed -e "s/<VARNAME>/$varname/g" > nml.sosie.$varname
    sosie.x -f nml.sosie.$varname
    rm nml.sosie.$varname

done

# merge together
mv NO3_NEWS-CCS1_clim.nc NEWS-CCS1_clim.nc

for varname in LDON SLDON SRDON NDET PO4 LDOP SLDOP SRDOP ; do
    ncks -A -v $varname ${varname}_NEWS-CCS1_clim.nc -o NEWS-CCS1_clim.nc
    rm ${varname}_NEWS-CCS1_clim.nc
done
