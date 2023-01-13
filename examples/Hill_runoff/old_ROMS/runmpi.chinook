#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=runoff-remap
#SBATCH --account=akwaters
#SBATCH --output=bdry.%j
#SBATCH --no-requeue
#SBATCH -p t1small

cd $SLURM_SUBMIT_DIR
. /usr/share/Modules/init/bash
module purge
# Gnu universe
module load toolchain/foss/2016b
#
#module load slurm/22.05.4
#module load foss/2022a
#module load netCDF-Fortran/4.5.4
module list
source activate snowdrift
#source activate snowy

#
#  Prolog
#
echo " "
echo "++++ Chinook ++++ $PGM_NAME began:    `date`"
echo "++++ Chinook ++++ $PGM_NAME hostname: `hostname`"
echo "++++ Chinook ++++ $PGM_NAME uname -a: `uname -a`"
echo " "
TBEGIN=`echo "print time();" | perl`

which python
#for year in {1993..2013..2}
#do
#  python make_river_file.py ${year}
#  python add_rivers.py Hill_rivers_${year}.nc
#  python make_river_clim.py discharge_*${year}_q_Hill_NWGOA.nc Hill_rivers_${year}.nc
#  python set_vshape.py Hill_rivers_${year}.nc
#done
python make_river_file.py 2014
python add_rivers.py Hill_rivers_2014.nc
python make_river_clim.py discharge_2014_q_Hill_NWGOA.nc Hill_rivers_2014.nc
python set_vshape.py Hill_rivers_2014.nc
#

#
#  Epilog
#
TEND=`echo "print time();" | perl`
echo " "
echo "++++ Chinook ++++ $PGM_NAME pwd:      `pwd`"
echo "++++ Chinook ++++ $PGM_NAME ended:    `date`"
echo "++++ Chinook ++++ $PGM_NAME walltime: `expr $TEND - $TBEGIN` seconds"
