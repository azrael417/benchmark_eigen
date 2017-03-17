#!/bin/bash
#SBATCH --partition=workq
#SBATCH --nodes=16
#SBATCH --ntasks=512
#SBATCH -C BW28
#SBATCH --time=01:30:00
#SBATCH -J PDSY_F90

source /opt/modules/default/init/bash
# module swap PrgEnv-cray PrgEnv-gnu
module list

PDSYROOT=/cray/css/users/pcarrier/Applications/Z_others/eigen/pdsy_samples/FINAL

ulimit -s unlimited
ulimit -a

NPROW=32
NPCOL=16
RANKS=$(($NPROW*$NPCOL))
SOLVER_TYPE="PDSYEVD"
# SOLVER_TYPE="ELPA"

echo '# of MPI ranks=' $RANKS

MATRIX_SIZE=8192

export OMP_NUM_THREADS=1
echo ======== start ==============
date
echo ======== start ==============
  time srun --ntasks $RANKS \
       $PDSYROOT/eigen.x $NPROW $NPCOL $MATRIX_SIZE $SOLVER_TYPE\
       1> eigen.log \
echo ======== end ==============
date
echo ======== end ==============

