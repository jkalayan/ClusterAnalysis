#!/bin/bash --login
#$ -cwd
#$ -S /bin/bash
#$ -N run-clusters
#$ -pe smp.pe 8    # Each task will use X cores


## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## load modules
module load apps/anaconda3/5.2.0
## assign slots 
echo "NSLOTS" $NSLOTS
export OMP_NUM_THREADS=$NSLOTS
## set path for code
EXE="$HOME/bin/ClusterAnalysis/runClusterAnalysis.py" # path to folder containing code
## make sure all .py files have the correct permissions
chmod +x ${EXE}
for dir in "$HOME/bin/ClusterAnalysis/ClusterAnalysis/*"; do
	chmod +x ${dir}/*.py
done
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


## set variables for code

DIR="/mnt/iusers01/rb01/m17206an/test_directory" # path to where topology and trajectoy files are

TOPOLOGY="Mol3tc.top" # topology file name
TRAJECTORY="equil.mdcrd" # trajectory file name

START=0 # starting frame for analysis
END=100 # end frame for analysis
STEP=1 # step in-between frames to analyse

CUTOFF=3 # cutoff distance in Angstrom
SOLVENT="WAT" # name of solvent molecules
EXCLUDE="Cl- Na+ DMS" # name of molecules to exclude from cluster analysis

## setup instructions for running code and execute
EXE2="${EXE} -top ${DIR}/${TOPOLOGY} -traj ${DIR}/${TRAJECTORY} -s ${START} -e ${END} -df ${STEP} -r ${CUTOFF} -solv ${SOLVENT} -ex ${EXCLUDE} -hb > file.log"
echo ${EXE2}
$(eval ${EXE2} )

