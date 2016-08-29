#!/bin/bash
clear 
echo Starting at `date`

go=1

#BIN_DIR=/home/sebastien/Programs/dftU/dfptu_git_03.07.2016/bin
BIN_DIR=/home/bienvenue/Documents/DFPT+U/dfptu_git_03.07.2016/bin
PSEUDO_DIR=.
TMP_DIR=./tmp



#################################
# ph.x
#################################
BIN_LIST="ph.x"
PH_COMMAND="mpirun -np $go $BIN_DIR/ph.x"
PH_COMMAND="$BIN_DIR/ph.x"


# Phonon calculation at Gamma
cat > ./inputs/LiCoO2.ph.in << EOF
phonons of LiCoO2
 &inputph
  prefix = 'LiCoO2',
  outdir='$TMP_DIR/',
  tr2_ph = 1.0d-14,
  amass(1) = 58.933194,
  amass(2) = 15.999,
  amass(3) =  6.94,
  ldisp = .true.,
  nq1=4, nq2=4, nq3=4
  fildyn='LiCoO2.dyn',
 /
EOF



$PH_COMMAND < ./inputs/LiCoO2.ph.in > ./outputs/LiCoO2.ph.out





echo finished at `date`
