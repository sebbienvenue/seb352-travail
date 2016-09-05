#!/bin/sh
clear 
echo Starting at `date`
echo "running for post processing Si grid"
go=8

PSEUDO_DIR=.
PSEUDO="Si.pbesol-n-rrkjus_psl.1.0.0.UPF"
TMP_DIR=./tmp
BIN="/home/sebastien/Programs/dftU/dfptu_git_03.07.2016/bin"
#RUN="mpirun -np $go $BIN"


echo 
echo "Running q2r.x ..." 
$BIN/q2r.x < si.q2r.in > si.q2r.out
echo "Finished q2r.x"
echo 
echo "Running matdyn.x ..." 
$BIN/matdyn.x < si.matdyn.in > si.matdyn.out
echo "matdyn.x finished"
echo "Running plotband.x ..." 
$BIN/plotband.x < si.plotband.in > si.plotband.out
echo
echo "plotband.x finished"
echo 
echo 
echo plot with gnuplot
gnuplot plot_dispersion.gnu
echo
echo "FINISHED CALCULATIONS"


