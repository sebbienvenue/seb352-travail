#!/bin/bash
clear
rm -r out/ 2>/dev/null

echo "code exo 4"

PATH_BIN2=/home/bienvenue/Programs/espresso-5.3.0/bin
PATH_BIN=/home/bienvenue/Documents/DFPT+U/dfptu_git_03.07.2016/bin
go="mpirun -np 6 "

echo "Running pw.x ..."

$go$PATH_BIN/pw.x < AlAs.scf.in > AlAs.scf.out

echo "Finished pw.x"
echo
echo "Running ph.x ..." 
$go$PATH_BIN/ph.x < AlAs.ph.in > AlAs.ph.out
echo "Finished ph.x"
echo 
echo "Running q2r.x ..." 
$go$PATH_BIN/q2r.x < AlAs.q2r.in > AlAs.q2r.out
echo "Finished q2r.x"
echo 
echo "Running matdyn.x ..." 
$go$PATH_BIN/matdyn.x < AlAs.matdyn.in > AlAs.matdyn.out
echo "matdyn.x finished"
echo "Running plotband.x ..." 
#$PATH_BIN/band_plot.x < AlAs.plotband.in > AlAs.plotband.out
$PATH_BIN/plotband.x < AlAs.plotband.in > AlAs.plotband.out
#$PATH_BIN2/plotband.x < AlAs.plotband.in > AlAs.plotband.out
echo
echo "plotband.x finished"
echo 
echo "FINISHED CALCULATIONS"

