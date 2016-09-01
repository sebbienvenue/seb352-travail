#!/bin/sh
clear 
echo "Starting at `date`"

go=8

echo "running for Si ..."


PSEUDO_DIR=.
PSEUDO="Si.pbesol-n-rrkjus_psl.1.0.0.UPF"
TMP_DIR=./tmp
BIN="/home/bienvenue/Documents/DFPT+U/dfptu_git_03.07.2016/bin"
RUN="mpirun -np $go $BIN"

#inputs
ECUT="30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 110.0 120.0 130.0 140.0"
RATIO="3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0"
LISTK="4"
LISTA="10.20"
echo "a    k    ecut    ratio    energy" > results.dat
echo " running the scf calculation for Si..."

for a in $LISTA
do
for ratio in $RATIO
do 
#echo "a    k    ecut    ratio    energy" > results.$ratio.dat
for j in $ECUT
do
for k in $LISTK
do 

rm -rf TMP_DIR/*

#ecutrho
RHO=`expr "$ratio * $j" | bc -l`

#self-consistent calculation
cat > ./inputs/si.scf.$a.$k.$j.$ratio.in << EOF
 &control
    calculation='scf',
    restart_mode='from_scratch',
    prefix='si'
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    tstress = .true.
    tprnfor = .true.
 /
 &system
    ibrav = 2, celldm(1) =$a, nat=  2, ntyp= 1,
    ecutwfc = $j,
    ecutrho = $RHO,
 /
 &electrons
    mixing_beta = 0.7
    diagonalization = 'david'
    mixing_mode = 'plain'
    conv_thr =  1.0d-8
 /
ATOMIC_SPECIES
 Si  28.086  Si.pbesol-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25
K_POINTS {automatic}
 $k $k $k 0 0 0
EOF

echo " running a=$a k=$k ecut=$j ratio=$ratio ..."
$RUN/pw.x < ./inputs/si.scf.$a.$k.$j.$ratio.in > ./outputs/si.scf.$a.$k.$j.$ratio.out

energy=`grep "!    total energy" ./outputs/si.scf.$a.$k.$j.$ratio.out | awk '{print $5}'`
echo "$a $k $j $ratio $energy" >> results.$ratio.dat

done
done
done
done

echo "Ended at `date`"
