#!/bin/sh
clear 
echo Starting at `date`

go=8

echo running for Si
echo Starting at `date`

ECUT="60 70 80 90 100 110 120 130 140"
LISTK="6 8 10 12 16"
LISTA="10.20"
echo "a,    k,    ecut,    energy" > results.dat
echo " running the scf calculation for Si..."

for a in $LISTA
do
for j in $ECUT
do 
for k in $LISTK
do 

PSEUDO_DIR=.
PSEUDO="Si.pbesol-n-rrkjus_psl.1.0.0.UPF"
TMP_DIR=./tmp
BIN="/home/bienvenue/Documents/DFPT+U/dfptu_git_03.07.2016/bin"
RUN="mpirun -np $go $BIN"

rm -rf TMP_DIR/*

#self-consistent calculation
cat > ./inputs/si.scf.$a.$k.$j.in << EOF
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
    ecutwfc = 120.0,
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

echo " running a=$a k=$k ecut=$j ..."
$RUN/pw.x < ./inputs/si.scf.$a.$k.$j.in > ./outputs/si.scf.$a.$k.$j.out

energy=`grep "!    total energy" ./outputs/si.scf.$a.$k.$j.out | awk '{print $5}'`
echo "$a, $k, $j, $energy" >> results.dat

done
done
done

echo Ended at `date`
