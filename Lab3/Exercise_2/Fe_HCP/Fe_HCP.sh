#!/bin/bash

echo "a a/c V E" > result-hcp-opti.dat
PW_LAUNCH="/home/sebastien/Programs/espresso/espresso-5.3.0/bin/pw.x"
PREFIX="Fe_HCP"
spanvol=0.5
span=0.01
ca0=1.58
vol0=138

for j in {-15..15}
do
for i in {0..0}	#{-15..15}
do

#vol=`echo "$spanvol * $i + $vol0" | bc -l`
c_over_a=`echo "$span * $j + $ca0" |  bc -l`
#a=`echo "$vol $c_over_a" | awk '{printf "%.4f", (2*$1/sqrt(3)/$2)^(1.0/3)}'`
a=4.6548
vol=`echo "0.5 * ($a)^3 * sqrt(3) * $c_over_a" | bc -l`

INFILE="$PREFIX.scf.$vol.$c_over_a.$a.in"
OUTFILE="$PREFIX.scf.$vol.$c_over_a.$a.out"
rm -f $INFILE
rm -f $OUTFILE

cat > $INFILE << EOF
 &control
    calculation = 'scf'
    restart_mode='from_scratch'
    prefix='$PREFIX'
    tstress = .true.
    tprnfor = .true.
    outdir = './temp/'
    pseudo_dir = '../PP/'
 /        
 &system    
    ibrav=  4
    celldm(1) = $a
    celldm(3) = $c_over_a
    nat=  2
    ntyp= 1
    ecutwfc = 55.0
    ecutrho = 240.0
    occupations = 'smearing'
    degauss = 0.03
    smearing = 'm-v'
 /
 &electrons
    mixing_beta = 0.7 
    conv_thr =  1.0d-8
 /
 ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-n-rrkjus_psl.0.2.4.UPF
 ATOMIC_POSITIONS crystal
  Fe   0.3333333333  0.6666666667  0.25 
  Fe   0.6666666667  0.3333333333  0.75
 K_POINTS automatic
    12 12 6 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE

E=`grep ! $OUTFILE | awk '{print $5}'`
echo "$a $c_over_a $vol $E" >> result-hcp-opti.dat

done
done
rm -r *.out *.in temp/*
echo "DONE"
