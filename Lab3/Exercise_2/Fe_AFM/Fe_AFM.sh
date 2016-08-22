#!/bin/bash -f

clear
PW_LAUNCH="/home/sebastien/Programs/espresso/espresso-5.3.0/bin/pw.x"

a0=5.27033
span=0.01

# Prefix name
PREFIX="Fe_AFM"
echo "a E" > result-afm.dat
# Loop over lattice constants
for i in {-20..20}
do

a=`echo "$a0 + $i * $span" | bc -l`

INFILE="$PREFIX.scf.$a.in"
OUTFILE="$PREFIX.scf.$a.out"

cat > $INFILE << EOF
 &control
    calculation = 'scf'
    restart_mode='from_scratch'
    prefix='$PREFIX'
    tstress = .true.
    tprnfor = .true.
    outdir = './temp1/'
    pseudo_dir = '../PP/'
 /        
 &system    
    ibrav=  1
    celldm(1) = $a
    nat=  2
    ntyp= 2
    ecutwfc = 55
    ecutrho = 240
    occupations = 'smearing'
    degauss = 0.03
    smearing = 'm-v'
    nspin = 2
    starting_magnetization(1) =  0.7
    starting_magnetization(2) =  -0.7
 /
 &electrons
    mixing_beta = 0.7 
    conv_thr =  1.0d-8
 /
 ATOMIC_SPECIES
  FeU  55.845  Fe.pbe-n-rrkjus_psl.0.2.4.UPF
  FeD  55.845  Fe.pbe-n-rrkjus_psl.0.2.4.UPF

 ATOMIC_POSITIONS crystal
  FeU  0.0 0.0 0.0
  FeD  0.5 0.5 0.5
 K_POINTS automatic
    12 12 12 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE
E=`grep ! $OUTFILE | awk '{print $5}'`
echo "$a $E" >> result-afm.dat
done
rm -r *.save *.in *.out temp1/*