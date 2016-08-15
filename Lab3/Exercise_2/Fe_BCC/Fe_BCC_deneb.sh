#!/bin/bash -f
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --cpus-per-task 1
#SBATCH --time=00:30:00
#SBATCH --mem=14000
#SBATCH --partition=mse-468
#SBATCH --account=mse-468
#SBATCH --reservation=mse-468
 
#These are the libraries necessary for our code to run
module purge
module load  intel/15.0.0
module load  intelmpi/5.0.1
module load  fftw/3.3.4/intel-15.0.0

cd /scratch/spbienve/fcc/Fe_BCC/
echo "lattice Energy" > results-bcc.dat

##### QE-Compilers ######
PW_LAUNCH="srun /scratch/spbienve/executables/pw.x"

# Input data
#
# List of values for the lattice paramter
PREFIX="Fe_BCC"
a0=5.4
span=0.01

# Loop over lattice constants
for i in {-10..10}
do
a=`echo "$a0 + $i * $span" | bc -l`

INFILE="$PREFIX.scf.$a.in"
OUTFILE="$PREFIX.scf.$a.out"
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
    ibrav=  3
    celldm(1) = $a
    nat=  1
    ntyp= 1
    ecutwfc = 55.0
    ecutrho = 240.0
    starting_magnetization(1) = 0.7
    occupations = 'smearing'
    degauss = 0.03
    smearing = 'm-v'
    nspin = 2
 /
 &electrons
    mixing_beta = 0.7 
    conv_thr =  1.0d-8
 /
 ATOMIC_SPECIES
  Fe  55.845  Fe.pbe-n-rrkjus_psl.0.2.4.UPF
 ATOMIC_POSITIONS crystal
  Fe 0.00 0.00 0.00 
 K_POINTS automatic
    12 12 12 0 0 0 
EOF

$PW_LAUNCH < $INFILE > $OUTFILE
E=`grep ! $OUTFILE | awk '{print $5}'`
echo "$a $E" >> results-bcc.dat

done
rm -r *.save *.in *.out temp/*
