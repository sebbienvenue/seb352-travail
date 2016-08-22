#!/bin/bash
clear 

echo Starting at `date`


#BIN_DIR=/home/spbienve/dft_U/dfptu_git_03.07.2016/bin
BIN_DIR=/home/bienvenue/Documents/DFPT+U/dfptu_git_03.07.2016/bin
PSEUDO_DIR=.
TMP_DIR=./tmp

#clean the tmp/
rm -rf $TMP_DIR/* 
rm input_tmp.in

#################################
# pw.x
#################################
BIN_LIST="pw.x"

PSEUDO_LIST="co_pbesol_v1.2.uspp.F.UPF O.pbesol-n-rrkjus_psl.0.1.UPF li_pbesol_v1.4.uspp.F.UPF"

# Self-consistent calculation
cat > ./inputs/LiCoO2.scf.in << EOF
 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='LiCoO2',
    pseudo_dir = '$PSEUDO_DIR/',
    outdir='$TMP_DIR/'
    verbosity='high'
 /
 &system
    ibrav = 0,
    celldm(1) = 5.6565
    nat = 4, 
    ntyp = 3,
    nspin = 1,
    ecutwfc = 80.0,
    ecutrho = 800.0,
    lda_plus_u = .true.
    lda_plus_u_kind = 0,
    U_projection_type = 'atomic',
    Hubbard_U(1) = 8.9288
 /
 &electrons
    diagonalization='david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-12
 /
ATOMIC_SPECIES
 Co  58.933194 co_pbesol_v1.2.uspp.F.UPF
 O   15.999    O.pbesol-n-rrkjus_psl.0.1.UPF
 Li   6.94     li_pbesol_v1.4.uspp.F.UPF
ATOMIC_POSITIONS {crystal}
 Co       0.000002777  -0.000009040   0.000000001
 Li       0.500019952   0.499924316   0.500000034
 O        0.260990277   0.260969573   0.217073943
 O        0.739015074   0.739013121   0.782926071
CELL_PARAMETERS
  0.469726941   0.813605960   0.001548470
 -0.469743434   0.813594472   0.001572713
  0.000017760   0.539401256   1.538453479
K_POINTS {automatic}
 8 8 6 0 0 0
EOF

PW_COMMAND="mpirun -np 4 $BIN_DIR/pw.x"


$PW_COMMAND < ./inputs/LiCoO2.scf.in > ./outputs/LiCoO2.scf.out


#################################
# ph.x
#################################



echo finished at `date`
