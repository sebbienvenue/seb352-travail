#!/bin/bash

#SBATCH --no-requeue
#SBATCH --job-name="aiida-None"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00


module purge
module load intel/15.0.0
module load intelmpi/5.0.1
module load fftw/3.3.4/intel-15.0.0





module load intel/16.0.3
module load intelmpi/5.1.3
module load fftw/3.3.4-mpi
module load mkl/11.3.3

'mpirun' '-np' '8' '/home/spbienve/scratch/executables/pw_new.x' '-in' 'aiida.in'  > 'aiida.out' 
