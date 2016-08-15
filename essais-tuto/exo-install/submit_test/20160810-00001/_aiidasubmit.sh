#!/bin/bash

#SBATCH --no-requeue
#SBATCH --job-name="aiida-None"
#SBATCH --get-user-env
#SBATCH --output=_scheduler-stdout.txt
#SBATCH --error=_scheduler-stderr.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00


'mpirun' '-np' '1' '/home/sebastien/Programs/espresso/espresso-5.3.0/bin/pw.x' '-in' 'aiida.in'  > 'aiida.out' 
