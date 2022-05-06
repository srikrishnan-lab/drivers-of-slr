#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=julia_test
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1               
#SBATCH --output=julia_test.out
#SBATCH --error=julia_test.err
#SBATCH --ntasks-per-node=16


module load julia

julia Model_Ensemble.jl  
