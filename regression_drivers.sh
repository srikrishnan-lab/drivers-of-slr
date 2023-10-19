#!/bin/bash
#SBATCH -J slr_drivers
#SBATCH -p normal
#SBATCH -t 10:00:00
#SBATCH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH -n 1
#SBATCH -o scen_disc.out
#SBATCH -e scen_disc.err

echo "starting at `date` on `hostname`"

module purge

julia +1.9 src/regression_drivers.jl

echo "ended at `date` on `hostname`"
exit 0
