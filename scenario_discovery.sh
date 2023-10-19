#!/bin/bash
#SBATCH -J SLR_SCEN_DISC
#SBATCH -p normal
#SBATCH -t 10:00:00
#SBATCH --exclusive
#SBATCH -n 1
#SBATCH -o scen_disc.out
#SBATCH -e scen_disc.err

echo "starting at `date` on `hostname`"

module purge
julia

julia src/scenario_discovery.jl

echo "ended at `date` on `hostname`"
exit 0