for T in 1 2
do
SLURM="#!/bin/bash\n\

#SBATCH -t 100:00:00\n\
#SBATCH --job-name=julia_batch_${T}\n\
#SBATCH -p normal\n\
#SBATCH --export=ALL\n\
#SBATCH --exclusive\n\
#SBATCH --nodes=1\n\              
#SBATCH --output=julia_batch_${T}.out\n\
#SBATCH --error=julia_batch${T}.err\n\
#SBATCH --ntasks-per-node=16\n\
   
   
module load julia\n\
julia Model_Ensemble.jl $T"

echo -e $SLURM | sbatch
sleep 0.5
done
