#!/bin/sh

#SBATCH -n 12
#SBATCH --time=72:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=flu8@illinois.edu


module purge
module use /data/cigi/common/cigi-modules
module add GNU610
module add GPU
module load anaconda2

python -u run.py $1 $2 15
