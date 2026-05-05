#!/bin/bash -l
#SBATCH --job-name=make_gif
#SBATCH --account=m3739
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=cpu
#SBATCH --cpus-per-task=128
#SBATCH --qos=debug
#SBATCH --time=00:30:00

# Initialize conda for compute node
module load conda
conda activate pygkyl

# Run Flan run
srun --cpu-bind=cores python nt_pt_plots.py

