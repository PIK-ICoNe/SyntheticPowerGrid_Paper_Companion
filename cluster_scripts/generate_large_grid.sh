#!/bin/bash

#SBATCH --qos=medium
#SBATCH --partition=standard
#SBATCH --account=icone
#SBATCH --error=%x-%j-%N.err
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --tasks-per-node=10
#SBATCH --job-name=generate_large_grid
#SBATCH --mem=64GB

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.8.5
julia --threads=$SLURM_NTASKS generate_large_grid.jl
