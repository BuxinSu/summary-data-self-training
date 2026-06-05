#!/bin/bash

#SBATCH -A your_account
#SBATCH --time=04:00:00
#SBATCH --job-name=summery-step1
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=ALL
#SBATCH --mem=30g                # 50GB memory per job
#SBATCH --nodes=1                # One node per job
#SBATCH --array=1-71            # Submit 186 jobs in parallel (trait IDs 1 to 186)

# Load the R module (if necessary)
module load r/4.4.1

# Run the R script with the array task ID
Rscript single-ancestry-step1-DXA.R --trait=${SLURM_ARRAY_TASK_ID}
