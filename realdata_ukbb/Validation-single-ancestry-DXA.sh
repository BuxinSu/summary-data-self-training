#!/bin/bash

#SBATCH -A your_account
#SBATCH --time=1-00:00:00
#SBATCH --job-name=pumas-validation
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=ALL
#SBATCH --mem=50g                # 50GB memory per job
#SBATCH --nodes=1                # One node per job

# Load the R module (if necessary)
module load r/4.1.2

# Run the R script with the array task ID
Rscript Validation-single-ancestry-DXA.R
