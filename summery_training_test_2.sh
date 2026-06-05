#!/bin/bash

#SBATCH -A standby
#SBATCH --job-name=summery_training     # Job name
#SBATCH --output=/path/to/summary_training/summery_n_%A_%a.out # Updated output file with job array ID
#SBATCH --ntasks=1                   # Number of tasks (1 task = 1 core)
#SBATCH --nodes=1                    # Number of nodes (1 node)
#SBATCH --mail-type=end
#SBATCH --mail-user=your_email@example.com
#SBATCH --time=4:00:00
#SBATCH --mem 70g

module load anaconda 

# Activate virtual environment if needed
source ~/venv/bin/activate

# Define n_values and kappa_values
n_values=(5000 10000 50000)
kappa_values=(0.1 0.2 0.5)

# Job array index (SLURM_ARRAY_TASK_ID)
index_n=$(($SLURM_ARRAY_TASK_ID / ${#kappa_values[@]}))
index_k=$(($SLURM_ARRAY_TASK_ID % ${#kappa_values[@]}))

n=${n_values[$index_n]}
kappa=${kappa_values[$index_k]}

# Run the Python script with the current n and kappa values
python summery_train_test_2.py --n_value $n --kappa_value $kappa
