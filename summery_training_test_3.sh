#!/bin/bash

#SBATCH -A standby
#SBATCH --job-name=summery_training
#SBATCH --output=/path/to/summary_training/summery_n_%A_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36  # Adjust this number based on how many CPUs you want to allocate per task
#SBATCH --nodes=1
#SBATCH --mail-type=end
#SBATCH --mail-user=your_email@example.com
#SBATCH --time=4:00:00
#SBATCH --mem=70g

module load anaconda
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
