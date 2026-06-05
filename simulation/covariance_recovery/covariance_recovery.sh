#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/submit_cov_recovery.out

mkdir -p /path/to/summary_training/results/covariance_recovery
mkdir -p job_script

p_values=(1000 2000 5000)

for p in "${p_values[@]}"; do
    job_script="job_script/job_cov_recovery_p${p}.sh"

    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --job-name=cov_rec_p${p}
#SBATCH --output=/path/to/summary_training/logs/output_cov_recovery_p${p}.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/covariance_recovery.py --p ${p} --rho 0.6 --n_reps 50 --output_dir /path/to/summary_training/results/covariance_recovery/ --seed 42
EOT

    sbatch "$job_script"
done
