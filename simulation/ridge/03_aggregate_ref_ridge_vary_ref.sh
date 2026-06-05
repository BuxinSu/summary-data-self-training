#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/submit_aggregate_vary_ref.out

mkdir -p job_script

job_script="job_script/job_aggregate_vary_ref.sh"

cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=01:00:00
#SBATCH --mem=8g
#SBATCH --nodes=1
#SBATCH --job-name=aggregate_vary_ref
#SBATCH --output=/path/to/summary_training/logs/output_aggregate_vary_ref.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/03_aggregate_ref_ridge_vary_ref.py
EOT

sbatch "$job_script"
