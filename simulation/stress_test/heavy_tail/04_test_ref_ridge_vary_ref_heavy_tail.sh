#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/submit_test_vary_ref_heavy_tail.out

mkdir -p job_script

job_script="job_script/job_test_vary_ref_heavy_tail.sh"

cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=02:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
#SBATCH --job-name=test_vary_ref_heavy_tail
#SBATCH --output=/path/to/summary_training/logs/output_test_vary_ref_heavy_tail.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/04_test_ref_ridge_vary_ref_heavy_tail.py
EOT

sbatch "$job_script"
