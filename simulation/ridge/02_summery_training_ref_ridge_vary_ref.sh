#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/submit_ref_ridge_vary_ref.out

mkdir -p /path/to/summary_training/results/ref_ridge
mkdir -p job_script

kappa_values=(0.05 0.5 0.9)
noise_std_values=(0.5 2.0 3.0)
p=5000
n_values=(7500 12500)
n_w_values=(200 1000 20000)

for iteration in {1..2}; do
    for n in "${n_values[@]}"; do
        for n_w in "${n_w_values[@]}"; do
            for kappa in "${kappa_values[@]}"; do
                for noise_std in "${noise_std_values[@]}"; do
                    # Create a unique job script for each combination
                    job_script="job_script/job_vary_ref_kappa${kappa}_noise${noise_std}_p${p}_n${n}_nw${n_w}_iter${iteration}.sh"

                    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=16g
#SBATCH --nodes=1
#SBATCH --job-name=vref_kappa${kappa}_noise${noise_std}_p${p}_n${n}_nw${n_w}_iter${iteration}
#SBATCH --output=/path/to/summary_training/logs/output_vary_ref_kappa${kappa}_noise${noise_std}_p${p}_n${n}_nw${n_w}_iter${iteration}.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/02_summery_training_ref_ridge_vary_ref.py --kappa ${kappa} --noise_std ${noise_std} --iteration ${iteration} --p ${p} --n ${n} --n_w ${n_w}
EOT

                    # Submit the job script
                    sbatch "$job_script"
                done
            done
        done
    done
done
