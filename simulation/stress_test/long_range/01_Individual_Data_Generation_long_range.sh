#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/Individual_Data_Generation_long_range.out

mkdir -p Individual_Data_long_range
mkdir -p job_script

kappa_values=(0.05 0.5 0.9)
noise_std_values=(0.5 2.0 3.0)
p=5000
n_values=(7500 12500)
rank=5
alpha=0.1

for n in "${n_values[@]}"; do
    for kappa in "${kappa_values[@]}"; do
        for noise_std in "${noise_std_values[@]}"; do
            # Generate a deterministic seed from (n, kappa, noise_std)
            seed=$(echo "${n}_${kappa}_${noise_std}" | cksum | awk '{print $1 % 100000}')

            # Create a unique job script for each combination
            job_script="job_script/job_individual_lr_n${n}_kappa${kappa}_noise${noise_std}.sh"

            cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=19G
#SBATCH --nodes=1
#SBATCH --job-name=lr_n${n}_k${kappa}_ns${noise_std}
#SBATCH --output=/path/to/summary_training/logs/output_lr_n${n}_kappa${kappa}_noise${noise_std}.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/01_Individual_Data_Generation_long_range.py --p ${p} --n ${n} --rank ${rank} --alpha ${alpha} --kappa ${kappa} --noise_std ${noise_std} --output_dir Individual_Data_long_range/ --seed ${seed}
EOT

            # Submit the job script
            sbatch "$job_script"
        done
    done
done
