#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/Individual_Data_Generation_heavy_tail.out

mkdir -p Individual_Data_heavy_tail
mkdir -p job_script

kappa_values=(0.05 0.5 0.9)
noise_std_values=(0.5 2.0 3.0)
p=5000
n_values=(7500 12500)
df_values=(3 5)

for df in "${df_values[@]}"; do
    for n in "${n_values[@]}"; do
        for kappa in "${kappa_values[@]}"; do
            for noise_std in "${noise_std_values[@]}"; do
                # Generate a deterministic seed from (n, df, kappa, noise_std)
                seed=$(echo "${n}_${df}_${kappa}_${noise_std}" | cksum | awk '{print $1 % 100000}')

                # Create a unique job script for each combination
                job_script="job_script/job_individual_ht_n${n}_df${df}_kappa${kappa}_noise${noise_std}.sh"

                cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=19G
#SBATCH --nodes=1
#SBATCH --job-name=ht_n${n}_df${df}_k${kappa}_ns${noise_std}
#SBATCH --output=/path/to/summary_training/logs/output_ht_n${n}_df${df}_kappa${kappa}_noise${noise_std}.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/01_Individual_Data_Generation_heavy_tail.py --p ${p} --n ${n} --df ${df} --kappa ${kappa} --noise_std ${noise_std} --output_dir Individual_Data_heavy_tail/ --seed ${seed}
EOT

                # Submit the job script
                sbatch "$job_script"
            done
        done
    done
done
