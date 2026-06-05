#!/bin/bash
#SBATCH -A your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --output=/path/to/summary_training/logs/Reference_Panel_Generation.out

mkdir -p reference_panels
mkdir -p job_script

# Parameters
p_values=(5000)
n_w_values=(200 1000 20000)

for p in "${p_values[@]}"; do
    for n_w in "${n_w_values[@]}"; do
        # Seed depends only on (p, n_w) — same W across all chunks
        seed=$(echo "${p}_${n_w}" | cksum | awk '{print $1 % 100000}')

        for chunk_index in $(seq 0 9); do
            # Create a unique job script for each (p, n_w, chunk_index)
            job_script="job_script/job_reference_p${p}_nw${n_w}_chunk${chunk_index}.sh"

            cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH --account=your_account
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --job-name=ref_p${p}_nw${n_w}_c${chunk_index}
#SBATCH --output=/path/to/summary_training/logs/output_p${p}_nw${n_w}_chunk${chunk_index}.txt
#SBATCH --chdir=/path/to/summary_training

module purge
module load conda
source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate /path/to/summary_training/myenv

python /path/to/summary_training/code/Reference_Panel_Generation.py --p ${p} --n_w ${n_w} --rho 0.6 --output_dir reference_panels/ --seed ${seed} --chunk_index ${chunk_index}
EOT

            # Submit the job script
            sbatch "$job_script"
        done
    done
done
