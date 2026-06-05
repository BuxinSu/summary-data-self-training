#!/bin/bash

rm -rf results_reference_ridge
rm -rf job_script
rm -rf output
rm -rf error

mkdir -p results_reference_ridge
mkdir -p job_script

kappa_values=(0.05 0.5 0.9)
noise_std_values=(0.5 2 4)
n_values=(100 500 2500 5000)

for iteration in {1..20}; do
    for kappa in "${kappa_values[@]}"; do
        for noise_std in "${noise_std_values[@]}"; do
            for n in "${n_values[@]}"; do
                # Create a unique job script for each combination
                job_script="job_script/job_kappa${kappa}_noise${noise_std}_n${n}_iter${iteration}.sh"

                # Ensure output and error directories exist
                mkdir -p output/output_kappa${kappa}_noise${noise_std}_n${n}
                mkdir -p error/error_kappa${kappa}_noise${noise_std}_n${n}

                cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH -A standby
#SBATCH --job-name=sim_kappa${kappa}_noise${noise_std}_n${n}_iter${iteration}
#SBATCH --output=output/output_kappa${kappa}_noise${noise_std}_n${n}/output_kappa${kappa}_noise${noise_std}_n${n}_iter${iteration}.txt
#SBATCH --error=error/error_kappa${kappa}_noise${noise_std}_n${n}/error_kappa${kappa}_noise${noise_std}_n${n}_iter${iteration}.txt
#SBATCH --nodes=1 
#SBATCH --time=04:00:00
#SBATCH --mem=64g

module load anaconda

# Activate virtual environment if needed
source /path/to/summary_training/myenv/bin/activate

python summery_training_ref_ridge_12_22.py --kappa ${kappa} --noise_std ${noise_std} --n ${n} --iteration ${iteration}
EOT

                # Submit the job script
                sbatch "$job_script"
            done
        done
    done
done
