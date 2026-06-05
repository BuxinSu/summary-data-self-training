#!/bin/bash
rm -rf results_reference_ridge_01_09
rm -rf job_script
rm -rf output
rm -rf error

mkdir -p results_reference_ridge_01_09
mkdir -p job_script

kappa_values=(0.05 0.5 0.9)
noise_std_values=(0.5 2 4)
p_values=(5000 10000)

for iteration in {1..20}; do
    for kappa in "${kappa_values[@]}"; do
        for noise_std in "${noise_std_values[@]}"; do
            for p in "${p_values[@]}"; do
                # Create a unique job script for each combination
                job_script="job_script/job_kappa${kappa}_noise${noise_std}_p${p}_iter${iteration}.sh"

                # Ensure output and error directories exist
                mkdir -p output/output_kappa${kappa}_noise${noise_std}_p${p}
                mkdir -p error/error_kappa${kappa}_noise${noise_std}_p${p}

                cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH -A standby
#SBATCH --job-name=sim_kappa${kappa}_noise${noise_std}_p${p}_iter${iteration}
#SBATCH --output=output/output_kappa${kappa}_noise${noise_std}_p${p}/output_kappa${kappa}_noise${noise_std}_p${p}_iter${iteration}.txt
#SBATCH --error=error/error_kappa${kappa}_noise${noise_std}_p${p}/error_kappa${kappa}_noise${noise_std}_p${p}_iter${iteration}.txt
#SBATCH --nodes=1 
#SBATCH --time=04:00:00
#SBATCH --mem=16g

module load anaconda

# Activate virtual environment if needed
source /path/to/summary_training/myenv/bin/activate

python summery_training_ref_ridge_01_09.py --kappa ${kappa} --noise_std ${noise_std} --iteration ${iteration} --p ${p}
EOT

                # Submit the job script
                sbatch "$job_script"
            done
        done
    done
done
