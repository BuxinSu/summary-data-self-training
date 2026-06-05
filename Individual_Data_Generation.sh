#!/bin/bash
# rm -rf Individual_Data
# rm -rf job_script
# rm -rf output
# rm -rf error

# mkdir -p Individual_Data
# mkdir -p job_script

kappa_values=(0.05 0.5)
noise_std_values=(0.5 4)
p_values=(5000)

for p in "${p_values[@]}"; do
    for kappa in "${kappa_values[@]}"; do
        for noise_std in "${noise_std_values[@]}"; do
            # Create a unique job script for each combination
            job_script="job_script/job_p${p}_kappa${kappa}_noise${noise_std}.sh"

            # Ensure output and error directories exist
            mkdir -p output/output_p${p}_kappa${kappa}_noise${noise_std}
            mkdir -p error/error_p${p}_kappa${kappa}_noise${noise_std}

            cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH -A standby
#SBATCH --job-name=sim_p${p}_kappa${kappa}_noise${noise_std}
#SBATCH --output=output/output_p${p}_kappa${kappa}_noise${noise_std}/output_p${p}_kappa${kappa}_noise${noise_std}.txt
#SBATCH --error=error/error_p${p}_kappa${kappa}_noise${noise_std}/error_p${p}_kappa${kappa}_noise${noise_std}.txt
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem=64g

module load anaconda

# Activate virtual environment if needed
source /path/to/summary_training/myenv/bin/activate

python Individual_Data_Generation.py --p ${p} --kappa ${kappa} --noise_std ${noise_std} --output_dir Individual_Data/ --seed 42
EOT

            # Submit the job script
            sbatch "$job_script"
        done
    done
done
