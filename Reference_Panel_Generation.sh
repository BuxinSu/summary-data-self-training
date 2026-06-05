#!/bin/bash

# Directories for results and job scripts
mkdir -p reference_panels
mkdir -p job_scripts

# Parameters
p_values=(5000 10000)

for p in "${p_values[@]}"; do
    # Create a unique job script for each p value
    job_script="job_scripts/job_p${p}.sh"

    # Ensure output and error directories exist
    mkdir -p output/output_p${p}
    mkdir -p error/error_p${p}

    cat <<EOT > "$job_script"
#!/bin/bash
#SBATCH -A standby
#SBATCH --job-name=sim_p${p}
#SBATCH --output=output/output_p${p}/output_p${p}.txt
#SBATCH --error=error/error_p${p}/error_p${p}.txt
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem=40g
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=ALL

module load anaconda

# Activate virtual environment if needed
source /path/to/summary_training/myenv/bin/activate

python Reference_Panel_Generation.py --p ${p} --rho 0.6 --output_dir reference_panels/
EOT

    # Submit the job script
    sbatch "$job_script"
done
