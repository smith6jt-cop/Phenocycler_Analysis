#!/bin/bash
#SBATCH --job-name=phenocycler_preprocess
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu  # UPDATE: Use your UFL email
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/preprocess_%j.out
#SBATCH --error=logs/preprocess_%j.err
#SBATCH --qos=your_qos        # UPDATE: Find with 'sacctmgr show assoc user=$USER'
#SBATCH --account=your_account  # UPDATE: Find with 'sacctmgr show assoc user=$USER'

# Phenocycler Preprocessing - HiPerGator SLURM Script

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"

# Load modules
module load conda

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phenocycler_analysis

# Change to project directory
cd $SLURM_SUBMIT_DIR/../..

# Run notebook
echo "Running preprocessing notebook..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir=data/processed \
    notebooks/01_preprocessing.ipynb

echo "Job completed at: $(date)"
