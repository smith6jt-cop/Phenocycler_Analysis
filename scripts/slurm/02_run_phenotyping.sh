#!/bin/bash
#SBATCH --job-name=phenocycler_phenotype
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu  # UPDATE: Use your UFL email
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/phenotype_%j.out
#SBATCH --error=logs/phenotype_%j.err
#SBATCH --qos=your_qos        # UPDATE: Find with 'sacctmgr show assoc user=$USER'
#SBATCH --account=your_account  # UPDATE: Find with 'sacctmgr show assoc user=$USER'

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"

module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phenocycler_analysis

cd $SLURM_SUBMIT_DIR/../..

echo "Running phenotyping notebook..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=86400 \
    --output-dir=data/processed \
    notebooks/02_phenotyping.ipynb

echo "Job completed at: $(date)"
