#!/bin/bash
#SBATCH --job-name=phenocycler_spatial
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_%j.out
#SBATCH --error=logs/spatial_%j.err
#SBATCH --qos=your_qos
#SBATCH --account=your_account

echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"

module load conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phenocycler_analysis

cd $SLURM_SUBMIT_DIR/../..

echo "Running spatial analysis notebook..."
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=172800 \
    --output-dir=data/processed \
    notebooks/03_spatial_analysis.ipynb

echo "Job completed at: $(date)"
