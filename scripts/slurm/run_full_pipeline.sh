#!/bin/bash
#SBATCH --job-name=phenocycler_pipeline
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@ufl.edu  # UPDATE: Use your UFL email
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=96:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --qos=your_qos        # UPDATE: Find with 'sacctmgr show assoc user=$USER'
#SBATCH --account=your_account  # UPDATE: Find with 'sacctmgr show assoc user=$USER'

# Full Phenocycler Analysis Pipeline

# Configurable timeouts (in seconds)
TIMEOUT_PREPROCESS=86400   # 24 hours
TIMEOUT_PHENOTYPE=86400    # 24 hours
TIMEOUT_SPATIAL=172800     # 48 hours
TIMEOUT_GROUPCOMP=43200    # 12 hours
TIMEOUT_TISSUE=86400       # 24 hours

echo "=========================================="
echo "Phenocycler Analysis Pipeline - Full Run"
echo "=========================================="
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo ""

# Load modules
module load conda

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate phenocycler_analysis

# Change to project directory
cd $SLURM_SUBMIT_DIR/../..

# Function to run notebook and check success
run_notebook() {
    local notebook=$1
    local timeout=$2
    
    echo ""
    echo "=========================================="
    echo "Running: $notebook"
    echo "Timeout: $timeout seconds"
    echo "Started at: $(date)"
    echo "=========================================="
    
    jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=$timeout \
        --output-dir=data/processed \
        notebooks/$notebook
    
    if [ $? -eq 0 ]; then
        echo "✓ $notebook completed successfully"
        return 0
    else
        echo "✗ $notebook failed!"
        return 1
    fi
}

# Run pipeline sequentially
run_notebook "01_preprocessing.ipynb" $TIMEOUT_PREPROCESS || exit 1
run_notebook "02_phenotyping.ipynb" $TIMEOUT_PHENOTYPE || exit 1
run_notebook "03_spatial_analysis.ipynb" $TIMEOUT_SPATIAL || exit 1
run_notebook "04_group_comparisons.ipynb" $TIMEOUT_GROUPCOMP || exit 1
run_notebook "05_tissue_comparisons.ipynb" $TIMEOUT_TISSUE || exit 1

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "Finished at: $(date)"
echo "=========================================="
