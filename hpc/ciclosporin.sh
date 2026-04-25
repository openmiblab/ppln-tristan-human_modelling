#!/bin/bash
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --time=95:00:00
#SBATCH --mail-user=s.sourbron@sheffield.ac.uk
#SBATCH --mail-type=FAIL,END
#SBATCH --comment=ciclosporin
#SBATCH --job-name=ciclosporin
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

unset SLURM_CPU_BIND
export SLURM_EXPORT_ENV=ALL

module load Anaconda3/2024.02-1
module load Python/3.10.8-GCCcore-12.2.0

# Initialize Conda for this non-interactive shell
eval "$(conda shell.bash hook)"

conda activate tristan

# Get the current username
USERNAME=$(whoami)
BASE_DIR="/mnt/parscratch/users/$USERNAME/scripts/tristan-human-stage-2-modelling"
DATA="$BASE_DIR/data"
BUILD="$BASE_DIR/wip/build"

# Safer to use explicit path of python installation
PYTHON="/users/md1spsx/.conda/envs/tristan/bin/python"

# Run scripts
srun $PYTHON "$BASE_DIR/wip/tristan_ciclosporin_effect.py" --data="$DATA" --build="$BUILD" --scans=2
srun $PYTHON "$BASE_DIR/wip/tristan_ciclosporin_effect.py" --data="$DATA" --build="$BUILD" --scans=1