#!/bin/bash
# ============================================================
# PRP49 multi-seed parallel: --array=1-4%4 = 4 independent 1-GPU jobs
# Each task builds its own config (seed / output dir)
# ============================================================
#SBATCH --job-name=prp49_seed
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
#SBATCH --array=1-4%4
#SBATCH --output=%x-%A_%a.out
#SBATCH --error=%x-%A_%a.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49

PRP49_ROOT=$HOME/LassoPep
cd $PRP49_ROOT/PRP49

SEED=$((40 + SLURM_ARRAY_TASK_ID))
CONF=config_seed${SEED}.yaml
sed -e "s/^  seed: .*/  seed: ${SEED}/" \
    -e "s|^  output_dir: .*|  output_dir: ./runs_seed${SEED}|" \
    -e "s|^  checkpoint_dir: .*|  checkpoint_dir: ./runs_seed${SEED}/checkpoints|" \
    config.yaml > ${CONF}

echo "=== task ${SLURM_ARRAY_TASK_ID} seed=${SEED} conf=${CONF} ==="
python -m prp49.train --config ${CONF} --device cuda
echo "SEED_DONE ${SEED}"
