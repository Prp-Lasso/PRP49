#!/bin/bash
# ============================================================
# PRP49 S1 warmup: Bernett 20k pairs, a100, 1 GPU (~1-2 h)
# Both ESM-2 35M encoders frozen; only BAN + head train.
# Output: runs_warmup/checkpoints for later curriculum-learning init
# ============================================================
#SBATCH --job-name=prp49_warmup
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00-03:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49

PRP49_ROOT=$HOME/LassoPep
cd $PRP49_ROOT/PRP49

echo "=== node: $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv
echo "=== pairs: $(wc -l ../mvp_cpu/bernett_pairs_20k.csv) ==="

python -m prp49.train --config config_warmup.yaml --device cuda
echo "WARMUP_DONE"
