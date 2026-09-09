#!/bin/bash
# ============================================================
# PRP49 L2 fine-tuning: Siyuan-1 a100 queue, 1 GPU
# cpus-per-task=16 is the official a100 pairing (16 cores per GPU)
# ============================================================
#SBATCH --job-name=prp49_train
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=2-00:00:00
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
python -c "import torch, transformers; print('torch', torch.__version__, torch.cuda.is_available(), 'tf', transformers.__version__)"

python -m prp49.train --config config.yaml --device cuda
echo "TRAIN_DONE"
