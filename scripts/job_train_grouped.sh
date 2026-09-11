#!/bin/bash
# PRP49 leakage-safe verification: identical to improved run but grouped CV by peptide
#SBATCH --job-name=prp49_grp
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

cd $HOME/LassoPep/PRP49
echo "=== node: $(hostname) ==="
nvidia-smi --query-gpu=name,memory.total --format=csv

python -m prp49.train --config config_improved_grouped.yaml --device cuda \
    --init_ckpt "runs_warmup/checkpoints/fold{fold}_best.pt"
echo "TRAIN_GROUPED_DONE"
