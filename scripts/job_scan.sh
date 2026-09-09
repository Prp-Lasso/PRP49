#!/bin/bash
# ============================================================
# PRP49 10-peptide x 10-target scan: a100, 1 GPU
# usage: sbatch job_scan.sh <runs/checkpoints/fold0_best.pt>
# ============================================================
#SBATCH --job-name=prp49_scan
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00-02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

CKPT=${1:?usage: sbatch job_scan.sh <runs/checkpoints/fold0_best.pt>}

source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49

cd $HOME/LassoPep/PRP49
mkdir -p ../results
python -m prp49.scan --config config.yaml --ckpt ${CKPT} \
    --peptides ../mvp_cpu/scan_peptides.csv \
    --targets ../mvp_cpu/scan_targets.fasta \
    --device cuda
echo "SCAN_DONE"
