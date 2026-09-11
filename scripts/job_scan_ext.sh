#!/bin/bash
# PRP49 extended scan: 13 peptides (incl. 3 homologs) x 10 targets
#SBATCH --job-name=prp49_scanx
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=00-01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -e
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1

CKPT=${1:?usage: sbatch job_scan_ext.sh <ckpt>}

source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate prp49

cd $HOME/LassoPep/PRP49
mkdir -p ../results
python -m prp49.scan --config config.yaml --ckpt ${CKPT} \
    --peptides ../mvp_cpu/scan_peptides_ext.csv \
    --targets ../mvp_cpu/scan_targets.fasta \
    --device cuda
mv ../results/scan_matrix.csv ../results/scan_matrix_ext.csv
echo "SCAN_EXT_DONE"
