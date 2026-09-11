#!/bin/bash
# ============================================================
# PRP49 docking array: 8 lasso peptides x 10 human targets + control
# Vina 1.2.3 (autodock_vina env), 1 CPU per task, 81 tasks
# Receptor pdbqt: obabel (add H + Gasteiger); ligand: meeko (rigid macrocycle)
# ============================================================
#SBATCH --job-name=prp49_dock
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00-01:30:00
#SBATCH --array=1-81%16
#SBATCH --output=dock_%a.out
#SBATCH --error=dock_%a.err

set -e
source /usr/share/lmod/lmod/init/profile
module load miniconda3/4.10.3
source activate autodock_vina

cd $HOME/LassoPep/docking
mkdir -p pdbqt out

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tasks.csv)
IFS=',' read -r TID PEP REC CX CY CZ SX SY SZ <<< "$LINE"

REC_PDBQT="pdbqt/${REC}.pdbqt"
PEP_PDBQT="pdbqt/${PEP}.pdbqt"

# receptor: obabel add H + charges, then strip ROOT/BRANCH tags (vina rigid format)
if [ ! -f "$REC_PDBQT" ]; then
  obabel -ipdb "rec/${REC}.pdb" -opdbqt -O "pdbqt/${REC}_full.pdbqt" -p 7.4 2>>"out/${TID}.prep.err"
  python3 rec_clean.py "pdbqt/${REC}_full.pdbqt" "$REC_PDBQT" 2>>"out/${TID}.prep.err"
fi
# ligand: obabel pdbqt (add H + charges), then rigidify (lasso peptide as rigid body)
if [ ! -f "$PEP_PDBQT" ]; then
  obabel -ipdb "pep/${PEP}.pdb" -opdbqt -O "pdbqt/${PEP}_full.pdbqt" -p 7.4 2>>"out/${TID}.prep.err"
  python3 rigidify.py "pdbqt/${PEP}_full.pdbqt" "$PEP_PDBQT" 2>>"out/${TID}.prep.err"
fi

vina --receptor "$REC_PDBQT" --ligand "$PEP_PDBQT" \
     --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
     --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
     --exhaustiveness 8 --num_modes 5 --cpu 1 \
     --out "out/${TID}_${PEP}_${REC}.pdbqt" 2>&1 | tail -3

grep "REMARK VINA RESULT" "out/${TID}_${PEP}_${REC}.pdbqt" | head -1 > "out/${TID}.score"
cat "out/${TID}.score"
echo "DOCK_DONE $TID $PEP $REC"
