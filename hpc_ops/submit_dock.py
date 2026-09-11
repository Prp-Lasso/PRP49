"""Upload rigidify.py + updated job script; test rigid ligand prep; submit array."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()

put_file(cli, r'D:\deepseek_harness\prp49\docking\rigidify.py',
         f'{home}/LassoPep/docking/rigidify.py', progress=False)
put_file(cli, r'D:\deepseek_harness\prp49\hpc_deploy\scripts\job_dock_array.sh',
         f'{home}/LassoPep/scripts/job_dock_array.sh', progress=False)

# test on login node (file conversion only, no docking)
st, out, err = run(cli, 'cd ~/LassoPep/docking && '
                        'source /usr/share/lmod/lmod/init/profile; module load miniconda3/4.10.3; '
                        'source activate autodock_vina; '
                        'obabel -ipdb pep/MccJ25.pdb -opdbqt -O pdbqt/MccJ25_full.pdbqt -p 7.4 2>&1 | tail -1; '
                        'python3 rigidify.py pdbqt/MccJ25_full.pdbqt pdbqt/MccJ25.pdbqt; '
                        'head -2 pdbqt/MccJ25.pdbqt; tail -2 pdbqt/MccJ25.pdbqt')
print(out)
print(err[-300:] if err else '')

# submit array
st, out, err = run(cli, 'cd ~/LassoPep && sbatch scripts/job_dock_array.sh; '
                        'squeue -u $USER -o "%.10i %.14j %.8T %.10M"')
print('submit:', out, err)
cli.close()
print('DOCK_SUBMITTED')
