"""Submit model scan (10 peptides x 10 targets) on a100 with fold0 checkpoint."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

cli = connect()
st, out, err = run(cli, 'ls ~/LassoPep/PRP49/runs/checkpoints/ 2>/dev/null; '
                        'ls ~/LassoPep/mvp_cpu/scan_peptides.csv ~/LassoPep/mvp_cpu/scan_targets.fasta 2>&1')
print('=== assets ===')
print(out)
st, out, err = run(cli, 'cd ~/LassoPep && sbatch scripts/job_scan.sh runs/checkpoints/fold0_best.pt')
print('scan submit:', out.strip(), err.strip())
st, out, err = run(cli, 'squeue -u $USER -o "%.10i %.14j %.8T %.10M" | head -6')
print(out)
cli.close()
