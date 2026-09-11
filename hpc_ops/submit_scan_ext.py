"""Upload extended scan inputs and submit scan job."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()
put_file(cli, r'D:\deepseek_harness\prp49\mvp_cpu\scan_peptides_ext.csv',
         f'{home}/LassoPep/mvp_cpu/scan_peptides_ext.csv', progress=False)
put_file(cli, r'D:\deepseek_harness\prp49\hpc_deploy\scripts\job_scan_ext.sh',
         f'{home}/LassoPep/scripts/job_scan_ext.sh', progress=False)
st, out, err = run(cli, 'cd ~/LassoPep && chmod +x scripts/job_scan_ext.sh; '
                        'sbatch scripts/job_scan_ext.sh runs/checkpoints/fold0_best.pt')
print('scan_ext:', out.strip(), err.strip())
st, out, err = run(cli, 'squeue -u $USER -o "%.10i %.14j %.8T %.10M" | head -8')
print(out)
cli.close()
