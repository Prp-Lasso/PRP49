"""Download scan matrix."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run

cli = connect()
sftp = cli.open_sftp()
sftp.get('/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/results/scan_matrix.csv',
         r'D:\deepseek_harness\prp49\results\scan_matrix.csv')
sftp.get('/dssg/home/acct-clswxl/clswxl-ccmbi1/LassoPep/PRP49/runs/checkpoints/fold0_best.pt',
         r'D:\deepseek_harness\prp49\mvp_cpu\fold0_best_650M.pt')
print('downloaded scan_matrix.csv + fold0 ckpt')
sftp.close()
cli.close()
