"""Upload docking directory (pep/rec/tasks + job script) to server."""
import os
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME')
home = out.strip()
sftp = cli.open_sftp()

local_dir = r'D:\deepseek_harness\prp49\docking'
for sub in ['pep', 'rec']:
    run(cli, f'mkdir -p {home}/LassoPep/docking/{sub}')
    for fn in os.listdir(os.path.join(local_dir, sub)):
        sftp.put(os.path.join(local_dir, sub, fn), f'{home}/LassoPep/docking/{sub}/{fn}')
        print('up', sub, fn)
sftp.put(os.path.join(local_dir, 'tasks.csv'), f'{home}/LassoPep/docking/tasks.csv')
print('up tasks.csv')
sftp.close()
put_file(cli, r'D:\deepseek_harness\prp49\hpc_deploy\scripts\job_dock_array.sh',
         f'{home}/LassoPep/scripts/job_dock_array.sh', progress=False)
st, out, err = run(cli, 'ls ~/LassoPep/docking/pep ~/LassoPep/docking/rec | head -25')
print(out)
cli.close()
print('DOCK_UPLOAD DONE')
