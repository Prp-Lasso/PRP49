"""Upload bundle to ~/LassoPep/ (SFTP, resumable not needed if stable)."""
import sys
sys.path.insert(0, r'D:\deepseek_harness\prp49\hpc_deploy\ops')
from opslib import connect, run, put_file

cli = connect()
st, out, err = run(cli, 'echo $HOME; ls -la ~/LassoPep | head -5')
home = out.strip().splitlines()[0]
print('home:', home)
remote_path = f'{home}/LassoPep/prp49_bundle.tar.gz'
put_file(cli, r'D:\deepseek_harness\prp49\hpc_deploy\prp49_bundle.tar.gz', remote_path)
st, out, err = run(cli, f'ls -la {remote_path}; sha256sum {remote_path}')
print(out)
cli.close()
print('UPLOAD DONE')
