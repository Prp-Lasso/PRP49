"""Template for the HPC remote-ops helpers used in this project.

Credentials are read from environment variables - never hard-code them.
    setx SJTU_HPC_USER clswxl-xxxxxx      (or export on Linux)
    setx SJTU_HPC_PASSWORD ...
"""
import os
import paramiko

HOST = os.environ.get('SJTU_HPC_HOST', 'sylogin.hpc.sjtu.edu.cn')
USER = os.environ['SJTU_HPC_USER']
PASSWORD = os.environ['SJTU_HPC_PASSWORD']
BASE = 'LassoPep'   # operate strictly inside the user's own directory


def connect(timeout=30):
    cli = paramiko.SSHClient()
    cli.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    cli.connect(HOST, port=22, username=USER, password=PASSWORD,
                timeout=timeout, banner_timeout=60, auth_timeout=60,
                look_for_keys=False, allow_agent=False)
    return cli


def run(cli, cmd, timeout=600, get_pty=False):
    stdin, stdout, stderr = cli.exec_command(cmd, timeout=timeout, get_pty=get_pty)
    out = stdout.read().decode('utf-8', 'replace')
    err = stderr.read().decode('utf-8', 'replace')
    return stdout.channel.recv_exit_status(), out, err
