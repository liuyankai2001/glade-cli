import os
import signal
import subprocess
import sys
from pathlib import Path


def worker_command(job_dir):
    return [sys.executable, '-B', '-m', 'app.worker', str(job_dir)]


def worker_environment():
    env = os.environ.copy()
    # app is directly under /opt/bactabolize in the image. Also support
    # importing the service by its repository namespace on development hosts.
    package_root = str(Path(__file__).resolve().parents[1])
    env['PYTHONPATH'] = package_root + os.pathsep + env.get('PYTHONPATH','')
    return env


def terminate_group(process, grace=1):
    if os.name == 'nt':
        subprocess.run(['taskkill','/PID',str(process.pid),'/T','/F'],capture_output=True)
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            pass
        # The group leader can exit before children. Always kill the group.
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    process.wait(timeout=10)
