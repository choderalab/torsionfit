__author__ = 'Chaya D. Stern'

from celery import Celery
import yaml
import os
from shutilwhich import which
from fnmatch import fnmatch
from celery import group

config_file = open(os.environ['CELERY_CONFIG'], 'r')
config = yaml.load(config_file)
config_file.close()

app = Celery('psi4_jobs',
             broker=config['broker'],
             backend=config['backend'],
             include=['torsionfit.qmscan.qmtasks'])


@app.task
def start_psi4_calculation(path, input_file):
    print(input_file)
    output_file = input_file.replace('dat', 'out')
    input_file = os.path.join(path, input_file)
    output_file = os.path.join(path, output_file)
    psi4_binary = which('psi4', mode=os.X_OK)
    cmd = psi4_binary + ' ' + input_file + ' -o ' + output_file + '2>&1'
    return os.system(cmd)


def run_psi4_distributed(directory):
    to_submit = []
    pattern = "*.dat"
    for path, subdir, files in os.walk(directory):
        for name in files:
            if fnmatch(name, pattern):
                path = os.path.join(os.getcwd(), path)
                to_submit.append((path, name))
    exits = group(start_psi4_calculation.apply_async(args=[path, input_file]) for path, input_file in to_submit)()
    return exits