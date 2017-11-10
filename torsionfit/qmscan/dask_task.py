from dask import delayed, compute
import os
from shutilwhich import which
from fnmatch import fnmatch
import subprocess

# To run a dask for distributed qm
# 1) write a python script that calls run_psi4_distributed. This script should call a client with the scheduler.json file
# 2) Write a qsub PBS script that calls the dask scheduler and writes a scheduler.json file
# 3) Write a bash script that qsubs the PBS script with array job using -t 0-N%m (N - number of jobs, m, how many jobs
# at a time) and then calls the python script that that has the client and python jobs.
def start_psi4_calculation(path, input_file, threads=4):
    print(input_file)
    output_file = input_file.replace('dat', 'out')
    input_file = os.path.join(path, input_file)
    output_file = os.path.join(path, output_file)
    psi4_binary = which('psi4', mode=os.X_OK)
    cmd = psi4_binary + ' ' + input_file + ' -o ' + output_file + ' -n' + threads
    process = subprocess.Popen(cmd, stderr=subprocess.PIPE, shell=True)
    output = process.communicate()
    return(cmd, process.wait(), output)

def run_psi4_distributed(directory, threads):
    to_submit = []
    pattern = "*.dat"
    for path, subdir, files in os.walk(directory):
        for name in files:
            if fnmatch(name, pattern):
                path = os.path.join(os.getcwd(), path)
                to_submit.append((path, name, threads))
    delayed_tasks = [delayed(start_psi4_calculation)(path, f, threads) for (path, f, threads) in to_submit]
    output = compute(delayed_tasks)
    return output
