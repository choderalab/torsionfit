#!/bin/tcsh
#  Batch script for single thread CPU job.
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=46:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
#PBS -k oe
#
# specify queue
#PBS -q batch
#
# nodes: number of nodes
#   ppn: number of processes per node
#PBS -l nodes=1:ppn=1
#PBS -l mem=100g
#
# export all my environment variables to the job
#PBS -V
#
# job name (default = name of script file)
#PBS -N 360_all_mult_off_ch_b_an_off
#
# mail settings (one or more characters)
# email is sent to local user, unless another email address is specified with PBS -M option 
# n: do not send mail
# a: send mail if job is aborted
# b: send mail when job begins execution
# e: send mail when job terminates
#PBS -M chaya.stern@choderalab.org
#PBS -m abe
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
#PBS -o /cbio/jclab/home/chayas/src/ChayaSt/torsionfit/examples/butane/360/all_torsion/eliminate_phase_mult_off/charge_off/bond_angles_off 

# Change to working directory used for job submission
cd $PBS_O_WORKDIR

# Launch my program.
python /cbio/jclab/home/chayas/src/ChayaSt/torsionfit/examples/butane/360/all_torsion/eliminate_phase_mult_off/charge_off/bond_angles_off/sampling.py
