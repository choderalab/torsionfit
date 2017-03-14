from torsionfit.qmscan.torsion_scan import generate_scan_input
import torsionfit.qmscan.qmtasks as qmtasks

# To generate structures for dihedral scan, use generate_dihedral.py in pymol
# To use generate_dihedral in pymol, open pymol.
# import script into PyMOL
# run /path_to_torsionfit/torsionfit/qm_scan/generate_dihedral.py
# load butane.pdb
# call torsion_drive function
# torsion_drive C1, C2, C3, C4, 10, butane, path_to_torsion_scan_inputs, butane

# Generate psi4 input files
generate_scan_input('B3LYP_torsion_scan/', filetype='pdb', mol_name='butane', dihedral='4 7 10 14',
                    method=['B3LYP'], basis_set=['cc-PVDZ'], symmetry='C1', mem='1000 mb')

# Submit tasks
# set environment variable CELERY_CONFIG which points to the config file (when running on the cluster the config file
# has to be edited for the location of the redis server
# First open redis by calling redis-server
# Start celey workers with this command:
# celery -A qmtasks worker -l info (Run this command in directory where qmtasks is (torsionfit/qmscan). The default
# concurrency is 8 which means that 8 jobs will run at a time. This option can be changed).
# Then run this script
qmtasks.run_psi4_distributed('B3LYP_torsion_scan/')
