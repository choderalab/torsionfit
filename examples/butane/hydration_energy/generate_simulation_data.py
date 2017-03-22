from torsionfit.hydration_energy import energy
from parmed.charmm import CharmmParameterSet
try:
    import cPickle as pickle
except:
    import pickle

# Load CHARMM parameters files
params = CharmmParameterSet('../../data/charmm_ff/top_all36_cgenff.rtf', '../../data/charmm_ff/par_all36_cgenff.prm',
                           '../../data/charmm_ff/toppar_water_ions.str')

db = energy.create_openmm_system('butane', '../structure', params)
db_solv = energy.generate_simulation_data(db, params, n_steps=4000, n_iter=500)

pickle.dump(db_solv, open('solvated.pickle', 'w'))
