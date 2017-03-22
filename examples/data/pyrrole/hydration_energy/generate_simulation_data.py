from torsionfit.hydration_energy import energy
from parmed.charmm import CharmmParameterSet
try:
    import cPickle as pickle
except:
    import pickle

# Load CHARMM parameters files
params = CharmmParameterSet('../../charmm_ff/top_all36_cgenff.rtf', '../../charmm_ff/par_all36_cgenff.prm',
                           '../../charmm_ff/toppar_water_ions.str', '../pyrrol.str')

db = energy.create_openmm_system('pyrrole', '../structure', params)
db_solv = energy.generate_simulation_data(db, params, n_steps=3000, n_iter=400)

# remove system and context because pickle is giving me problems pickling them

db_new = {}
for key in db_solv['solvated'].keys():
    if key == 'system' or key == 'context':
	continue
    db_new[key] = db_solv['solvated'][key]

pickle.dump(db_new, open('solvated.pickle', 'w'))
pickle.dump(db_solv, open('solvated_test.pickle', 'w'))

