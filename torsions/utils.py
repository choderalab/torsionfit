
"""
"""

from chemistry.charmm import CharmmPsfFile
from chemistry.charmm.parameters import CharmmParameterSet
import simtk.openmm as mm
import simtk.openmm.app as app
from cclib.parser import Gaussian
import re
import mdtraj as md
import numpy as np
from copy import deepcopy

def to_optimize(param, stream, penalty = 10):
    ''' returns a list of dihedrals to optimize and updates CharmmParameterSet
    with stream files
    
    Parameters
    ----------
    param : CharmmParameterSet
    stream: list of stream files 
    penalty: int for CGenFF penalty cutoff (Default = 10)
    
    Returns list of tuples containing dihedrals to optimize
    
    '''
    keys = [i for i in param.dihedral_types.iterkeys()]
    for j in stream:
        param.read_stream_file(j)
    return [k for k in param.dihedral_types.iterkeys() 
    if k not in keys and param.dihedral_types[k].penalty >= penalty]
        


def compute_energy_from_positions(param, structure, positions, platform=None):
    '''
    Computes energy for a given structure with a given parameter set
    
    Parameters
    ----------
    param: chemistry.charmm.CharmmParameterSet
    structure: chemistry.structure 
    positions: simtk.unit.Quantity wrapping [natoms,ndim] numpy array
    platform: simtk.openmm.Platform to evaluate energy on (if None, will select automatically)
    '''
    integrator = mm.VerletIntegrator(0.004)
    system = structure.createSystem(param, )
    if platform != None:
        context = mm.Context(system, integrator, platform)
    else:
        context = mm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy = True)
    del context
    return state.getPotentialEnergy()
    
# def parse_Glog(logFile, top):
#     '''
#     Extracts coordinates and energy from a list of Gaussian 09 log files
#
#     Paramerers
#     ----------
#     logFile:  list of str file names of Gaussian09 log files
#     top : {str, Trajectory, Topology} Pass in either the path to a RCSB PDB file,
#     a trajectory with  one frame,  or a topology to supply this information
#
#     Returns
#     -------
#     traj: mdtraj.Trajectory with all orientations in g09 log files
#     energy: np.array energies corresponding to each frame in traj
#     steps: list of lists [step_number, scan_point, total_scan_points ]
#     step_number is the step in geom opt of that scan point
#     '''
#     traj = md.load(top)
#     steps = []
#     index = (2,12,-1)
#     log = Gaussian(logFile)
#     data = log.parse()
#     for fi in logFile:
#         f = open(fi, 'r')
#         for line in f:
#             if re.search('Step', line):
#                 try:
#                     steps.append([line.rsplit()[j] for j in index])
#                 except:
#                     pass
#         f.close
#     # convert angstroms to nm (mdtraj Trajectory coordinates are in nm, Gaussian
#     # are in angstroms)
#     coords = data.atomcoords*0.1
#     t = deepcopy(traj)
#     for k in range(len(coords)):
#         if k == 0:
#             traj.xyz = coords[0]
#         else:
#             t.xyz = coords[k]
#             traj = traj.join(t)
#     return (traj, data.scfenergies, steps)
#
# def extract_opt(trajectory, energy, steps):
#     '''
#     extract optimzed geomoetry and corresponding energies
#
#     Parameters
#     ----------
#     trajectory: mdtraj.Trajectory
#     energy: numpy.array with energies correspoding to positions in mdtraj.Trajectory
#     steps: list of lists [n_step, scan_point, total_scan_points]
#
#     Returns
#     -------
#     traj: mdtraj.Trajectory with geometry optimized orientations
#     energy_opt: numpy.array of energies corresponding to geometry optimized orientations
#     steps: list of lists of [step, scan_point, total_scan_points]
#     '''
#     positions = np.ndarray((0,trajectory.n_atoms,3))
#     energy_opt = np.ndarray((0))
#     points = []
#     traj = trajectory[0]
#     # check if all are of the same size
#     try:
#         assert(trajectory.n_frames == len(energy) == len(steps))
#     except AssertionError:
#         print 'energy and steps must be the same length as trajectory.n_frames'
#
#     for i, j in enumerate(steps):
#         try:
#             if j[1] != steps[i+1][1]:
#                 positions = np.append(positions, trajectory[i].xyz, axis = 0)
#                 energy_opt = np.append(energy_opt, energy[i])
#                 points.append(steps[i])
#         except IndexError:
#             positions = np.append(positions, trajectory[i].xyz, axis = 0)
#             energy_opt = np.append(energy_opt, energy[i])
#             points.append(steps[i])
#
#     t = deepcopy(traj)
#     for k in range(len(positions)):
#         if k == 0:
#             traj.xyz = positions[0]
#         else:
#             t.xyz = positions[k]
#             traj = traj.join(t)
#
#     return (traj, energy_opt, points)
#
# def get_g09_filename(filename):
#     return filename.split('/')[-1]
#
# def get_frag_name(filename):
#     f = filename.split('/')[-1]
#     return f.split('.')[0]
#
# def get_torsion(filename):
#     f = open(filename, 'r')
#     for line in f:
#         if re.search('   Scan   ', line):
#             torsion = line.split()[2].split(',')
#             torsion[0] = torsion[0][-1]
#             torsion[-1] = torsion[-1][0]
#             for i in range(len(torsion)):
#                 torsion[i] = (int(torsion[i]) - 1)
#     f.close()
#     return torsion
#

                
     
        
        
    
       
    
            
        

