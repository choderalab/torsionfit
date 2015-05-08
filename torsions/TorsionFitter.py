# -*- coding: utf-8 -*-
"""
Created on Thu May  7 19:32:22 2015

@author: sternc1
"""
from collections import OrderedDict
import pandas as pd
from cclib.parser import Gaussian
import re
import mdtraj as md

#def load_scan(filenames, top):
#    """load Gaussian torsio scan log files from disk
#    
#    Parameters
#    ---------
#    filename: list of str
#    top: str. Pass in the path to a RCB PDB file
#    
#    Return
#    ------
#    torsionscanset: torsion.TorsionScanSet
#    """
#    traj = md.load(top)
#    for fi in filenames:
#        filename = get_g09_filename(fi)
#        fragment = get_frag_name(fi)
#        torsion = get_torsion(fi)
#        steps = get_steps(fi)
#        log = Gaussian(fi)
#        data = log.parse()
        
        
    
    
def get_g09_filename(filename):
    return filename.split('/')[-1]

def get_frag_name(filename):
    f = filename.split('/')[-1]
    return f.split('.')[0]

def get_torsion(filename):
    f = open(filename, 'r')
    for line in f:
        if re.search('   Scan   ', line):
            torsion = line.split()[2].split(',')
            torsion[0] = torsion[0][-1]
            torsion[-1] = torsion[-1][0]
            for i in range(len(torsion)):
                torsion[i] = (int(torsion[i]) - 1)
    f.close()
    return torsion
    
def get_steps(filename):
    steps = []
    index = (2,12,-1)
    f = open(filename, 'r')
    for line in f:
        if re.search('Step', line):
            try:
                steps.append([line.rsplit()[j] for j in index])
            except:
                pass
    f.close
    return steps
        

class TorsionScanSet(object):
    """container object for torsion scan
    
    A TorsionScanSet should be constructed by loading Gaussian log files from the 
    scan with a topoloyg - usually a PDB file
    
    Attributes
    ----------
    
    positions: dict((str, str, str, str) : unit.Quantity(np.array[natoms, ndim], u.nm))
    Dictionary mapping the 4 element tuple of mol_name, torsion, directions, and step
    to the positions
    torsion_set: pandas data frame. 
    """
    
    def __init__(self):
        """Create new TorsionScanSet object"""
        self._positions = OrderedDict()
        self._torsion_set = pd.DataFrame()
        
        
    def extract_opt(self):
        """Returns TorsionScanSet object with only optimized geometries"""
        
    def compute_energy(self, frag,  param, platform=None):
        """computes energy for a given structure with a given paramere set"""
        
        
        
    
    
    
    