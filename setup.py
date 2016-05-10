""" torsionfit: A toolkit for Bayesian torsion parameterization for molecular mechanics forcefields.

This project provides tools for using Bayesian parameterization techniques and Markov chain Monte Carlo
to fit molecular mechanics torsion Fourier terms to quantum chemical data.
"""

from __future__ import print_function

import os
from os.path import relpath, join

import numpy
import versioneer

from setuptools import setup, Extension, find_packages

DOCLINES = __doc__.split("\n")

########################
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)
Programming Language :: Python
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

################################################################################
# USEFUL SUBROUTINES
################################################################################
def find_package_data(data_root, package_root):
    files = []
    for root, dirnames, filenames in os.walk(data_root):
        for fn in filenames:
            files.append(relpath(join(root, fn), package_root))
    return files


################################################################################
# SETUP
################################################################################

setup(
    name='torsionfit',
    author='Chaya Stern',
    author_email='chaya.stern@choderalab.org',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='LGPL',
    url='https://github.com/choderalab/torsions',
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=['torsionfit', 'torsionfit.tests'],
    package_data={ 'torsionfit.tests': ['reference/*']},
    zip_safe=False,
    install_requires=[
        'numpy',
        'pymc',
        'pandas',
        'cclib',
        'openmm>=6.3',
        'parmed',
        'netcdf4'
        ],
    )

