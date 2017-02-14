import mdtraj as md
import os
from fnmatch import fnmatch


def pdb_to_psi4(pdb, mol_name, method, basis_set, charge=0, multiplicity=1, symmetry=None, geom_opt=True,
                sp_energy=False, fixed_dih=None, mem=None):
    """

    :param pdb: str
        path to pdb file
    :param method: list of str
        QM method (see psi4 website for options)
        If length 2, first one will be used for geom opt and second for spe.
    :param basis_set: str
        specification of basis set
    :param symmetry: str
        symmetry of molecule. Default is None.
    :param geom_opt: bool
        if True, will generate input file for geometry optimization
    :param sp_energy: bool
        if True, will run a single point energy calculation (if geom_opt also true, SPE calculation will occur after
        geom opt
    :param fixed_dih: str
        string of dihedral that should be fixed at specified angle. Format: "4 7 10 14 90.00"
        default: None - will not fix dihedral
        Beware:
        ------
        Because of a bug in psi4, dihedral angle can't be exactly 0 (same would apply for 180) so use 0.001 instead

    :param mem: int
        memory allocation for calculation
    :param outfile: str
        if specified, will save file there
    :return:
        psi4 input string. If outfile, save file to specified path
    """

    input_string = ""

    if mem is not None:
        input_string += "\nmemory {}\n".format(mem)

    input_string += "\nmolecule {}".format(mol_name)
    input_string += " {\n"
    if symmetry is not None:
        input_string += "  symmetry {}\n".format(symmetry)
    input_string += "  {} {} \n".format(charge, multiplicity)

    mol = md.load(pdb)
    for i, atom in enumerate(mol.topology.atoms):
        element = atom.element.symbol
        # Convert to Angstroms
        xyz = mol.xyz[0]*10
        input_string += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(element, xyz[i][0], xyz[i][1], xyz[i][2])

    input_string += "  units Angstrom\n"
    input_string += "}\n"

    if fixed_dih is not None:
        input_string += '\ndih_string = "{}"'.format(fixed_dih)
        # ToDo add string because that's the only thing that seems to work
        input_string += '\nset optking fixed_dihedral = $dih_string\n'

    if geom_opt:
        input_string += "\noptimize('{}/{}')\n".format(method[0], basis_set[0])

    if sp_energy:
        input_string += "\nenergy('{}/{}')\n".format(method[-1], basis_set[-1])

    return input_string


def generate_scan_input(root, filetype, mol_name, dihedral, method, basis_set, charge=0, multiplicity=1, symmetry=None,
                        geom_opt=True, sp_energy=False, mem=None):
    """
    This function takes a directory and writes out psi4 input files for all files that match the filetype specified

    :param root: str
        path to files
    :param filetype: str
        input filetypes
    :param mol_name: str
        molecule name
    :param dihedral: str
        index of atoms that should remain fixed. format '1  2  3  4'
    :param method: list of str
        QM method (see psi4 website for options)
    :param basis_set: list of str
        see psi4 website for options
    :param charge: int
        default 0
    :param multiplicity: int
        default 1
    :param symmetry: str
        symmetry of molecule. default None
    :param geom_opt: bool
        if True, run geometry optimization
    :param sp_energy: bool
        if True, run a single point energy calculation after geomoetry optimization
    :param mem: str
        memory allocation

    """

    input_files = []
    pattern = "*.{}".format(filetype)
    for path, subdir, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                input_files.append(os.path.join(path, name))

    for f in input_files:
        fixed_dih_angle = f.split('/')[-2]
        if fixed_dih_angle == '0':
            fixed_dih_angle = '0.001'
        if fixed_dih_angle == '180':
            fixed_dih_angle = '180.001'
        dihedral_string = dihedral + ' ' + fixed_dih_angle
        output = pdb_to_psi4(pdb=f, mol_name=mol_name, method=method, basis_set=basis_set, charge=charge,
                             multiplicity=multiplicity, symmetry=symmetry, geom_opt=geom_opt, sp_energy=sp_energy,
                             fixed_dih=dihedral_string, mem=mem)

        filename = f.replace(filetype, 'dat')
        psi4_input = open(filename, 'w')
        psi4_input.write(output)
        psi4_input.close()




