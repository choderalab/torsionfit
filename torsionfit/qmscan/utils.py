from openeye import oechem, oeiupac, oedepict
import os
import numpy as np
from torsionfit.utils import logger
import time
import json
import psi4

def write_oedatabase(moldb, ofs, mlist, size):
    """
    This function writes out a new oedatabase from an existing database

    Parameters
    ----------
    moldb: OpenEye database object
    ofs: output stream
    mlist: list of indices of molecules in moldb for new database
    size: size of new database

    """
    for molidx in mlist[:size:]:
        moldb.WriteMolecule(ofs, molidx)


def to_smi(smiles, path, base, return_fname=False):
    """
    This function writes out an .smi file for a list of SMILES
    Parameters
    ----------
    smiles: list of SMILES.
        The list can also contain strings that include name for SMILES separated by a space. ("SMILES Name")
    path: str
        path to output file
    base: str
        base name for output file

    """
    fname = os.path.join(path, base + '.smi')
    outf = open(fname, 'w')
    smiles_list = map(lambda x: x+"\n", list(smiles))
    outf.writelines(smiles_list)
    outf.close()
    if return_fname:
        return fname


def create_oedatabase_idxfile(ifname):
    """
    This function creates an index file associated with a given molecule filename. It write out the file in the same
    directory the parent molecule file is and adds an .idx extension to parent molecule filename.

    From OpenEye's documentation:
    The index file of a molecule dtabase stores file position offsets of the molecules in the file. Generating an index
    file can be expensive, but it can be created only once and then it can speed up the handling of large molecule files
    significantly

    Parameters
    ----------
    ifname: str
        absolute path to molecule file
    """
    idx_fname = oechem.OEGetMolDatabaseIdxFileName(ifname)

    if os.path.exists(idx_fname):
        oechem.OEThrow.Warning("{} index file already exists".format(idx_fname))
    elif not oechem.OECreateMolDatabaseIdx(ifname):
        oechem.OEThrow.Warning("Unable to create {} molecule index file".format(idx_fname))


def new_output_stream(outname):
    """
    This function creates a new oechem.oemolostream.
    Parameters
    ----------
    outname: str
        name of outputfile.

    Returns
    -------
    ofs: oechem.oemolostream

    """
    ofs = oechem.oemolostream()
    if not ofs.open(outname):
        oechem.OEThrow.Fatal("Unable to open {} for writing".format(outname))
    return ofs


def to_oemol(filename, title=None):
    """Create OEMol from file. If more than one mol in file, return list of OEMols.

    Parameters
    ----------
    filename: str
        absolute path to
    title: str
        title for molecule. If None, IUPAC name will be given as title.

    Returns
    -------
    mollist: list
        list of OEMol for multiple molecules. OEMol if file only has one molecule.
    """
    ifs = oechem.oemolistream(filename)
    moldb = oechem.OEMolDatabase(ifs)
    mollist = []

    for mol in moldb.GetOEMols():
        molecule = oechem.OEMol(mol)
        mollist.append(normalize_molecule(molecule, title))

    if len(mollist) <= 1:
        mollist = mollist[0]

    ifs.close()

    return mollist


def normalize_molecule(molecule, title=None):
    """Normalize a copy of the molecule by checking aromaticity, adding explicit hydrogens and renaming by IUPAC name
    or given title

    Parameters
    ----------
    molecule: OEMol
        The molecule to be normalized:
    title: str
        Name of molecule. If none, will use IUPAC name

    Returns
    -------
    molcopy: OEMol
        A (copied) version of the normalized molecule

    """
    molcopy = oechem.OEMol(molecule)

    # Assign aromaticity.
    oechem.OEAssignAromaticFlags(molcopy, oechem.OEAroModelOpenEye)

    # Add hydrogens.
    oechem.OEAddExplicitHydrogens(molcopy)

    # Set title to IUPAC name.
    name = title
    if not name:
        name = oeiupac.OECreateIUPACName(molcopy)
    molcopy.SetTitle(name)

    # Check for any missing atom names, if found reassign all of them.
    if any([atom.GetName() == '' for atom in molcopy.GetAtoms()]):
        oechem.OETriposAtomNames(molcopy)

    return molcopy


def png_atoms_labeled(smiles, fname):
    """Write out png file of molecule with atoms labeled with their index.

    Parameters
    ----------
    smiles: str
        SMILES
    fname: str
        absolute path and filename for png

    """

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smiles)
    oedepict.OEPrepareDepiction(mol)

    width, height = 300, 200

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
    opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


def png_wiberg_labels(mol, fname, width=600, height=400):
    """
    Generate png figure of molecule. Bonds are labeled with Wiberg bond order

    Parameters
    ----------
    mol: OpenEye OEMol
    fname: str
        filename for png
    width: int
    height: int

    Returns
    -------
    bool:
    """

    oedepict.OEPrepareDepiction(mol)


    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    # opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
    # opts.SetAtomPropLabelFont(oedepict.OEFont(oechem.OEDarkGreen))

    bondlabel = LabelBondOrder()
    opts.SetBondPropertyFunctor(bondlabel)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    return oedepict.OERenderMolecule(fname, disp)


class LabelBondOrder(oedepict.OEDisplayBondPropBase):
    def __init__(self):
        oedepict.OEDisplayBondPropBase.__init__(self)

    def __call__(self, bond):
        bondOrder = bond.GetData('WibergBondOrder')
        label = "{:.2f}".format(bondOrder)
        return label

    def CreateCopy(self):
        copy = LabelBondOrder()
        return copy.__disown__()


def mol2_to_psi4json(infile):
    """

    Parameters
    ----------
    infile

    Returns
    -------

    """
    pass


def create_mapped_smiles(mol):
    """
    Generate an index-tagged explicit hydrogen SMILES.
    Exmaple:
    SMILES string for carbon monoxide "CO"
    With index-tagged explicit hydrogen SMILES this becomes
    '[H:3][C:1]([H:4])([H:5])[O:2][H:6]'

    Parameters
    ----------
    mol: OEMOl

    Returns
    -------
    index-tagged explicit hydrogen SMILES str

    """
    # Check if molecule already has explicit hydrogens
    HAS_HYDROGENS = oechem.OEHasExplicitHydrogens(mol)
    if not HAS_HYDROGENS:
        # Add explicit hydrogens
        oechem.OEAddExplicitHydrogens(mol)
    for atom in mol.GetAtoms():
        atom.SetMapIdx(atom.GetIdx() + 1)

    return oechem.OEMolToSmiles(mol)


def mol_to_tagged_smiles(infile, outfile):
    """
    Generate .smi from input mol with index-tagged explicit hydrogen SMILES
    Parameters
    ----------
    infile: str
        input molecule file
    outfile: str
        output smi file. Must be smi or ism

    """
    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open {} for reading".format(infile))

    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open {} for writing".format(outfile))
    if ofs.GetFormat() not in [oechem.OEFormat_ISM, oechem.OEFormat_SMI]:
        oechem.OEThrow.Fatal("Output format must be SMILES")

    for mol in ifs.GetOEMols():
        smiles = create_mapped_smiles(mol)
        oechem.OEWriteMolecule(ofs, mol)


def get_atom_map(tagged_smiles, molecule=None):
    """
    Maps tag index in SMILES to atom index in OEMol.
    Parameters
    ----------
    tagged_smiles: str
        index-tagged explicit Hydrogen SMILES string
    molecule: OEMOl
        The OEMol to get the map for. Default is None.
        When None, a new OEMol will be generated from SMILES. If an OEMol is provided, the map will be to that OEMol.

    Returns
    -------
    map: dict
        maps tag to index in returned OEMol. {map_idx: atom_idx}
    molecule: OEMol. Only if no OEMol was provided.

    """
    target = molecule
    if not molecule:
        target = oechem.OEMol()
        oechem.OESmilesToMol(target, tagged_smiles)

    # create maximum common substructure object
    mcss = oechem.OEMCSSearch(tagged_smiles, oechem.OEMCSType_Approximate)
    mcss.SetMaxMatches(1)
    mcss.SetMinAtoms(target.GetMaxAtomIdx())
    # Scoring function
    mcss.SetMCSFunc(oechem.OEMCSMaxAtoms())

    atom_map = {}
    unique = True
    t1 = time.time()
    for count, match in enumerate(mcss.Match(target, unique)):
        for ma in match.GetAtoms():
            atom_map[ma.pattern.GetMapIdx()] = ma.target.GetIdx()
    t2 = time.time()
    logger().info('MCSS took {} seconds'.format(t2-t1))

    m = oechem.OEGraphMol()
    oechem.OESubsetMol(m, match, True)
    logger().info("Match SMILES: {}".format(oechem.OEMolToSmiles(m)))

    if molecule:
        return atom_map
    else:
        return atom_map, target


def to_mapped_xyz(molecule, atom_map, conformer=None, xyz_format=False):
    """
    Generate xyz coordinates for molecule in the order given by the atom_map. atom_map is a dictionary that maps the
    tag on the SMILES to the atom idex in OEMol.
    Parameters
    ----------
    molecule
    atom_map
    conformer

    Returns
    -------
    str: XYZ file format

    """
    xyz = ""
    for k, mol in enumerate(molecule.GetConfs()):
        if k == conformer or conformer is None:
            if xyz_format:
                xyz += "{}\n".format(mol.GetMaxAtomIdx())
                xyz += "{}\n".format(mol.GetTitle())
            coords = oechem.OEFloatArray(mol.GetMaxAtomIdx() * 3)
            mol.GetCoords(coords)
            for mapping in range(1, len(atom_map)+1):
                idx = atom_map[mapping]
                atom = mol.GetAtom(oechem.OEHasAtomIdx(idx))
                syb = oechem.OEGetAtomicSymbol(atom.GetAtomicNum())

                xyz += "  {}      {:05.3f}   {:05.3f}   {:05.3f}\n".format(syb,
                                                                           coords[idx * 3],
                                                                           coords[idx * 3 + 1],
                                                                           coords[idx * 3 + 2])
    return xyz


# def run_psi4_json(tagged_smiles, molecule, driver, method, basis, properties=None, return_output=True, xyz_file=True):
#     """
#
#     Parameters
#     ----------
#     tagged_smiles
#     xyz
#     driver
#     method
#     basis
#     properties
#     return_output
#     xyz_file
#
#     Returns
#     -------
#
#     """
#     json_data = {}
#     json_data["tagged_smiles"] = tagged_smiles
#     json_data["molecule"] = molecule
#     json_data["driver"] = "property"
#     json_data["kwargs"] = {"properties": properties}
#     json_data["method"] = method
#     json_data["options"] = {"BASIS": basis}
#     json_data["return_output"] = return_output
#
#     name = molecule.split("\n")[2]
#     if xyz_file:
#         file = open("{}.xyz".format(name), 'w')
#         file.write(molecule)
#         file.close()
#
#     j = json.dumb(json_data, indent=4, sort_keys=True)
#     f = open("{}.input.json".format(name), 'w')
#     f.close()
#
#     psi4.json_wrapper.run_json(json_data)
#     j = json.dump(json_data, indent=4, sort_keys=True)
    f = open("{}.output.json".format(name))
