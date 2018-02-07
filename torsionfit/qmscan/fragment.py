#!/usr/bin/env python3
# (C) 2017 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

#############################################################################
# Depicts connected fragment combinations
#############################################################################


from itertools import combinations
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye import oemedchem

from openmoltools import openeye

import yaml
import os
from pkg_resources import resource_filename
import copy
import itertools

from torsionfit.qmscan import utils


def generate_fragments(inputf, output_dir, pdf=False, combinatorial=True, MAX_ROTORS=2):
    """

    Parameters
    ----------
    inputf:
    output_dir:
    pdf:
    combinatorial:
    MAX_ROTORS:

    Returns
    -------

    """
    ifs = oechem.oemolistream()

    fragments = {}
    smiles_unique = set()

    if ifs.open(inputf):
        for mol in ifs.GetOEGraphMols():
            print('fragmenting {}...'.format(mol.GetTitle()))
            charged, frags = _generate_fragments(mol)
            if combinatorial:
                smiles = smiles_with_combined(frags, charged, MAX_ROTORS=MAX_ROTORS)
            else:
                smiles = frag_to_smiles(frags, charged)

            smiles_unique.update(list(smiles.keys()))
            if pdf:
                oname = '{}.pdf'.format(mol.GetTitle())
                ToPdf(charged, oname, frags)
            fragments[mol.GetTitle()] = smiles

    # Generate oedatabase for all fragments
    base, ext = inputf.split('.')
    ofname = base + '_frags'
    utils.to_smi(list(smiles_unique), output_dir, ofname)


def _generate_fragments(mol):
    """
    This function generates fragments from a molecule.

    Parameters
    ----------
    mol: OEMol

    Returns
    -------
    charged: charged OEMOl
    frags: dict of AtomBondSet mapped to rotatable bond index the fragment was built up from.
    """

    charged = openeye.get_charges(mol, keep_confs=1)

    tagged_rings, tagged_fgroups = tag_molecule(charged)

    # Iterate over bonds
    frags = {}
    for bond in charged.GetBonds():
        if bond.IsRotor():
            atoms, bonds = _build_frag(bond=bond, mol=charged, tagged_fgroups=tagged_fgroups, tagged_rings=tagged_rings)
            atom_bond_set = _to_AtomBondSet(charged, atoms, bonds)
            frags[bond.GetIdx()] = atom_bond_set

    return charged, frags


def _tag_fgroups(mol, fgroups_smarts=None):
    """
    This function tags atoms and bonds of functional groups defined in fgroup_smarts. fgroup_smarts is a dictionary
    that maps functional groups to their smarts pattern. It can be user generated or from yaml file.

    Parameters
    ----------
    mol: Openeye OEMolGraph
    frgroups_smarts: dictionary of functional groups mapped to their smarts pattern.
        Default is None. It uses 'fgroup_smarts.yaml'

    Returns
    -------
    fgroup_tagged: dict
        a dictionary that maps indexed functional groups to corresponding atom and bond indices in mol

    """
    if not fgroups_smarts:
        # Load yaml file
        fn = resource_filename('torsionfit', os.path.join('qmscan', 'fgroup_smarts.yml'))
        fgroups_smarts = yaml.safe_load(open(fn, 'r'))
    fgroup_tagged = {}
    for f_group in fgroups_smarts:
        qmol = oechem.OEQMol()
        if not oechem.OEParseSmarts(qmol, fgroups_smarts[f_group]):
            print('OEParseSmarts failed')
        ss = oechem.OESubSearch(qmol)
        oechem.OEPrepareSearch(mol, ss)

        for i, match in enumerate(ss.Match(mol, True)):
            fgroup_atoms = set()
            for ma in match.GetAtoms():
                fgroup_atoms.add(ma.target.GetIdx())
                tag = oechem.OEGetTag('fgroup')
                ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))
            fgroup_bonds = set()
            for ma in match.GetBonds():
                #if not ma.target.IsInRing():
                fgroup_bonds.add(ma.target.GetIdx())
                tag =oechem.OEGetTag('fgroup')
                ma.target.SetData(tag, '{}_{}'.format(f_group, str(i)))

            fgroup_tagged['{}_{}'.format(f_group, str(i))] = (fgroup_atoms, fgroup_bonds)
    return fgroup_tagged


def _tag_rings(mol):
    """
    This function tags ring atom and bonds with ringsystem index

    Parameters
    ----------
    mol: OpenEye OEMolGraph

    Returns
    -------
    tagged_rings: dict
        maps ringsystem index to ring atom and bond indices

    """
    tagged_rings = {}
    nringsystems, parts = oechem.OEDetermineRingSystems(mol)
    for ringidx in range(1, nringsystems +1):
        ringidx_atoms = set()
        for atom in mol.GetAtoms():
            if parts[atom.GetIdx()] == ringidx:
                ringidx_atoms.add(atom.GetIdx())
                tag = oechem.OEGetTag('ringsystem')
                atom.SetData(tag, ringidx)
        # Find bonds in ring and tag
        ringidx_bonds = set()
        for a_idx in ringidx_atoms:
            atom = mol.GetAtom(oechem.OEHasAtomIdx(a_idx))
            for bond in atom.GetBonds():
                nbrAtom = bond.GetNbr(atom)
                nbrIdx = nbrAtom.GetIdx()
                if nbrIdx in ringidx_atoms and nbrIdx != a_idx:
                    ringidx_bonds.add(bond.GetIdx())
                    tag = oechem.OEGetTag('ringsystem')
                    bond.SetData(tag, ringidx)
        tagged_rings[ringidx] = (ringidx_atoms, ringidx_bonds)
    return tagged_rings


def _ring_fgroup_union(mol, tagged_rings, tagged_fgroups, wbo_threshold=1.2):
    """
    This function combines rings and fgroups that are conjugated (the bond between them has a Wiberg bond order > 1.2)

    Parameters
    ----------
    mol: OpenEye OEMolGraph
    tagged_rings: dict
        map of ringsystem indices to ring atom and bond indices
    tagged_fgroup: dict
        map of fgroup to fgroup atom and bond indices

    Returns
    -------
    tagged_fgroup: dict
        updated tagged_fgroup mapping with rings that shouldn't be fragmented from fgroups
    """
    ring_idxs = list(tagged_rings.keys())
    fgroups = list(tagged_fgroups.keys())
    tagged_fgroups = copy.deepcopy(tagged_fgroups)

    # Check if fgroups are overlapping. If they are - combine them.
    for func_group_1, func_group_2 in itertools.combinations(list(tagged_fgroups.keys()), 2):
        atoms_intersection = tagged_fgroups[func_group_1][0].intersection(tagged_fgroups[func_group_2][0])
        if len(atoms_intersection) > 1:
            # Combine fgroups
            atoms_union = tagged_fgroups[func_group_1][0].union(tagged_fgroups[func_group_2][0])
            bonds_union = tagged_fgroups[func_group_1][-1].union(tagged_fgroups[func_group_2][-1])
            tagged_fgroups[func_group_1] = (atoms_union, bonds_union)
            tagged_fgroups[func_group_2] = (atoms_union, bonds_union)
    for idx in ring_idxs:
        for fgroup in fgroups:
            atom_intersection = tagged_rings[idx][0].intersection(tagged_fgroups[fgroup][0])
            if len(atom_intersection) > 1:
                # Must include ring if including fgroup. Add ring atoms and bonds to fgroup
                atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                tagged_fgroups[fgroup] = (atoms_union, bonds_union)
            elif len(atom_intersection) > 0:
                # Check Wiberg bond order of bond
                # First find bond connectiong fgroup and ring
                atom = mol.GetAtom(oechem.OEHasAtomIdx(atom_intersection.pop()))
                for a in atom.GetAtoms():
                    if a.GetIdx() in tagged_fgroups[fgroup][0]:
                        bond = mol.GetBond(a, atom)
                        if bond.GetData('WibergBondOrder') > wbo_threshold:
                            # Don't cut off ring.
                            atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                            tagged_fgroups[fgroup] = (atoms_union, bonds_union)
                        # Should I also combine non-rotatable rings? This will pick up the alkyn in ponatinib?
                        if not bond.IsRotor():
                            atoms_union = tagged_rings[idx][0].union(tagged_fgroups[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(tagged_fgroups[fgroup][-1])
                            tagged_fgroups[fgroup] = (atoms_union, bonds_union)
                        # Do something when it's neither for edge cases (Afatinib)?
    return tagged_fgroups


def tag_molecule(mol, func_group_smarts=None):
    """
    Tags atoms and molecules in functional groups and ring systems. The molecule gets tagged and the function returns
    a 2 dictionaries that map
    1) ring system indices to corresponding atoms and bonds indices in molecule
    2) functional groups to corresponding atoms and bonds indices in molecule

    Parameters
    ----------
    mol: OEMol
    func_group_smarts: dict
        dictionary mapping functional groups to SMARTS. Default is None and uses shipped yaml file.

    Returns
    -------
    tagged_rings: dict
        mapping of ring system index to corresponding atoms and bonds in molecule. Each index maps to 2 sets. First set
        includes atom indices and second set includes bonds indices
    tagged_func_group: dict
        mapping of functional group to corresponding atoms and bonds in molecule. Each functional group maps to 2 sets.
        The first set is atom indices, the second set is bond indices.

    """
    tagged_func_group = _tag_fgroups(mol, func_group_smarts)
    tagged_rings = _tag_rings(mol)

    tagged_func_group = _ring_fgroup_union(mol=mol, tagged_fgroups=tagged_func_group, tagged_rings=tagged_rings)

    return tagged_rings, tagged_func_group


def _is_fgroup(fgroup_tagged, element):
    """
    This function checks if an atom or a bond is part of a tagged fgroup.

    Parameters
    ----------
    atom: Openeye Atom Base
    fgroup_tagged: dict of indexed functional group and corresponding atom and bond indices

    Returns
    -------
    atoms, bonds: sets of atom and bond indices if the atom is tagged, False otherwise

    """
    try:
        fgroup = element.GetData('fgroup')
        atoms, bonds = fgroup_tagged[fgroup]
        return atoms, bonds
    except ValueError:
        return False


def _to_AtomBondSet(mol, atoms, bonds):
    """
    Builds OpeneyeAtomBondet from atoms and bonds set of indices
    Parameters
    ----------
    mol: Openeye OEMolGraph
    atoms: Set of atom indices
    bonds: Set of bond indices

    Returns
    -------
    AtomBondSet: Openeye AtomBondSet of fragment
    """

    AtomBondSet = oechem.OEAtomBondSet()
    for a_idx in atoms:
        AtomBondSet.AddAtom(mol.GetAtom(oechem.OEHasAtomIdx(a_idx)))
    for b_idx in bonds:
        AtomBondSet.AddBond(mol.GetBond(oechem.OEHasBondIdx(b_idx)))
    return AtomBondSet


def _is_ortho(bond, rot_bond, next_bond):
    """
    This function checks if a bond is ortho to the rotatable bond
    Parameters
    ----------
    bond: OEBondBase
        current bond to check if it's ortho
    rot_bond: OEBondBase
        the rotatable bond the bond needs to be ortho to
    next_bond: OEBondBase
        The bond between rot_bond and bond if bond is ortho to rot_bond

    Returns
    -------
    bool: True if ortho, False if not
    """

    bond_attached = set()
    rot_attached = set()

    for bon in [bond.GetBgn(), bond.GetEnd()]:
        for b in bon.GetBonds():
            bond_attached.add(b.GetIdx())
    for bon in [rot_bond.GetBgn(), rot_bond.GetEnd()]:
        for b in bon.GetBonds():
            rot_attached.add(b.GetIdx())

    if not next_bond.IsInRing():
        next_attached = set()
        for bon in [next_bond.GetBgn(), next_bond.GetEnd()]:
            for b in bon.GetBonds():
                next_attached.add(b.GetIdx())

    intersection = (bond_attached & rot_attached)
    if not bool(intersection) and not next_bond.IsInRing():
        # Check if it's ortho to next bond
        intersection = (bond_attached & next_attached)
    return bool(intersection)


def _build_frag(bond, mol, tagged_fgroups, tagged_rings):
    """
    This functions builds a fragment around a rotatable bond. It grows out one bond in all directions
    If the next atoms is in a ring or functional group, it keeps that.
    If the next bond has a Wiberg bond order > 1.2, grow another bond and check next bond's Wiberg bond order.

    Parameters
    ----------
    bond: OpenEye bond
    mol: OpenEye OEMolGraph
    tagged_fgroups: dict
        maps functional groups to atoms and bond indices on mol
    tagged_rings: dict
        maps ringsystem index to atom and bond indices in mol

    Returns
    -------
    atoms, bonds: sets of atom and bond indices for fragment
    """


    atoms = set()
    bonds = set()
    b_idx = bond.GetIdx()
    bonds.add(b_idx)
    beg = bond.GetBgn()
    end = bond.GetEnd()
    beg_idx = beg.GetIdx()
    end_idx = end.GetIdx()

    atoms.add(beg_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol=mol, rotor_bond=bond, atom=beg, pair=end, fgroup_tagged=tagged_fgroups,
                                          tagged_rings=tagged_rings)
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    atoms.add(end_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol=mol, rotor_bond=bond, atom=end, pair=beg, fgroup_tagged=tagged_fgroups,
                                          tagged_rings=tagged_rings)
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    return atoms, bonds


def iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, i=0):
    """
    This function iterates over neighboring atoms and checks if it's part of a functional group, ring, or if the next
    bond has a Wiberg bond order > 1.2.

    Parameters
    ----------
    mol: Openeye OEMolGraph
    atom: Openeye AtomBase
        atom that will iterate over
    paired: OpeneEye AtomBase
        atom that's bonded to this atom in rotor_bond
    fgroup_tagged: dict
        map of functional group and atom and bond indices in mol
    tagged_rings: dict
        map of ringsystem index and atom and bond indices in mol
    rotor_bond: Openeye Bond base
        rotatable bond that the fragment is being built on

    Returns
    -------
    atoms, bonds: sets of atom and bond indices of the fragment

    """
    def _iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, atoms_2, bonds_2, i=0):


        for a in atom.GetAtoms():
            if a.GetIdx() == pair.GetIdx():
                continue
            a_idx = a.GetIdx()
            next_bond = mol.GetBond(a, atom)
            nb_idx = next_bond.GetIdx()
            atoms_2.add(a_idx)
            if nb_idx in bonds_2:
                FGROUP_RING = False
                try:
                    ring_idx = a.GetData('ringsystem')
                    fgroup = a.GetData('fgroup')
                    FGROUP_RING = True
                except ValueError:
                    try:
                        ring_idx = atom.GetData('ringsystem')
                        fgroup = atom.GetData('fgroup')
                        FGROUP_RING = True
                    except ValueError:
                        continue
                if FGROUP_RING:
                    # Add ring and continue
                    ratoms, rbonds = tagged_rings[ring_idx]
                    atoms_2 = atoms_2.union(ratoms)
                    bonds_2 = bonds_2.union(rbonds)
                    rs_atoms, rs_bonds = _ring_substiuents(mol=mol, bond=next_bond, rotor_bond=rotor_bond,
                                                          tagged_rings=tagged_rings, ring_idx=ring_idx,
                                                          fgroup_tagged=fgroup_tagged)
                    atoms_2 = atoms_2.union(rs_atoms)
                    bonds_2 = bonds_2.union(rs_bonds)
                    continue

            if i > 0:
                wiberg = next_bond.GetData('WibergBondOrder')
                if wiberg < 1.2:
                    continue

            bonds_2.add(nb_idx)
            if a.IsInRing():
                ring_idx = a.GetData('ringsystem')
                ratoms, rbonds = tagged_rings[ring_idx]
                atoms_2 = atoms_2.union(ratoms)
                bonds_2 = bonds_2.union(rbonds)
                # Find non-rotatable sustituents
                rs_atoms, rs_bonds = _ring_substiuents(mol=mol, bond=next_bond, rotor_bond=rotor_bond,
                                                      tagged_rings=tagged_rings, ring_idx=ring_idx,
                                                      fgroup_tagged=fgroup_tagged)
                atoms_2 = atoms_2.union(rs_atoms)
                bonds_2 = bonds_2.union(rs_bonds)
            fgroup = _is_fgroup(fgroup_tagged, element=a)
            if fgroup: # and i < 1:
                atoms_2 = atoms_2.union(fgroup[0])
                bonds_2 = bonds_2.union(fgroup[-1])
                # if something is in a ring - have a flag? Then use that to continue iterating and change the flag

            for nb_a in a.GetAtoms():
                nn_bond = mol.GetBond(a, nb_a)
                if (nn_bond.GetData('WibergBondOrder') > 1.2) and (not nn_bond.IsInRing()) and (not nn_bond.GetIdx() in bonds_2):
                    # Check the degree of the atoms in the bond
                    deg_1 = a.GetDegree()
                    deg_2 = nb_a.GetDegree()
                    if deg_1 == 1 or deg_2 == 1:
                        continue

                    atoms_2.add(nb_a.GetIdx())
                    bonds_2.add(nn_bond.GetIdx())
                    i += 1
                    _iterate_nbratoms(mol, nn_bond, nb_a, pair, fgroup_tagged, tagged_rings, atoms_2, bonds_2, i=i)
        return atoms_2, bonds_2
    return _iterate_nbratoms(mol, rotor_bond, atom, pair, fgroup_tagged, tagged_rings, atoms_2=set(), bonds_2=set(), i=0)


def _ring_substiuents(mol, bond, rotor_bond, tagged_rings, ring_idx, fgroup_tagged):
    """
    This function finds ring substituents that shouldn't be cut off

    Parameters
    ----------
    mol: OpeneEye OEMolGraph
    bond: OpenEye Bond Base
        current bond that the iterator is looking at
    rotor_bond: OpeneEye Bond Base
        rotatable bond that fragment is being grown on
    tagged_rings: dict
        mapping of ring index and atom and bonds indices
    ring_idx: int
        ring index
    fgroup_tagged: dict
        mapping of functional group and atom and bond indices

    Returns
    -------
    rs_atoms, rs_bonds: sets of ring substituents atoms and bonds indices

    """
    rs_atoms = set()
    rs_bonds = set()
    r_atoms, r_bonds = tagged_rings[ring_idx]
    for a_idx in r_atoms:
        atom = mol.GetAtom(oechem.OEHasAtomIdx(a_idx))
        for a in atom.GetAtoms():
            if a.GetIdx() in rs_atoms:
                continue
            fgroup = False
            rs_bond = mol.GetBond(atom, a)
            if not a.IsInRing():
                if not rs_bond.IsRotor():
                    rs_atoms.add(a.GetIdx())
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = _is_fgroup(fgroup_tagged, element=a)
                elif _is_ortho(rs_bond, rotor_bond, bond):
                    # Keep bond and attached atom.
                    rs_atoms.add(a.GetIdx())
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = _is_fgroup(fgroup_tagged, element=a)
                if fgroup:
                    rs_atoms = rs_atoms.union(fgroup[0])
                    rs_bonds = rs_bonds.union(fgroup[-1])
            else:
                # Check if they are in the same ring
                r_idx2 = a.GetData('ringsystem')
                if r_idx2 != ring_idx:
                    if _is_ortho(rs_bond, rotor_bond, bond):
                        # Add ring system
                        rs_bonds.add(rs_bond.GetIdx())
                        r2_atoms, r2_bonds = tagged_rings[r_idx2]
                        rs_atoms = rs_atoms.union(r2_atoms)
                        rs_bonds = rs_bonds.union(r2_bonds)

    return rs_atoms, rs_bonds


def frag_to_smiles(frags, mol):
    """
    Convert fragments (AtomBondSet) to smiles string
    Parameters
    ----------
    frags
    mol

    Returns
    -------
    smiles: list of smiles strings

    """
    smiles = {}
    for frag in frags:
        fragatompred = oechem.OEIsAtomMember(frag.GetAtoms())
        fragbondpred = oechem.OEIsBondMember(frag.GetBonds())

        fragment = oechem.OEGraphMol()
        adjustHCount = True
        oechem.OESubsetMol(fragment, mol, fragatompred, fragbondpred, adjustHCount)
        s = oechem.OEMolToSmiles(fragment)
        if s not in smiles:
            smiles[s] = []
        smiles[s].append(frag)
    return smiles


def smiles_with_combined(frags, mol, MAX_ROTORS=2):
    """
    Generates Smiles:frags mapping for fragments and fragment combinations with less than MAX_ROTORS rotatable bonds

    Parameters
    ----------
    frags: dict of rot band mapped to AtomBondSet
    mol: OpenEye Mol

    Returns
    -------
    smiles: dict of smiles sting to fragment

    """
    frag_list = list(frags.values())
    comb_list = GetFragmentAtomBondSetCombinations(frag_list, MAX_ROTORS=MAX_ROTORS)

    combined_list = comb_list + frag_list

    smiles = frag_to_smiles(combined_list, mol)

    return smiles


# def main(argv=[__name__]):
#
#     itf = oechem.OEInterface(InterfaceData)
#     oedepict.OEConfigure2DMolDisplayOptions(itf)
#     oedepict.OEConfigureReportOptions(itf)
#
#     if not oechem.OEParseCommandLine(itf, argv):
#         return 1
#
#     iname = itf.GetString("-in")
#     oname = itf.GetString("-out")
#
#     max_rotors = itf.GetInt("-max_rotors")
#     min_rotors = itf.GetInt("-min_rotors")
#
#     frags = itf.GetList("-frags")
#
#     pagebypage = itf.GetBool("-pagebypage")
#
#     # check input/output files
#
#     ifs = oechem.oemolistream()
#     if not ifs.open(iname):
#         oechem.OEThrow.Fatal("Cannot open input file!")
#
#     ext = oechem.OEGetFileExtension(oname)
#     if not pagebypage and not oedepict.OEIsRegisteredMultiPageImageFile(ext):
#         oechem.OEThrow.Warning("Report will be generated into separate pages!")
#         pagebypage = True
#
#     # read a molecule
#
#     mol = oechem.OEGraphMol()
#     if not oechem.OEReadMolecule(ifs, mol):
#         oechem.OEThrow.Fatal("Cannot read input file!")
#     oedepict.OEPrepareDepiction(mol)
#
#     # initialize fragmentation function
#
#     fragfunc = GetFragmentationFunction(itf)
#
#     # initialize multi-page report
#
#     ropts = oedepict.OEReportOptions()
#     oedepict.OESetupReportOptions(ropts, itf)
#     ropts.SetFooterHeight(25.0)
#     ropts.SetHeaderHeight(ropts.GetPageHeight() / 4.0)
#     report = oedepict.OEReport(ropts)
#
#     # setup depiction options
#
#     opts = oedepict.OE2DMolDisplayOptions()
#     oedepict.OESetup2DMolDisplayOptions(opts, itf)
#     cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
#     opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
#     opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
#     opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
#     opts.SetAtomLabelFontScale(1.2)
#
#     # depict molecule with fragment combinations
#
#     DepictMoleculeWithFragmentCombinations(report, mol, frags, opts, max_rotors, min_rotors)
#
#     if pagebypage:
#         oedepict.OEWriteReportPageByPage(oname, report)
#     else:
#         oedepict.OEWriteReport(oname, report)
#
#     return 0
#

def ToPdf(mol, oname, frags):#, fragcombs):
    """
    Parameters
    ----------
    mol: charged OEMolGraph
    oname: str
        Output file name
    Returns
    -------

    """
    itf = oechem.OEInterface()
    oedepict.OEPrepareDepiction(mol)

    ropts = oedepict.OEReportOptions()
    oedepict.OESetupReportOptions(ropts, itf)
    ropts.SetFooterHeight(25.0)
    ropts.SetHeaderHeight(ropts.GetPageHeight() / 4.0)
    report = oedepict.OEReport(ropts)

    # setup decpiction options
    opts = oedepict.OE2DMolDisplayOptions()
    oedepict.OESetup2DMolDisplayOptions(opts, itf)
    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    opts.SetTitleLocation(oedepict.OETitleLocation_Hidden)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)
    opts.SetAtomLabelFontScale(1.2)

    DepictMoleculeWithFragmentCombinations(report, mol, frags, opts)

    oedepict.OEWriteReport(oname, report)

    return 0


def OeMolToGraph(oemol):
    """
    Convert charged molecule to networkX graph and add WiberBondOrder as edge weight

    Parameters
    ----------
    mol: charged OEMolGraph

    Returns
    -------
    G: NetworkX Graph of molecule

    """
    import networkx as nx
    G = nx.Graph()
    for atom in oemol.GetAtoms():
        G.add_node(atom.GetIdx(), name=atom.GetName())
    for bond in oemol.GetBonds():
        try:
            fgroup = bond.GetData('fgroup')
        except:
            fgroup = False
        G.add_edge(bond.GetBgnIdx(), bond.GetEndIdx(), weight=bond.GetData("WibergBondOrder"), index=bond.GetIdx(),
                   aromatic=bond.IsAromatic(), in_ring=bond.IsInRing(), fgroup=fgroup)
    return G


def FragGraph(G, bondOrderThreshold=1.2):
    """
    Fragment all bonds with Wiberg Bond Order less than threshold

    Parameters
    ----------
    G: NetworkX graph
    bondOrderThreshold: int
        thershold for fragmenting graph. Default 1.2

    Returns
    -------
    subgraphs: list of subgraphs
    """
    import networkx as nx

    ebunch = []
    for node in G.edge:
        if G.degree(node) <= 1:
            continue
        for node2 in G.edge[node]:
            if G.edge[node][node2]['weight'] < bondOrderThreshold and G.degree(node2) >1 \
                    and not G.edge[node][node2]['aromatic'] and not G.edge[node][node2]['in_ring']\
                    and not G.edge[node][node2]['fgroup']:
                ebunch.append((node, node2))
    # Cut molecule
    G.remove_edges_from(ebunch)
    # Generate fragments
    subgraphs = list(nx.connected_component_subgraphs(G))
    return subgraphs


def subgraphToAtomBondSet(graph, subgraph, oemol):
    """
    Build Openeye AtomBondSet from subrgaphs for enumerating fragments recipe

    Parameters
    ----------
    graph: NetworkX graph
    subgraph: NetworkX subgraph
    oemol: Openeye OEMolGraph

    Returns
    ------
    atomBondSet: Openeye oechem atomBondSet
    """
    # Build openeye atombondset from subgraphs
    atomBondSet = oechem.OEAtomBondSet()
    for node in subgraph.node:
        atomBondSet.AddAtom(oemol.GetAtom(oechem.OEHasAtomIdx(node)))
    for node1, node2 in subgraph.edges():
        index = graph.edge[node1][node2]['index']
        atomBondSet.AddBond(oemol.GetBond(oechem.OEHasBondIdx(index)))
    return atomBondSet


def SmilesToFragments(smiles, fgroup_smarts, bondOrderThreshold=1.2, chargesMol=True):
    """
    Fragment molecule at bonds below Bond Order Threshold

    Parameters
    ----------
    smiles: str
        smiles string of molecule to fragment

    Returns
    -------
    frags: list of OE AtomBondSets

    """
    # Charge molecule
    mol = oechem.OEGraphMol()
    oemol = openeye.smiles_to_oemol(smiles)
    charged = openeye.get_charges(oemol, keep_confs=1)

    # Tag functional groups
    tag_fgroups(charged, fgroups_smarts=fgroup_smarts)

    # Generate fragments
    G = OeMolToGraph(charged)
    subraphs = FragGraph(G, bondOrderThreshold=bondOrderThreshold)

    frags = []
    for subraph in subraphs:
        frags.append(subgraphToAtomBondSet(G, subraph, charged))

    if chargesMol:
        return frags, charged
    else:
        return frags


def DepictMoleculeWithFragmentCombinations(report, mol, frags, opts): #fragcombs, opts):

    stag = "fragment idx"
    itag = oechem.OEGetTag(stag)
    for fidx, frag in enumerate(frags):
        for bond in frags[frag].GetBonds():
            bond.SetData(itag, fidx)

    # setup depiction styles

    nrfrags = len(frags)
    colors = [c for c in oechem.OEGetLightColors()]
    if len(colors) < nrfrags:
        colors = [c for c in oechem.OEGetColors(oechem.OEYellowTint, oechem.OEDarkOrange, nrfrags)]

    bondglyph = ColorBondByFragmentIndex(colors, itag)

    lineWidthScale = 0.75
    fadehighlight = oedepict.OEHighlightByColor(oechem.OEGrey, lineWidthScale)

    # depict each fragment combinations

    #for frag in fragcombs:
    for frag in frags:

        cell = report.NewCell()
        disp = oedepict.OE2DMolDisplay(mol, opts)

        fragatoms = oechem.OEIsAtomMember(frags[frag].GetAtoms())
        fragbonds = oechem.OEIsBondMember(frags[frag].GetBonds())

        notfragatoms = oechem.OENotAtom(fragatoms)
        notfragbonds = oechem.OENotBond(fragbonds)

        oedepict.OEAddHighlighting(disp, fadehighlight, notfragatoms, notfragbonds)

        bond = mol.GetBond(oechem.OEHasBondIdx(frag))

        atomBondSet = oechem.OEAtomBondSet()
        atomBondSet.AddBond(bond)
        atomBondSet.AddAtom(bond.GetBgn())
        atomBondSet.AddAtom(bond.GetEnd())

        hstyle = oedepict.OEHighlightStyle_BallAndStick
        hcolor = oechem.OEColor(oechem.OELightBlue)
        oedepict.OEAddHighlighting(disp, hcolor, hstyle, atomBondSet)

        #oegrapheme.OEAddGlyph(disp, bondglyph, fragbonds)

        oedepict.OERenderMolecule(cell, disp)

    # depict original fragmentation in each header

    cellwidth, cellheight = report.GetHeaderWidth(), report.GetHeaderHeight()
    opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

    bondlabel = LabelBondOrder()
    opts.SetBondPropertyFunctor(bondlabel)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    #oegrapheme.OEAddGlyph(disp, bondglyph, oechem.IsTrueBond())

    headerpen = oedepict.OEPen(oechem.OEWhite, oechem.OELightGrey, oedepict.OEFill_Off, 2.0)
    for header in report.GetHeaders():
        oedepict.OERenderMolecule(header, disp)
        oedepict.OEDrawBorder(header, headerpen)


class ColorBondByFragmentIndex(oegrapheme.OEBondGlyphBase):
    def __init__(self, colorlist, tag):
        oegrapheme.OEBondGlyphBase.__init__(self)
        self.colorlist = colorlist
        self.tag = tag

    def RenderGlyph(self, disp, bond):

        bdisp = disp.GetBondDisplay(bond)
        if bdisp is None or not bdisp.IsVisible():
            return False

        if not bond.HasData(self.tag):
            return False

        linewidth = disp.GetScale() / 2.0
        color = self.colorlist[bond.GetData(self.tag)]
        pen = oedepict.OEPen(color, color, oedepict.OEFill_Off, linewidth)

        adispB = disp.GetAtomDisplay(bond.GetBgn())
        adispE = disp.GetAtomDisplay(bond.GetEnd())

        layer = disp.GetLayer(oedepict.OELayerPosition_Below)
        layer.DrawLine(adispB.GetCoords(), adispE.GetCoords(), pen)

        return True

    def ColorBondByFragmentIndex(self):
        return ColorBondByFragmentIndex(self.colorlist, self.tag).__disown__()

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

def IsAdjacentAtomBondSets(fragA, fragB):

    for atomA in fragA.GetAtoms():
        for atomB in fragB.GetAtoms():
            if atomA.GetBond(atomB) is not None:
                return True
    return False

def CountRotorsInFragment(fragment):
    return sum([bond.IsRotor() for bond in fragment.GetBonds()])

def IsAdjacentAtomBondSetCombination(fraglist):

    parts = [0] * len(fraglist)
    nrparts = 0

    for idx, frag in enumerate(fraglist):
        if parts[idx] != 0:
            continue

        nrparts += 1
        parts[idx] = nrparts
        TraverseFragments(frag, fraglist, parts, nrparts)

    return (nrparts == 1)


def TraverseFragments(actfrag, fraglist, parts, nrparts):

    for idx, frag in enumerate(fraglist):
        if parts[idx] != 0:
            continue

        if not IsAdjacentAtomBondSets(actfrag, frag):
            continue

        parts[idx] = nrparts
        TraverseFragments(frag, fraglist, parts, nrparts)


def CombineAndConnectAtomBondSets(fraglist):

    # combine atom and bond sets

    combined = oechem.OEAtomBondSet()
    for frag in fraglist:
        for atom in frag.GetAtoms():
            combined.AddAtom(atom)
        for bond in frag.GetBonds():
            combined.AddBond(bond)

    # add connecting bonds

    for atomA in combined.GetAtoms():
        for atomB in combined.GetAtoms():
            if atomA.GetIdx() < atomB.GetIdx():
                continue

            bond = atomA.GetBond(atomB)
            if bond is None:
                continue
            if combined.HasBond(bond):
                continue

            combined.AddBond(bond)

    return combined


def GetFragmentAtomBondSetCombinations(fraglist, MAX_ROTORS=2, MIN_ROTORS=1):
    """
    Enumerate connected combinations from list of fragments
    Parameters
    ----------
    mol: OEMolGraph
    fraglist: list of OE AtomBondSet
    MAX_ROTORS: int
        min rotors in each fragment combination
    MIN_ROTORS: int
        max rotors in each fragment combination

    Returns
    -------
    fragcombs: list of connected combinations (OE AtomBondSet)
    """

    fragcombs = []

    nrfrags = len(fraglist)
    for n in range(1, nrfrags):

        for fragcomb in combinations(fraglist, n):

            if IsAdjacentAtomBondSetCombination(fragcomb):

                frag = CombineAndConnectAtomBondSets(fragcomb)

                if (CountRotorsInFragment(frag) <= MAX_ROTORS) and (CountRotorsInFragment(frag) >= MIN_ROTORS):
                    fragcombs.append(frag)

    return fragcombs


def GetFragmentationFunction(itf):

    fstring = itf.GetString("-fragtype")
    if fstring == "funcgroup":
        return oemedchem.OEGetFuncGroupFragments
    if fstring == "ring-chain":
        return oemedchem.OEGetRingChainFragments
    if fstring == 'bemis-murko':
        return oemedchem.OEGetBemisMurcko
    return oemedchem.OEGetRingLinkerSideChainFragments


# InterfaceData = '''
# !CATEGORY "input/output options"
#     !PARAMETER -in
#       !ALIAS -i
#       !TYPE string
#       !REQUIRED true
#       !KEYLESS 1
#       !VISIBILITY simple
#       !BRIEF Input filename
#     !END
#
#     !PARAMETER -out
#       !ALIAS -o
#       !TYPE string
#       !REQUIRED true
#       !KEYLESS 2
#       !VISIBILITY simple
#       !BRIEF Output filename
#     !END
# !END
#
# !CATEGORY "fragmentation options"
#     !PARAMETER -fragtype
#       !ALIAS -ftype
#       !TYPE string
#       !REQUIRED false
#       !KEYLESS 3
#       !DEFAULT funcgroup
#       !LEGAL_VALUE funcgroup
#       !LEGAL_VALUE ring-chain
#       !LEGAL_VALUE ring-linker-sidechain
#       !LEGAL_VALUE bemis-murcko
#       !VISIBILITY simple
#       !BRIEF Fragmentation type
#     !END
# !END
#
# !CATEGORY "report options"
#
#     !PARAMETER -pagebypage
#       !ALIAS -p
#       !TYPE bool
#       !REQUIRED false
#       !DEFAULT false
#       !VISIBILITY simple
#       !BRIEF Write individual numbered separate pages
#     !END
#
# !CATEGORY "max rotors option"
#     !PARAMETER -max_rotors
#         !ALIAS -m
#         !TYPE int
#         !REQUIRED false
#         !DEFAULT 3
#         !BRIEF max rotatable bonds
#     !END
#
# !CATEGORY "min rotors option"
#     !PARAMETER -min_rotors
#         !ALIAS -n
#         !TYPE int
#         !REQUIRED false
#         !DEFAULT 0
#         !BRIEF max rotatable bonds
#     !END
#
# !CATEGORY "fragments"
#     !PARAMETER -frags
#         !ALIAS -f
#         !TYPE list
#         !REQUIRED true
#         !BRIEF fragments to rebuild
#     !END
#
# !END
# '''
#
# if __name__ == "__main__":
#     sys.exit(main(sys.argv))
