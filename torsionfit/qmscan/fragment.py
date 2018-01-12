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



import sys
from itertools import combinations
from openeye import oechem
from openeye import oedepict
from openeye import oegrapheme
from openeye import oemedchem

import networkx as nx
from openmoltools import openeye


def tag_fgroups(mol, fgroups_smarts):
    """

    Parameters
    ----------
    mol: Openeye OEMolGraph
    frgroups_smarts: dictionary of functional groups mapped to their smarts pattern

    Returns
    -------
    fgroup_tagged: dict
        a dictionary that maps indexed functional groups to corresponding atom and bond indices in mol

    """

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


def tag_rings(mol):
    """

    Parameters
    ----------
    mol

    Returns
    -------

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


def ring_fgroup_union(mol, tagged_rings, fgroup_tagged):
    """

    Parameters
    ----------
    tagged_rings
    fgroup_tagged

    Returns
    -------

    """
    ring_idxs = list(tagged_rings.keys())
    fgroups = list(fgroup_tagged.keys())
    for idx in ring_idxs:
        for fgroup in fgroups:
            atom_intersection = tagged_rings[idx][0].intersection(fgroup_tagged[fgroup][0])
            if len(atom_intersection) > 1:
                # Must include ring if including fgroup. Add ring atoms and bonds to fgroup
                atoms_union = tagged_rings[idx][0].union(fgroup_tagged[fgroup][0])
                bonds_union = tagged_rings[idx][-1].union(fgroup_tagged[fgroup][-1])
                fgroup_tagged[fgroup] = (atoms_union, bonds_union)
            elif len(atom_intersection) > 0:
                # Check Wiberg bond order of bond
                # First find bond connectiong fgroup and ring
                atom = mol.GetAtom(oechem.OEHasAtomIdx(atom_intersection.pop()))
                for a in atom.GetAtoms():
                    if a.GetIdx() in fgroup_tagged[fgroup][0]:
                        bond = mol.GetBond(a, atom)
                        if bond.GetData('WibergBondOrder') > 1.2:
                            # Don't cut off ring.
                            atoms_union = tagged_rings[idx][0].union(fgroup_tagged[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(fgroup_tagged[fgroup][-1])
                            fgroup_tagged[fgroup] = (atoms_union, bonds_union)
                        # Should I also combine non-rotatable rings? This will pick up the alkyn in ponatinib?
                        if not bond.IsRotor():
                            atoms_union = tagged_rings[idx][0].union(fgroup_tagged[fgroup][0])
                            bonds_union = tagged_rings[idx][-1].union(fgroup_tagged[fgroup][-1])
                            fgroup_tagged[fgroup] = (atoms_union, bonds_union)
                        # Do something when it's neither for edge cases (Afatinib)?

    return fgroup_tagged


def is_fgroup(fgroup_tagged, atom=None, bond=None):
    """

    Parameters
    ----------
    atom: Openeye Atom Base
    fgroup_tagged: dict of indexed functional group and corresponding atom and bond indices

    Returns
    -------
    atoms, bonds: sets of atom and bond indices if the atom is tagged, False otherwise

    """
    if atom:
        try:
            fgroup = atom.GetData('fgroup')
            atoms, bonds = fgroup_tagged[fgroup]
            return atoms, bonds
        except ValueError:
            return False
    if bond:
        try:
            fgroup = bond.GetData('fgroup')
            atoms, bonds = fgroup_tagged[fgroup]
            return atoms, bonds
        except ValueError:
            return False


def is_hbond(bond):
    """
    Checks if the bond involves a hydrogen
    Parameters
    ----------
    bond: Openeye Bond Base

    Returns
    -------
    bond: Openeye Bond Base if bond involves hydrogen. False otherwise

    """
    hbond = None
    beg = bond.GetBgn()
    end = bond.GetEnd()
    if beg.IsHydrogen() or end.IsHydrogen():
        hbond = bond
    return hbond


def find_ring(a, mol, fgroup_tagged, ratoms_l = [], rbonds_l = [], rot_bond=None, next_bond=None, fg_ring=[],
              ratoms=set(), rbonds=set()):
    """
    This function finds the rest of the ring system that atom a is part of and returns the atom and bond indices of
    the ring atoms and bonds.

    Parameters
    ----------
    a: Openeye AtomBase
        The atom that's in a ring
    ratoms_l: list
        ring atom list of indices
    rbonds_l: list
        ring bonds list of indices
    rot_bond: Openeye BondBase
        the bond attached to a that is not in the ring. If not None, will check for ortho constituents instead of saving
        all non rotatable constituents
    next_bond: Openeye BondBase
        The bond next to the rotatable bond. This is used to check for ortho.

    Returns
    -------
    rbonds_l, ratoms_l: lists of ring atoms and bonds indices

    """
    if a.GetIdx() not in ratoms_l:
        ratoms_l.append(a.GetIdx())
        ratoms.add(a.GetIdx())
    for bond in a.GetBonds():
        if bond.IsInRing(): # or is_hbond(bond):

           # print('in ring: {}'.format(bond))
            beg = bond.GetBgn()
            end = bond.GetEnd()
            beg_index = beg.GetIdx()
            end_index = end.GetIdx()
            if bond.GetIdx() not in rbonds_l:
                rbonds_l.append(bond.GetIdx())
                rbonds.add(bond.GetIdx())
            if beg_index not in ratoms_l:
                ratoms_l.append(beg_index)
                ratoms.add(beg_index)
                find_ring(beg, mol, fgroup_tagged, ratoms_l, rbonds_l, rot_bond, next_bond)
            if end_index not in ratoms_l:
                ratoms_l.append(end_index)
                ratoms.add(end_index)
                find_ring(end, mol, fgroup_tagged, ratoms_l, rbonds_l, rot_bond, next_bond)
        else:
            # check for functional group
            #fgroup = is_fgroup(fgroup_tagged, atom=bond.GetBgn())
            #if not fgroup:
            #    fgroup = is_fgroup(fgroup_tagged, atom=bond.GetEnd())
            #if fgroup:
            #    for fa_idx in fgroup[0]:
            #        ratoms_l.append(fa_idx)

            #    for fb_idx in fgroup[-1]:
            #        rbonds_l.append(fb_idx)
            # Check if orhto
            if is_ortho(bond, rot_bond, next_bond):
                beg_idx = bond.GetBgn().GetIdx()
                end_idx = bond.GetEnd().GetIdx()
                if bond.GetIdx() not in rbonds_l:
                    rbonds_l.append(bond.GetIdx())
                    rbonds.add(bond.GetIdx())
                if beg_idx not in ratoms_l:
                    ratoms_l.append(beg_idx)
                    ratoms.add(beg_idx)
                if end_idx not in ratoms_l:
                    ratoms_l.append(end_idx)
                    ratoms.add(end_idx)
                # Check for functional group
                fgroup = is_fgroup(fgroup_tagged, atom=bond.GetBgn())
                if not fgroup:
                    fgroup = is_fgroup(fgroup_tagged, atom=bond.GetEnd())
                if fgroup:
                   # print(fgroup)
                    for fa_idx in fgroup[0]:
                        ratoms_l.append(fa_idx)
                        ratoms.add(fa_idx)
                        # Check if atom is in a ring. Add entire ring
                        # fg_atom = mol.GetAtom(oechem.OEHasAtomIdx(fa_idx))
                        # if fg_atom.IsInRing():
                        #     print('fg atom in ring: {}'.format(fg_atom))
                        #
                        #     find_ring(fg_atom, mol, fgroup_tagged, ratoms_l, rbonds_l, rot_bond, next_bond)

                    for fb_idx in fgroup[-1]:
                        fg_bond = mol.GetBond(oechem.OEHasBondIdx(fb_idx))
                        if fg_bond.IsInRing():
                            #print('fg bond in ring: {}'.format(fg_bond))
                            #print(fg_ring)
                            beg = fg_bond.GetBgn()
                            end = fg_bond.GetEnd()
                            beg_index = beg.GetIdx()
                            end_index = end.GetIdx()
                            #if fg_bond.GetIdx() not in rbonds_l:
                            #    rbonds_l.append(fg_bond.GetIdx())
                            if beg_index not in fg_ring:
                                fg_ring.append(beg_index)
                                find_ring(beg, mol, fgroup_tagged, ratoms_l=[], rbonds_l=[], rot_bond=rot_bond, next_bond=next_bond, fg_ring=fg_ring)
                            if end_index not in fg_ring:
                                fg_ring.append(end_index)
                                find_ring(end, mol, fgroup_tagged, ratoms_l=[], rbonds_l=[], rot_bond=rot_bond, next_bond=next_bond, fg_ring=fg_ring)
                        else:
                            rbonds_l.append(fb_idx)
                            rbonds.add(fb_idx)

            # Check if non rotatable
            if not bond.IsRotor():
                # Keep all constituents that are not rotatable (halogens, O, ch3..)
                beg_idx = bond.GetBgn().GetIdx()
                end_idx = bond.GetEnd().GetIdx()
                if bond.GetIdx() not in rbonds_l:
                    rbonds_l.append(bond.GetIdx())
                    rbonds.add(bond.GetIdx())
                if beg_idx not in ratoms_l:
                    ratoms_l.append(beg_idx)
                    ratoms.add(beg_idx)
                if end_idx not in ratoms_l:
                    ratoms_l.append(end_idx)
                    ratoms.add(end_idx)

    return ratoms, rbonds


def to_AtomBondSet(mol, atoms, bonds):
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


def is_ortho(bond, rot_bond, next_bond):
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
    beg_b = bond.GetBgn()
    end_b = bond.GetEnd()
    beg_r = rot_bond.GetBgn()
    end_r = rot_bond.GetEnd()
    for b in beg_b.GetBonds():
        bond_attached.add(b.GetIdx())
    for b in end_b.GetBonds():
        bond_attached.add(b.GetIdx())
    for b in beg_r.GetBonds():
        rot_attached.add(b.GetIdx())
    for b in end_r.GetBonds():
        rot_attached.add(b.GetIdx())
    if not next_bond.IsInRing():
        next_attached = set()
        beg_n = next_bond.GetBgn()
        end_n = next_bond.GetEnd()
        for b in beg_n.GetBonds():
            next_attached.add(b.GetIdx())
        for b in end_n.GetBonds():
            next_attached.add(b.GetIdx())

    intersection = (bond_attached & rot_attached)
    if not bool(intersection) and not next_bond.IsInRing():
        # Check if it's ortho to next bond
        intersection = (bond_attached & next_attached)
    return bool(intersection)


def build_frag(bond, mol, tagged_fgroups, tagged_rings):


    atoms = set()
    bonds = set()
    b_idx = bond.GetIdx()
    bonds.add(b_idx)
    beg = bond.GetBgn()
    end = bond.GetEnd()
    beg_idx = beg.GetIdx()
    end_idx = end.GetIdx()

    atoms.add(beg_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol, bond, beg, tagged_fgroups, tagged_rings, atoms_2=set(), bonds_2=set())
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    atoms.add(end_idx)
    atoms_nb, bonds_nb = iterate_nbratoms(mol, bond, end, tagged_fgroups, tagged_rings, atoms_2=set(), bonds_2=set())
    atoms = atoms.union(atoms_nb)
    bonds = bonds.union(bonds_nb)

    return atoms, bonds








# def build_frag(bond, mol, fgroup_tagged, atoms=set(), bonds=set(), i=0):
#     """
#
#     Parameters
#     ----------
#     bond
#     mol
#     fgroup_tagged
#     atoms
#     bonds
#     i
#
#     Returns
#     -------
#
#     """
#     b_idx = bond.GetIdx()
#     if b_idx not in bonds:
#         bonds.add(b_idx)
#     beg = bond.GetBgn()
#     end = bond.GetEnd()
#     beg_idx = beg.GetIdx()
#     end_idx = end.GetIdx()
#     #if beg_idx not in atoms:
#     atoms.add(beg_idx)
#     print('beg atom: {}'.format(beg))
#     iterate_nbratom(beg, mol, fgroup_tagged, bond, beg_idx, end_idx, atoms, bonds, i)
#     print('atoms after beg iteration: {}'.format(atoms))
#     atoms.add(end_idx)
#     print('end atom: {}'.format(end))
#     iterate_nbratom(end, mol, fgroup_tagged, bond, end_idx, beg_idx, atoms, bonds, i)
#     print('atoms after end iteration: {}'.format(atoms))
#     return atoms, bonds


# def iterate_nbratom(atom, mol, fgroup_tagged, bond, s_idx, o_idx, atoms, bonds, i, check_ortho=False):
#     """
#
#     Parameters
#     ----------
#     atom
#     mol
#     fgroup_tagged
#     bond
#     s_idx
#     o_idx
#     atoms
#     bonds
#     i
#
#     Returns
#     -------
#
#     """
#     for a in atom.GetAtoms():
#         a_idx = a.GetIdx()
#         next_bond = mol.GetBond(a, atom)
#         nb_idx = next_bond.GetIdx()
#         if next_bond.GetData('WibergBondOrder') <= 1.2 and i > 0:
#             continue
#         atoms.add(a_idx)
#         if nb_idx != bond.GetIdx() and nb_idx not in bonds:
#             bonds.add(nb_idx)
#             if a_idx != o_idx:
#                 if a.IsInRing():
#                     r_atoms, r_bonds = find_ring(a, mol, fgroup_tagged, ratoms_l=[], rbonds_l=[], rot_bond=bond,
#                                                  next_bond=next_bond, fg_ring=[], ratoms=set(), rbonds=set())
#                     for ra_idx in r_atoms:
#                         atoms.add(ra_idx)
#                     for rb_idx in r_bonds:
#                         bonds.add(rb_idx)
#                 fgroup = is_fgroup(fgroup_tagged, atom=a)
#                 if fgroup:
#                     for fa_idx in fgroup[0]:
#                         atoms.add(fa_idx)
#                         # Check if atom is in a ring. If it's in a ring, include the ring - we don't want to fragment rings
#                         fg_atom = mol.GetAtom(oechem.OEHasAtomIdx(fa_idx))
#                         if fg_atom.IsInRing():
#                             print('Atom in fg and ring: {}'.format(fa_idx))
#                             print('atoms before added ring: {}'.format(atoms))
#                             print('bonds before added ring: {}'.format(bonds))
#
#                             r_atoms, r_bonds = find_ring(a, mol, fgroup_tagged, ratoms_l=[], rbonds_l=[], rot_bond=bond,
#                                                          next_bond=next_bond, fg_ring=[], ratoms=set(), rbonds=set())
#                             for ra_idx in r_atoms:
#                                 atoms.add(ra_idx)
#                             for rb_idx in r_bonds:
#                                 bonds.add(rb_idx)
#                             print('atoms after added ring: {}'.format(atoms))
#                             print('bonds after added ring: {}'.format(bonds))
#
#                     for fb_idx in fgroup[-1]:
#                         bonds.add(fb_idx)
#
#                 else:
#                     for n_atom in a.GetAtoms():
#                         n_idx = n_atom.GetIdx()
#                         if n_idx not in (a_idx, s_idx):
#                             nn_bond = mol.GetBond(n_atom, a)
#                             if nn_bond.GetData('WibergBondOrder') >=1.2:
#                                 bonds.add(nn_bond.GetIdx())
#                                 atoms.add(n_idx)
#                                 i +=1
#                                 build_frag(nn_bond, mol, fgroup_tagged, atoms, bonds, i=i)



def iterate_nbratoms(mol, rotor_bond,  atom, fgroup_tagged, tagged_rings, atoms_2=set(), bonds_2=set()):
    """

    Parameters
    ----------
    mol
    atom
    fgroup_tagged
    tagged_rings
    rotor_bond

    Returns
    -------

    """
    #atoms_2 = set()
    #bonds_2 = set()
    for a in atom.GetAtoms():
        a_idx = a.GetIdx()
        #if a_idx in atoms_2:
        #    continue
        next_bond = mol.GetBond(a, atom)
        nb_idx = next_bond.GetIdx()
        atoms_2.add(a_idx)
        if nb_idx in bonds_2:
            FGROUP_RING = False
            # Check Data for both atoms in bond.
            try:
                a.GetData('ringsystem')
                a.GetData('fgroup')
                FGROUP_RING = True
            except ValueError:
                try:
                    atom.GetData('ringsystem')
                    atom.GetData('fgroup')
                    FGROUP_RING = True
                except ValueError:
                    pass
            if FGROUP_RING:
                # Add ring and then continue?
                continue


            #print('Already in list: {}'.format(next_bond))
            # if a.IsInRing() or atom.IsInRing():
            #     print('one of the atoms are in a ring')
            #     fgroup = is_fgroup(fgroup_tagged, atom=a)
            #     if not fgroup:
            #         fgroup = is_fgroup(fgroup_tagged, atom= a)
            #     if fgroup:
            #         print('one atom is in an fgroup')
            #         # Check that it's not the same atom that's in the fgroup and ring. I should probably write a function
            #         # That checks for one atom in bond is in a ring and another is in an fgroup and there's no resonance.
            # # Check if one atom in bond is in a ring and the other in an fgroup. Then we need to add the ring.
            # fgroup = is_fgroup(fgroup_tagged, bond=next_bond)
            # if not next_bond.IsInRing() and not fgroup:
            #     continue
            #     # Check if it's also in fgroup
            #
            #     # if not fgroup:
                #     continue
                # #print('In list but not in ring: {}'.format(next_bond))
                #continue
                # ring_idx = next_bond.GetData('ringsystem')
                # ratoms, rbonds = tagged_rings[ring_idx]
                #
                # if nb_idx not in rbonds:
                # # add to new list?
                # print('alredy in list and not in ring: {}'.format(next_bond))
                # continue
        bonds_2.add(nb_idx)
        if a.IsInRing():
            print('In ring: {}'.format(a))
            ring_idx = a.GetData('ringsystem')
            ratoms, rbonds = tagged_rings[ring_idx]
            atoms_2 = atoms_2.union(ratoms)
            bonds_2 = bonds_2.union(rbonds)
            # Find non-rotatable sustituents
            rs_atoms, rs_bonds = ring_substiuents(mol, next_bond, rotor_bond, tagged_rings, ring_idx, fgroup_tagged)
            atoms_2 = atoms_2.union(rs_atoms)
            bonds_2 = bonds_2.union(rs_bonds)
        fgroup = is_fgroup(fgroup_tagged, atom=atom)
        if fgroup:
            atoms_2 = atoms_2.union(fgroup[0])
            bonds_2 = bonds_2.union(fgroup[-1])
            # if something is in a ring - have a flag? Then use that to continue iterating and change the flag
        #Check for resonance here?
        for nb_a in a.GetAtoms():
            nn_bond = mol.GetBond(a, nb_a)
            if nn_bond.GetData('WibergBondOrder') > 1.2:
                #print(next_bond)
                atoms_2.add(nb_a.GetIdx())
                #bonds_2.add(nn_bond.GetIdx())
                iterate_nbratoms(mol, nn_bond, nb_a, fgroup_tagged, tagged_rings, atoms_2, bonds_2)
    return atoms_2, bonds_2

def ring_substiuents(mol, bond, rotor_bond, tagged_rings, ring_idx, fgroup_tagged):
    rs_atoms = set()
    rs_bonds = set()
    r_atoms, r_bonds = tagged_rings[ring_idx]
    for a_idx in r_atoms:
        atom = mol.GetAtom(oechem.OEHasAtomIdx(a_idx))
        for a in atom.GetAtoms():
            if a.GetIdx() in rs_atoms:
                continue
            fgroup = False
            if not a.IsInRing():
                rs_bond = mol.GetBond(atom, a)
                if not rs_bond.IsRotor():
                    rs_atoms.add(a.GetIdx())
                    #bond = mol.GetBond(atom, a)
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = is_fgroup(fgroup_tagged, atom=a)
                elif is_ortho(rs_bond, rotor_bond, bond):
                    # Keep bond and attached atom.
                    rs_atoms.add(a.GetIdx())
                    rs_bonds.add(rs_bond.GetIdx())
                    # Check for functional group
                    fgroup = is_fgroup(fgroup_tagged, atom=a)
                if fgroup:
                    rs_atoms = rs_atoms.union(fgroup[0])
                    rs_bonds = rs_bonds.union(fgroup[-1])
                #print('functional_group: {}'.format(fgroup))
                #print('ring idex: {}'.format(tagged_rings[ring_idx]))
                    # Check if ortho
                    #is
                    #if bond.IsInRing():
                        # bond is inside ring so we need to check rotor_bond
                    #    bond = rotor_bond


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


def GetFragmentAtomBondSetCombinations(fraglist, MAX_ROTORS=3, MIN_ROTORS=1):
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
