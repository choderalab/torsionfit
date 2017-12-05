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

def ToPdf(mol, oname, frags, fragcombs):
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

    DepictMoleculeWithFragmentCombinations(report, mol, frags, fragcombs, opts)

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
        G.add_edge(bond.GetBgnIdx(), bond.GetEndIdx(), weight=bond.GetData("WibergBondOrder"), index=bond.GetIdx())
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
            if G.edge[node][node2]['weight'] < bondOrderThreshold and G.degree(node2) >1:
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


def SmilesToFragments(smiles, bondOrderThreshold=1.2, chargesMol=True):
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


def DepictMoleculeWithFragmentCombinations(report, mol, frags, fragcombs, opts):

    stag = "fragment idx"
    itag = oechem.OEGetTag(stag)
    for fidx, frag in enumerate(frags):
        for bond in frag.GetBonds():
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

    for frag in fragcombs:

        cell = report.NewCell()
        disp = oedepict.OE2DMolDisplay(mol, opts)

        fragatoms = oechem.OEIsAtomMember(frag.GetAtoms())
        fragbonds = oechem.OEIsBondMember(frag.GetBonds())

        notfragatoms = oechem.OENotAtom(fragatoms)
        notfragbonds = oechem.OENotBond(fragbonds)

        oedepict.OEAddHighlighting(disp, fadehighlight, notfragatoms, notfragbonds)
        oegrapheme.OEAddGlyph(disp, bondglyph, fragbonds)

        oedepict.OERenderMolecule(cell, disp)

    # depict original fragmentation in each header

    cellwidth, cellheight = report.GetHeaderWidth(), report.GetHeaderHeight()
    opts.SetDimensions(cellwidth, cellheight, oedepict.OEScale_AutoScale)
    opts.SetAtomColorStyle(oedepict.OEAtomColorStyle_WhiteMonochrome)

    bondlabel = LabelBondOrder()
    opts.SetBondPropertyFunctor(bondlabel)
    disp = oedepict.OE2DMolDisplay(mol, opts)
    oegrapheme.OEAddGlyph(disp, bondglyph, oechem.IsTrueBond())

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
