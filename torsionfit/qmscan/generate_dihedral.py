from pymol import stored, cmd
import os
import errno


def torsion_drive(atom1, atom2, atom3, atom4, interval, selection, path, mol_name,):
    """
    This function generates input pdbs of dihedral angles selected of intervals specified with interval
    :param atom1: name of atom 1 of dihedral
    :param atom2: name of atom 2 of dihedral
    :param atom3: name of atom 3 of dihedral
    :param atom4: name of atom 4 of dihedral
    :param interval: int or float (in degrees) of intervals to generate torsion scan for
    :param selection: name of selection for molecule
    :param path: path to where pdb files should be saved
    :param mole_name: name of molecule to append to filenamen
    """

    atom1 = selection + " and name " + atom1
    atom2 = selection + " and name " + atom2
    atom3 = selection + " and name " + atom3
    atom4 = selection + " and name " + atom4
    for angle in range(0, 360 + int(interval), int(interval)):
        try:
            os.makedirs('%s/%i' % (path, angle))
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

        cmd.set_dihedral(atom1, atom2, atom3, atom4, angle)
        filename = '%s/%i/%s_%i.pdb' % (path, angle, mol_name, angle)
        cmd.save(filename, selection, 1)

cmd.extend("torsion_drive", torsion_drive)