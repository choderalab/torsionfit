from pymol import stored, cmd
import os
import errno


def torsion_drive(atom1, atom2, atom3, atom4, interval, selection, path, mol_name,):
    """
    This function generates input pdbs of dihedral angles selected of intervals specified with interval
    :param atom1:
    :param atom2:
    :param atom3:
    :param atom4:
    :param interval:
    :param selection
    :param path:
    :param mole_name:
    :return:
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