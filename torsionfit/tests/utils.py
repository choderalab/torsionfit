import os
from pkg_resources import resource_filename


def get_fun(name):
    """Get the full path to one of the reference files shipped for testing

        These files are in torsionfit/testing/reference

    :param
        name: str
            Name of file to load

    :returns
        fn : str
            full path to file
    """

    fn = resource_filename('torsionfit', os.path.join('tests', 'reference', name))

    if not os.path.exists(fn):
        raise ValueError('%s does not exist. If you just added it you will have to re install' % fn)

    return fn


def get_files(fnames):
    """Get the full path to a list of reference files shipped for testing
        These file are in torsionfit/testing/reference
    :param
        list: list
             list of files
    :return:
        fnl: list of str
            full path to all files in the list
    """
    fnl = []
    for name in fnames:
        fn = resource_filename('torsionfit', os.path.join('tests', name))
        if not os.path.exists(fn):
            raise ValueError('%s does not exist. If you just added it you will have to re install' % fn)
        fnl.append(fn)

    return fnl
