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


