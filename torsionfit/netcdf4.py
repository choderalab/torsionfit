import netCDF4 as netcdf
from pymc.database import base
from pymc import six
import pymc
import numpy as np

__all__ = ['Trace', 'Database']

class Trace(base.Trace):

    """netCDF4 Trace class"""

    def _initialize(self, chain, length):
        """Create a netCDF4 variable.
        """

        if self._getfunc is None:
            self._getfunc = self.db.model._funs_to_tally[self.name]




        # Determine size
        # try:
        #     self._shape = np.shape(self._getfunc())
        # except TypeError:
        #     self._shape = None
        #
        # self._vstr = ', '.join(var_str(self._shape))
        #
        # # If the table already exists, exit now.
        # if chain != 0:
        #     return
        #
        # # Create the variable name strings.
        # vstr = ', '.join(v + ' FLOAT' for v in var_str(self._shape))
        # query = """CREATE TABLE IF NOT EXISTS [%s]
        #              (recid INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
        #               trace  int(5), %s)""" % (self.name, vstr)
        # self.db.cur.execute(query)


class Database(base.Database):

    """netCDF4 database"""

    def __init__(self, dbname, dbmode='w'):
        """ Open or create netCDF4 database

        :param dbname: string
            name of database file
        :param dbmode: {'a' or 'w'}
            File mode. Use 'a' to append and 'w' to overwrite an existing file

        """
        self.__name__ = 'netcdf'
        self.dbname = dbname
        self.trace_names = []
        self._traces = {} # dictionary of trace objects
        self.__Trace__ = Trace

        # open netCDF4 file for writing
        self.ncfile = netcdf.Dataset(dbname, dbmode, version= 'NETCDF4')

        # Set global attributes
        setattr(self.ncfile, 'title', self.dbname)

        # # create netCDF4 variable
        # for name, fun in six.iteritems()
        # self.ncvar = self.ncfile.createVariable(self.name, dtype, ('nsamples', self._getfunc()))

    def _initialize(self, funs_to_tally, length=None):

        # set dimensions
        self.ncfile.createDimension('nsamples', 0) # unlimited number of iterations

        for name, fun in six.iteritems(funs_to_tally):
            self.ncfile.createVariable(name, np.asarray(fun()).dtype.str, ('nsamples',))


    def connect_model(self, model):
        """Link the Database to the Model instance.
        In case a new database is created from scratch, ``connect_model``
        creates Trace objects for all tallyable pymc objects defined in
        `model`.
        If the database is being loaded from an existing file, ``connect_model``
        restore the objects trace to their stored value.
        :Parameters:
        model : pymc.Model instance
          An instance holding the pymc objects defining a statistical
          model (stochastics, deterministics, data, ...)
        """
        # Changed this to allow non-Model models. -AP
        # We could also remove it altogether. -DH
        if isinstance(model, pymc.Model):
            self.model = model
        else:
            raise AttributeError('Not a Model instance.')

        # Restore the state of the Model from an existing Database.
        # The `load` method will have already created the Trace objects.
        if hasattr(self, '_state_'):
            names = set()
            for morenames in self.trace_names:
                names.update(morenames)
            for name, fun in six.iteritems(model._funs_to_tally):
                if name in self._traces:
                    self._traces[name]._getfunc = fun
                    names.discard(name)
            # if len(names) > 0:
            # print_("Some objects from the database have not been assigned a
            # getfunc", names)

        # Create a fresh new state.
        # We will be able to remove this when we deprecate traces on objects.
        else:
            for name, fun in six.iteritems(model._funs_to_tally):
                if name not in self._traces:
                    self._traces[name] = self.__Trace__(name=name, getfunc=fun, db=self)



