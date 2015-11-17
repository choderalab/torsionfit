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

    def tally(self, index, chain=-1):
        """ Add current value to chain
        :param chain:
        :return:
        """

        value = self._getfunc()
        index = self.model.db._current_iter
        self.db.ncfile.variables[self.name][index] = value




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


    def _initialize(self, funs_to_tally, length=None):

        for name, fun in six.iteritems(funs_to_tally):
            if name not in self._traces:
                self._traces[name] = self.__Trace__(name=name, getfunc=fun, db=self)

        # set dimensions
        self.ncfile.createDimension('nsamples', 0) # unlimited number of iterations

        for name, fun in six.iteritems(funs_to_tally):
            if not np.asarray(fun()).shape == ():
                self.ncfile.createDimension(name, np.asarray(fun()).shape[0])
                self.ncfile.createVariable(name, np.asarray(fun()).dtype.str, ('nsamples', name))
            else:
                self.ncfile.createVariable(name, np.asarray(fun()).dtype.str, ('nsamples',))


        self.trace_names.append(list(funs_to_tally.keys()))

        self.chains += 1


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

    def tally(self, chain=-1):
        """Append the current value of all tallyable object.
       :Parameters:
       chain : int
         The index of the chain to append the values to. By default, the values
         are appended to the last chain.
        """

        chain = range(self.chains)[chain]
        for name in self.trace_names[chain]:
            try:
                self._traces[name].tally(chain)
            except:
                cls, inst, tb = sys.exc_info()
                warnings.warn("""
Error tallying %s, will not try to tally it again this chain.
Did you make all the samevariables and step methods tallyable
as were tallyable last time you used the database file?
Error:
%s""" % (name, ''.join(traceback.format_exception(cls, inst, tb))))
                self.trace_names[chain].remove(name)



