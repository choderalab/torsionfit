import netCDF4 as netcdf
from pymc.database import base
from pymc import six
import pymc
import numpy as np
import os
import sys
import warnings
import traceback
try:
    import cPickle as pickle
except ImportError:
    import pickle
import codecs

__all__ = ['Trace', 'Database']


class Trace(base.Trace):

    """netCDF4 Trace class"""

    def _initialize(self):

        if self._getfunc is None:
            self._getfunc = self.db.model._funs_to_tally[self.name]

    def tally(self, index, chain=-1):
        """ Add current value to chain
        :param index: tall index
        :param chain: int
        """

        value = self._getfunc()
        self.db.ncfile['Chain#%d' % chain].variables[self.name][index] = value

    def gettrace(self, burn=0, thin=1, chain=-1, slicing=None):
        """Return the trace (last by default).
        :Parameters:
        burn : integer
          The number of transient steps to skip.
        thin : integer
          Keep one in thin.
        chain : integer
          The index of the chain to fetch. If None, return all chains. The
          default is to return the last chain.
        slicing : slice object
          A slice overriding burn and thin assignement.
        """

        if chain is not None:
            if chain == 0:
                chain = 0
            elif chain == -1:
                chain = range(self.db.chains)[chain]
            arr = self.db.ncfile['Chain#%d' % chain][self.name][:]
        elif chain is None:
            arr = self.db.ncfile['Chain#0'][self.name][:]
            for i in range(self.db.chains)[1:]:
                arr = np.append(arr, self.db.ncfile['Chain#%d' % i][self.name][:])

        if slicing is None:
            slicing = slice(burn, None, thin)

        return arr[slicing].squeeze()

    __call__ = gettrace

    def __getitem__(self, i):
        """Return the trace corresponding to item (or slice) i for the chain
        defined by self._chain.
        """
        if isinstance(i, slice):
            return self.gettrace(slicing=i, chain=self._chain)
        else:
            return self.gettrace(slicing=slice(i, i + 1), chain=self._chain)

    def length(self, chain=-1):
        """
        :param chain: int or None
            The chain index. If None return all chains. default is last chain
        :return: int length of chain
        """

        return len(self.gettrace(chain=chain))


class Database(base.Database):

    """netCDF4 database"""

    def __init__(self, dbname, dbmode='w'):
        """ Open or create netCDF4 database
        creates a group for every chain. Creates ncvar for every parameter with unlimited dimensions (and dimension of
        pymc parameter) to allow for unlimited iterations. The state of the chain is saved as a pickled string in ncvar
        state for each chain group.

        :param dbname: string
            name of database file
        :param dbmode: {'a' or 'w'}
            File mode. Use 'a' to append and 'w' to overwrite an existing file

        """
        self.__name__ = 'netcdf'
        self.dbname = dbname
        self.mode = dbmode
        self.trace_names = []
        self._traces = {} # dictionary of trace objects
        self.__Trace__ = Trace
        self._chains ={} # dictionary of states for each chain

        # check if database exists
        db_exists = os.path.exists(self.dbname)

        if db_exists and self.mode == 'w':
            print('overwriting file %s' % self.dbname)
            # overwrite
            os.remove(self.dbname)

        # open netCDF4 file
        self.ncfile = netcdf.Dataset(dbname, self.mode, version= 'NETCDF4')
        if 'state' not in self.ncfile.dimensions:
            self.ncfile.createDimension('state', 1)

        # Set global attributes
        if self.mode == 'w':
            setattr(self.ncfile, 'title', self.dbname)

        # Assign self.chain to last chain
        self.chains = len(self.ncfile.groups)
        for (i, group) in enumerate(self.ncfile.groups):
            self._chains[i] = pickle.loads(codecs.decode(self.ncfile[group]['state'][0].encode(), "base64"))

        # Load existing data
        existing_chains = self.ncfile.groups
        for chain in existing_chains:
            names = []
            for var in self.ncfile[chain].variables:
                if var not in self._traces:
                    self._traces[var] = Trace(name=var, db=self)

                names.append(var)
            self.trace_names.append(names)

    def _initialize(self, funs_to_tally, length=None):

        for name, fun in six.iteritems(funs_to_tally):
            if name not in self._traces:
                self._traces[name] = self.__Trace__(name=name, getfunc=fun, db=self)
            # if db is loaded from disk, it might not have its tallied step method
            self._traces[name]._initialize()

        i = self.chains
        self.ncfile.createGroup("Chain#%d" % i)

        # set dimensions
        if 'nsamples' not in self.ncfile['Chain#%d' %i].dimensions:
            self.ncfile['Chain#%d' % i].createDimension('nsamples', 0)  # unlimited number of iterations

        # sanity check that nsamples is unlimited
        if self.ncfile['Chain#%d' % i].dimensions['nsamples'].isunlimited():
            pass

        # create variables for pymc variables
        for name, fun in six.iteritems(funs_to_tally):
            if not np.asarray(fun()).shape == () and name not in self.ncfile['Chain#%d' % i].variables:
                # create ncvar with dimensions of pymc parameter array and nsamples
                self.ncfile['Chain#%d' % i].createDimension(name, np.asarray(fun()).shape[0])
                self.ncfile['Chain#%d' % i].createVariable(name, np.asarray(fun()).dtype.str, ('nsamples', name))
            elif name not in self.ncfile['Chain#%d' % i].variables:
                # all other ncvar than only need nsamples dimension
                self.ncfile['Chain#%d' % i].createVariable(name, np.asarray(fun()).dtype.str, ('nsamples',))

        if len(self.trace_names) < len(self.ncfile.groups):
            try:
                self.trace_names.append(list(self.ncfile['Chain#%d' % self.chains].variables))
            except IndexError:
                self.trace_names.append(list(funs_to_tally.keys()))
        self.tally_index = len(self.ncfile['Chain#%d' % self.chains].dimensions['nsamples'])
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
            if len(names) > 0:
                print("Some objects from the database have not been assigned a getfunc", names)

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
                self._traces[name].tally(self.tally_index, chain)
            except:
                cls, inst, tb = sys.exc_info()
                warnings.warn("""
Error tallying %s, will not try to tally it again this chain.
Did you make all the same variables and step methods tallyable
as were tallyable last time you used the database file?
Error:
%s""" % (name, ''.join(traceback.format_exception(cls, inst, tb))))
                self.trace_names[chain].remove(name)

        self.tally_index += 1

    def savestate(self, state, chain=-1):

        self._state_ = state

        # pickle state
        chain = range(self.chains)[chain]
        state_pickle = codecs.encode(pickle.dumps(state), "base64").decode()

        # save pickled state in ncvar in group for current chain
        if 'state' not in self.ncfile['Chain#%d' % chain].variables:
            self.ncfile['Chain#%d' % chain].createVariable('state', str, ('state',))
            self.ncfile['Chain#%d' % chain]['state'][0] = state_pickle

    def getstate(self, chain=-1):

        if self.chains == 0:
            chain = 0
        else:
            chain = range(self.chains)[chain]

        if len(self._chains) == 0:
            return {}
        else:
            return self._chains[chain]

    def close(self):
        self.ncfile.close()


def load(dbname, dbmode='a'):
    """ Load an existing netcdf database

    :param dbname: name of netcdf file to open
    :param dbmode: 'a': append, 'r': read-only
    :return: database
    """
    if dbmode == 'w':
        raise AttributeError("dbmode='w' not allowed for load")

    db = Database(dbname, dbmode=dbmode)

    return db




