from numpy import ma, squeeze
from pupynere import netcdf_file, _unmaskandscale
from dap.client import open as open_remote
from dap.dtypes import ArrayType, GridType, typemap

_typecodes = dict([[_v,_k] for _k,_v in typemap.items()])

class _RemoteFile(object):
    """A NetCDF file reader. API is the same as Scientific.IO.NetCDF."""

    def __init__(self, file, maskandscale=False, cache=None,\
                 username=None, password=None, verbose=False):
        self._buffer = open_remote(file,cache=cache,\
                       username=username,password=password,verbose=verbose)
        self._maskandscale = maskandscale
        self._parse()

    def read(self, size=-1):
        """Alias for reading the file buffer."""
        return self._buffer.read(size)

    def _parse(self):
        """Initial parsing of the header."""
        # Read header info.
        self._dim_array()
        self._gatt_array()
        self._var_array()

    def _dim_array(self):
        """Read a dict with dimensions names and sizes."""
        self.dimensions = {}
        self._dims = []
        for k,d in self._buffer.iteritems():
            if (isinstance(d, ArrayType) or isinstance(d, GridType)) and len(d.shape) == 1 and k == d.dimensions[0]:
                name = k
                length = len(d)
                self.dimensions[name] = length
                self._dims.append(name)  # preserve dim order

    def _gatt_array(self):
        """Read global attributes."""
        self.__dict__.update(self._buffer.attributes)

    def _var_array(self):
        """Read all variables."""
        # Read variables.
        self.variables = {}
        for k,d in self._buffer.iteritems():
            if isinstance(d, GridType) or isinstance(d, ArrayType):
                name = k
                self.variables[name] = _RemoteVariable(d,self._maskandscale)

    def close(self):
        # this is a no-op provided for compatibility
        pass

class _RemoteVariable(object):
    def __init__(self, var, maskandscale):
        self._var = var
        self._maskandscale = maskandscale
        self.dtype = var.type
        self.shape = var.shape
        self.dimensions = var.dimensions
        self.__dict__.update(var.attributes)

    def __getitem__(self, index):
        datout = squeeze(self._var.__getitem__(index))
        # automatically
        # - remove singleton dimensions
        # - create a masked array using missing_value or _FillValue attribute
        # - apply scale_factor and add_offset to packed integer data
        if self._maskandscale:
            return _unmaskandscale(self,datout)
        else:
            return datout

    def __len__(self):
       return self.shape[0]

    def typecode(self):
        return _typecodes[self.dtype]
