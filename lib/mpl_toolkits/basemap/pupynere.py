"""NetCDF reader.

Pupynere implements a PUre PYthon NEtcdf REader.

Copyright (c) 2003-2006 Roberto De Almeida <rob@pydap.org>

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject
to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"


import struct
import itertools
import mmap

from numpy import ndarray, empty, array, ma, squeeze, zeros
import numpy

from dap.client import open as open_remote
from dap.dtypes import ArrayType, GridType, typemap

has_pynio = True
try:
    from PyNGL import nio
except ImportError:
    has_pynio = False

ABSENT       = '\x00' * 8
ZERO         = '\x00' * 4
NC_BYTE      = '\x00\x00\x00\x01' 
NC_CHAR      = '\x00\x00\x00\x02'
NC_SHORT     = '\x00\x00\x00\x03'
NC_INT       = '\x00\x00\x00\x04'
NC_FLOAT     = '\x00\x00\x00\x05'
NC_DOUBLE    = '\x00\x00\x00\x06'
NC_DIMENSION = '\x00\x00\x00\n'
NC_VARIABLE  = '\x00\x00\x00\x0b'
NC_ATTRIBUTE = '\x00\x00\x00\x0c'

_typecodes = dict([[_v,_k] for _k,_v in typemap.items()])
# default _FillValue for netcdf types (apply also to corresponding
# DAP types).
_default_fillvals = {'c':'\0',
                     'S':"",
                     'b':-127,
                     'B':-127,
                     'h':-32767,
                     'H':65535,
                     'i':-2147483647L,
                     'L':4294967295L,
                     'q':-2147483647L,
                     'f':9.9692099683868690e+36,
                     'd':9.9692099683868690e+36}

def NetCDFFile(file, maskandscale=True):
    """NetCDF File reader.  API is the same as Scientific.IO.NetCDF.
    If 'file' is a URL that starts with 'http', it is assumed
    to be a remote OPenDAP dataset, and the python dap client is used
    to retrieve the data. Only the OPenDAP Array and Grid data
    types are recognized.  If file does not start with 'http', it
    is assumed to be a local file.  If possible, the file will be read 
    with a pure python NetCDF reader, otherwise PyNIO 
    (http://www.pyngl.ucar.edu/Nio.shtml) will be used (if it is installed).
    PyNIO supports NetCDF version 4, GRIB1, GRIB2, HDF4 and HDFEOS2 files.
    Data read from OPenDAP and NetCDF version 3 datasets will 
    automatically be converted to masked arrays if the variable has either
    a 'missing_value' or '_FillValue' attribute, and some data points
    are equal to the value specified by that attribute.  In addition,
    variables stored as integers that have the 'scale_factor' and
    'add_offset' attribute will automatically be rescaled to floats when
    read. If PyNIO is used, neither of the automatic conversions will
    be performed.  To suppress these automatic conversions, set the
    maskandscale keyword to False. 
    """
    if file.startswith('http'):
        return _RemoteFile(file,maskandscale)
    else:
        # use pynio if it is installed and the file cannot
        # be read with the pure python netCDF reader.  This allows
        # netCDF version 4, GRIB1, GRIB2, HDF4 and HDFEOS files
        # to be read.
        if has_pynio:
            try:
                f = _LocalFile(file,maskandscale)
            except:
                f = nio.open_file(file)
        # otherwise, use the pupynere netCDF 3 pure python reader.
        # (will fail if file is not a netCDF version 3 file).
        else:
            f = _LocalFile(file,maskandscale)
        return f
 
def _maskandscale(var,datout):
    totalmask = zeros(datout.shape,numpy.bool)
    fillval = None
    if hasattr(var, 'missing_value') and (datout == var.missing_value).any():
        fillval = var.missing_value
        totalmask += datout==fillval
    if hasattr(var, '_FillValue') and (datout == var._FillValue).any():
        if fillval is None:
            fillval = var._FillValue
        totalmask += datout==var._FillValue
    elif (datout == _default_fillvals[var.typecode()]).any():
        if fillval is None:
            fillval = _default_fillvals[var.typecode()]
        totalmask += datout==_default_fillvals[var.dtype]
    if fillval is not None:
        datout = ma.masked_array(datout,mask=totalmask,fill_value=fillval)
    try:
        datout = var.scale_factor*datout + var.add_offset
    except:
        pass
    return datout

class _RemoteFile(object):
    """A NetCDF file reader. API is the same as Scientific.IO.NetCDF."""

    def __init__(self, file, maskandscale):
        self._buffer = open_remote(file)
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
            return _maskandscale(self,datout)
        else:
            return datout

    def typecode(self):
        return _typecodes[self.dtype]


class _LocalFile(object):
    """A NetCDF file reader. API is the same as Scientific.IO.NetCDF."""

    def __init__(self, file, maskandscale):
        self._buffer = open(file, 'rb')
        self._maskandscale = maskandscale
        self._parse()

    def read(self, size=-1):
        """Alias for reading the file buffer."""
        return self._buffer.read(size)

    def _parse(self):
        """Initial parsing of the header."""
        # Check magic bytes.
        assert self.read(3) == 'CDF'

        # Read version byte.
        byte = self.read(1)
        self.version_byte = struct.unpack('>b', byte)[0]

        # Read header info.
        self._numrecs()
        self._dim_array()
        self._gatt_array()
        self._var_array()

    def _numrecs(self):
        """Read number of records."""
        self._nrecs = self._unpack_int()

    def _dim_array(self):
        """Read a dict with dimensions names and sizes."""
        assert self.read(4) in [ZERO, NC_DIMENSION]
        count = self._unpack_int()

        self.dimensions = {}
        self._dims = []
        for dim in range(count):
            name = self._read_string()
            length = self._unpack_int()
            if length == 0: length = None # record dimension
            self.dimensions[name] = length
            self._dims.append(name)  # preserve dim order

    def _gatt_array(self):
        """Read global attributes."""
        self.attributes = self._att_array()

        # Update __dict__ for compatibility with S.IO.N
        self.__dict__.update(self.attributes)

    def _att_array(self):
        """Read a dict with attributes."""
        assert self.read(4) in [ZERO, NC_ATTRIBUTE]
        count = self._unpack_int()

        # Read attributes.
        attributes = {}
        for attribute in range(count):
            name = self._read_string()
            nc_type = self._unpack_int()
            n = self._unpack_int()

            # Read value for attributes.
            attributes[name] = self._read_values(n, nc_type)

        return attributes

    def _var_array(self):
        """Read all variables."""
        assert self.read(4) in [ZERO, NC_VARIABLE]

        # Read size of each record, in bytes.
        self._read_recsize()

        # Read variables.
        self.variables = {}
        count = self._unpack_int()
        for variable in range(count):
            name = self._read_string()
            self.variables[name] = self._read_var()

    def _read_recsize(self):
        """Read all variables and compute record bytes."""
        pos = self._buffer.tell()
        
        recsize = 0
        count = self._unpack_int()
        for variable in range(count):
            name = self._read_string()
            n = self._unpack_int()
            isrec = False
            for i in range(n):
                dimid = self._unpack_int()
                name = self._dims[dimid]
                dim = self.dimensions[name]
                if dim is None and i == 0:
                    isrec = True
            attributes = self._att_array()
            nc_type = self._unpack_int()
            vsize = self._unpack_int()
            begin = [self._unpack_int, self._unpack_int64][self.version_byte-1]()

            if isrec: recsize += vsize

        self._recsize = recsize
        self._buffer.seek(pos)

    def _read_var(self):
        dimensions = []
        shape = []
        n = self._unpack_int()
        isrec = False
        for i in range(n):
            dimid = self._unpack_int()
            name = self._dims[dimid]
            dimensions.append(name)
            dim = self.dimensions[name]
            if dim is None and i == 0:
                dim = self._nrecs
                isrec = True
            shape.append(dim)
        dimensions = tuple(dimensions)
        shape = tuple(shape)

        attributes = self._att_array()
        nc_type = self._unpack_int()
        vsize = self._unpack_int()
        
        # Read offset.
        begin = [self._unpack_int, self._unpack_int64][self.version_byte-1]()

        return _LocalVariable(self._buffer.fileno(), nc_type, vsize, begin, shape, dimensions, attributes, isrec, self._recsize, maskandscale=self._maskandscale)

    def _read_values(self, n, nc_type):
        bytes = [1, 1, 2, 4, 4, 8]
        typecodes = ['b', 'c', 'h', 'i', 'f', 'd']
        
        count = n * bytes[nc_type-1]
        values = self.read(count)
        padding = self.read((4 - (count % 4)) % 4)
        
        typecode = typecodes[nc_type-1]
        if nc_type != 2:  # not char 
            values = struct.unpack('>%s' % (typecode * n), values)
            values = array(values, dtype=typecode) 
        else:
            # Remove EOL terminator.
            if values.endswith('\x00'): values = values[:-1]

        return values

    def _unpack_int(self):
        return struct.unpack('>i', self.read(4))[0]
    _unpack_int32 = _unpack_int

    def _unpack_int64(self):
        return struct.unpack('>q', self.read(8))[0]

    def _read_string(self):
        count = struct.unpack('>i', self.read(4))[0]
        s = self.read(count)
        # Remove EOL terminator.
        if s.endswith('\x00'): s = s[:-1]
        padding = self.read((4 - (count % 4)) % 4)
        return s

    def close(self):
        self._buffer.close()


class _LocalVariable(object):
    def __init__(self, fileno, nc_type, vsize, begin, shape, dimensions, attributes, isrec=False, recsize=0, maskandscale=True):
        self._nc_type = nc_type
        self._vsize = vsize
        self._begin = begin
        self.shape = shape
        self.dimensions = dimensions
        self.attributes = attributes  # for ``dap.plugins.netcdf``
        self.__dict__.update(attributes)
        self._is_record = isrec
        self._maskandscale = maskandscale

        # Number of bytes and type.
        self._bytes = [1, 1, 2, 4, 4, 8][self._nc_type-1]
        type_ = ['i', 'S', 'i', 'i', 'f', 'f'][self._nc_type-1]
        dtype = '>%s%d' % (type_, self._bytes)
        bytes = self._begin + self._vsize 

        if isrec:
            # Record variables are not stored contiguosly on disk, so we 
            # need to create a separate array for each record.
            self.__array_data__ = empty(shape, dtype)
            bytes += (shape[0] - 1) * recsize
            for n in range(shape[0]):
                offset = self._begin + (n * recsize)
                mm = mmap.mmap(fileno, bytes, access=mmap.ACCESS_READ)
                self.__array_data__[n] = ndarray.__new__(ndarray, shape[1:], dtype=dtype, buffer=mm, offset=offset, order=0)
        else:
            # Create buffer and data.
            mm = mmap.mmap(fileno, bytes, access=mmap.ACCESS_READ)
            self.__array_data__ = ndarray.__new__(ndarray, shape, dtype=dtype, buffer=mm, offset=self._begin, order=0)

        # N-D array interface
        self.__array_interface__ = {'shape'  : shape,
                                    'typestr': dtype,
                                    'data'   : self.__array_data__,
                                    'version': 3,
                                   }

    def __getitem__(self, index):
        datout = squeeze(self.__array_data__.__getitem__(index))
        # automatically
        # - remove singleton dimensions
        # - create a masked array using missing_value or _FillValue attribute
        # - apply scale_factor and add_offset to packed integer data
        if self._maskandscale:
            return _maskandscale(self,datout)
        else:
            return datout

    def getValue(self):
        """For scalars."""
        return self.__array_data__.item()

    def typecode(self):
        return ['b', 'c', 'h', 'i', 'f', 'd'][self._nc_type-1]
