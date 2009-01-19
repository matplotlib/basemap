"""
NetCDF reader/writer module.

This module implements the Scientific.IO.NetCDF API to read and create
NetCDF files. The same API is also used in the PyNIO and pynetcdf
modules, allowing these modules to be used interchangebly when working
with NetCDF files. The major advantage of ``scipy.io.netcdf`` over other 
modules is that it doesn't require the code to be linked to the NetCDF
libraries as the other modules do.

The code is based on the `NetCDF file format specification
<http://www.unidata.ucar.edu/software/netcdf/guide_15.html>`_. A NetCDF
file is a self-describing binary format, with a header followed by
data. The header contains metadata describing dimensions, variables
and the position of the data in the file, so access can be done in an
efficient manner without loading unnecessary data into memory. We use
the ``mmap`` module to create Numpy arrays mapped to the data on disk,
for the same purpose.

The structure of a NetCDF file is as follows:

    C D F <VERSION BYTE> <NUMBER OF RECORDS>
    <DIMENSIONS> <GLOBAL ATTRIBUTES> <VARIABLES METADATA>
    <NON-RECORD DATA> <RECORD DATA>

Record data refers to data where the first axis can be expanded at
will. All record variables share a same dimension at the first axis,
and they are stored at the end of the file per record, ie

    A[0], B[0], ..., A[1], B[1], ..., etc,
    
so that new data can be appended to the file without changing its original
structure. Non-record data are padded to a 4n bytes boundary. Record data
are also padded, unless there is exactly one record variable in the file,
in which case the padding is dropped.  All data is stored in big endian
byte order.

The Scientific.IO.NetCDF API allows attributes to be added directly to
instances of ``netcdf_file`` and ``netcdf_variable``. To differentiate
between user-set attributes and instance attributes, user-set attributes
are automatically stored in the ``_attributes`` attribute by overloading
``__setattr__``. This is the reason why the code sometimes uses
``obj.__dict__['key'] = value``, instead of simply ``obj.key = value``;
otherwise the key would be inserted into userspace attributes.

To create a NetCDF file::

    >>> import time
    >>> f = netcdf_file('simple.nc', 'w')
    >>> f.history = 'Created for a test'
    >>> f.createDimension('time', 10)
    >>> time = f.createVariable('time', 'i', ('time',))
    >>> time[:] = range(10)
    >>> time.units = 'days since 2008-01-01'
    >>> f.close()

To read the NetCDF file we just created::

    >>> f = netcdf_file('simple.nc', 'r')
    >>> print f.history
    Created for a test
    >>> time = f.variables['time']
    >>> print time.units
    days since 2008-01-01
    >>> print time.shape
    (10,)
    >>> print time[-1]
    9
    >>> f.close()

TODO: properly implement ``_FillValue``.
"""

__all__ = ['netcdf_file', 'netcdf_variable']


from operator import mul
from mmap import mmap, ACCESS_READ

from numpy import fromstring, ndarray, dtype, empty, array, asarray, squeeze,\
                  zeros, ma
from numpy import little_endian as LITTLE_ENDIAN


ABSENT       = '\x00\x00\x00\x00\x00\x00\x00\x00' 
ZERO         = '\x00\x00\x00\x00'
NC_BYTE      = '\x00\x00\x00\x01'
NC_CHAR      = '\x00\x00\x00\x02'
NC_SHORT     = '\x00\x00\x00\x03'
NC_INT       = '\x00\x00\x00\x04'
NC_FLOAT     = '\x00\x00\x00\x05'
NC_DOUBLE    = '\x00\x00\x00\x06'
NC_DIMENSION = '\x00\x00\x00\n'
NC_VARIABLE  = '\x00\x00\x00\x0b'
NC_ATTRIBUTE = '\x00\x00\x00\x0c'


TYPEMAP = { NC_BYTE:   ('b', 1),
            NC_CHAR:   ('c', 1),
            NC_SHORT:  ('h', 2),
            NC_INT:    ('i', 4),
            NC_FLOAT:  ('f', 4),
            NC_DOUBLE: ('d', 8) }

REVERSE = { 'b': NC_BYTE,
            'c': NC_CHAR,
            'h': NC_SHORT,
            'i': NC_INT,
            'f': NC_FLOAT,
            'd': NC_DOUBLE,

            # these come from asarray(1).dtype.char and asarray('foo').dtype.char,
            # used when getting the types from generic attributes.
            'l': NC_INT,
            'S': NC_CHAR }


class netcdf_file(object):
    """
    A ``netcdf_file`` object has two standard attributes: ``dimensions`` and
    ``variables``. The values of both are dictionaries, mapping dimension
    names to their associated lengths and variable names to variables,
    respectively. Application programs should never modify these
    dictionaries.

    All other attributes correspond to global attributes defined in the
    NetCDF file. Global file attributes are created by assigning to an
    attribute of the ``netcdf_file`` object.

    """
    def __init__(self, filename, mode='r', mmap=True, version=1,\
                 maskandscale=False):
        self.filename = filename
        self.use_mmap = mmap
        self.version_byte = version
        self._maskandscale = maskandscale

        assert mode in 'rw', "Mode must be either 'r' or 'w'."
        self.mode = mode

        self.dimensions = {}
        self.variables = {}

        self._dims = []
        self._recs = 0
        self._recsize = 0

        self.fp = open(self.filename, '%sb' % mode)

        self._attributes = {}

        if mode is 'r':
            self._read()

    def __setattr__(self, attr, value):
        # Store user defined attributes in a separate dict,
        # so we can save them to file later.
        try:
            self._attributes[attr] = value
        except AttributeError:
            pass
        self.__dict__[attr] = value

    def close(self):
        if not self.fp.closed:
            try:
                self.flush()
            finally:
                self.fp.close()
    __del__ = close

    def createDimension(self, name, length):
        self.dimensions[name] = length
        self._dims.append(name)

    def createVariable(self, name, type, dimensions):
        shape = tuple([self.dimensions[dim] for dim in dimensions]) 
        shape_ = tuple([dim or 0 for dim in shape])  # replace None with 0 for numpy

        if isinstance(type, basestring): type = dtype(type)
        typecode, size = type.char, type.itemsize
        dtype_ = '>%s' % typecode
        if size > 1: dtype_ += str(size)

        data = empty(shape_, dtype=dtype_)
        self.variables[name] = netcdf_variable(data, typecode, shape,\
                dimensions, maskandscale=self._maskandscale)
        return self.variables[name]

    def flush(self):
        if self.mode is 'w':
            self._write()
    sync = flush

    def _write(self):
        self.fp.write('CDF')
        self.fp.write(array(self.version_byte, '>b').tostring())

        # Write headers and data.
        self._write_numrecs()
        self._write_dim_array()
        self._write_gatt_array()
        self._write_var_array()

    def _write_numrecs(self):
        # Get highest record count from all record variables.
        for var in self.variables.values():
            if var.isrec and len(var.data) > self._recs:
                self.__dict__['_recs'] = len(var.data)
        self._pack_int(self._recs)

    def _write_dim_array(self):
        if self.dimensions:
            self.fp.write(NC_DIMENSION)
            self._pack_int(len(self.dimensions))
            for name in self._dims:
                self._pack_string(name)
                length = self.dimensions[name]
                self._pack_int(length or 0)  # replace None with 0 for record dimension
        else:
            self.fp.write(ABSENT)

    def _write_gatt_array(self):
        self._write_att_array(self._attributes)

    def _write_att_array(self, attributes):
        if attributes:
            self.fp.write(NC_ATTRIBUTE)
            self._pack_int(len(attributes))
            for name, values in attributes.items():
                self._pack_string(name)
                self._write_values(values)
        else:
            self.fp.write(ABSENT)

    def _write_var_array(self):
        if self.variables:
            self.fp.write(NC_VARIABLE)
            self._pack_int(len(self.variables))

            # Sort variables non-recs first, then recs.
            variables = self.variables.items()
            variables.sort(key=lambda (k, v): v._shape and not v.isrec)
            variables.reverse()
            variables = [k for (k, v) in variables]

            # Set the metadata for all variables.
            for name in variables:
                self._write_var_metadata(name)
            # Now that we have the metadata, we know the vsize of
            # each record variable, so we can calculate recsize.
            self.__dict__['_recsize'] = sum([
                    var._vsize for var in self.variables.values()
                    if var.isrec])
            # Set the data for all variables.
            for name in variables:
                self._write_var_data(name)
        else:
            self.fp.write(ABSENT)

    def _write_var_metadata(self, name):
        var = self.variables[name]

        self._pack_string(name)
        self._pack_int(len(var.dimensions))
        for dimname in var.dimensions:
            dimid = self._dims.index(dimname)
            self._pack_int(dimid)

        self._write_att_array(var._attributes)

        nc_type = REVERSE[var.typecode()]
        self.fp.write(nc_type)

        if not var.isrec:
            vsize = var.data.size * var.data.itemsize
            vsize += -vsize % 4
        else:  # record variable
            try:
                vsize = var.data[0].size * var.data.itemsize
            except IndexError:
                vsize = 0
            rec_vars = len([var for var in self.variables.values()
                    if var.isrec])
            if rec_vars > 1:
                vsize += -vsize % 4
        self.variables[name].__dict__['_vsize'] = vsize
        self._pack_int(vsize)

        # Pack a bogus begin, and set the real value later.
        self.variables[name].__dict__['_begin'] = self.fp.tell()
        self._pack_begin(0)

    def _write_var_data(self, name):
        var = self.variables[name]
        
        # Set begin in file header.
        the_beguine = self.fp.tell()
        self.fp.seek(var._begin)
        self._pack_begin(the_beguine)
        self.fp.seek(the_beguine)

        # Write data.
        if not var.isrec:
            self.fp.write(var.data.tostring())    
            count = var.data.size * var.data.itemsize
            self.fp.write('0' * (var._vsize - count))
        else:  # record variable
            # Handle rec vars with shape[0] < nrecs.
            if self._recs > len(var.data):
                shape = (self._recs,) + var.data.shape[1:]
                var.data.resize(shape)

            pos0 = pos = self.fp.tell()
            for rec in var.data:
                # Apparently scalars cannot be converted to big endian. If we
                # try to convert a ``=i4`` scalar to, say, '>i4' the dtype
                # will remain as ``=i4``.
                if not rec.shape and (rec.dtype.byteorder == '<' or
                        (rec.dtype.byteorder == '=' and LITTLE_ENDIAN)):
                    rec = rec.byteswap()
                self.fp.write(rec.tostring())
                # Padding
                count = rec.size * rec.itemsize
                self.fp.write('0' * (var._vsize - count))
                pos += self._recsize
                self.fp.seek(pos)
            self.fp.seek(pos0 + var._vsize)

    def _write_values(self, values):
        values = asarray(values) 
        values = values.astype(values.dtype.newbyteorder('>'))

        nc_type = REVERSE[values.dtype.char]
        self.fp.write(nc_type)

        if values.dtype.char == 'S':
            nelems = values.itemsize
        else:
            nelems = values.size
        self._pack_int(nelems)

        if not values.shape and (values.dtype.byteorder == '<' or
                (values.dtype.byteorder == '=' and LITTLE_ENDIAN)):
            values = values.byteswap()
        self.fp.write(values.tostring())
        count = values.size * values.itemsize
        self.fp.write('0' * (-count % 4))  # pad

    def _read(self):
        # Check magic bytes and version
        magic = self.fp.read(3)
        assert magic == 'CDF', "Error: %s is not a valid NetCDF 3 file" % self.filename
        self.__dict__['version_byte'] = fromstring(self.fp.read(1), '>b')[0]

        # Read file headers and set data.
        self._read_numrecs()
        self._read_dim_array()
        self._read_gatt_array()
        self._read_var_array()

    def _read_numrecs(self):
        self.__dict__['_recs'] = self._unpack_int()

    def _read_dim_array(self):
        header = self.fp.read(4)
        assert header in [ZERO, NC_DIMENSION]
        count = self._unpack_int()

        for dim in range(count):
            name = self._unpack_string()
            length = self._unpack_int() or None  # None for record dimension
            self.dimensions[name] = length
            self._dims.append(name)  # preserve order

    def _read_gatt_array(self):
        for k, v in self._read_att_array().items():
            self.__setattr__(k, v)

    def _read_att_array(self):
        header = self.fp.read(4)
        assert header in [ZERO, NC_ATTRIBUTE]
        count = self._unpack_int()

        attributes = {}
        for attr in range(count):
            name = self._unpack_string()
            attributes[name] = self._read_values()
        return attributes

    def _read_var_array(self):
        header = self.fp.read(4)
        assert header in [ZERO, NC_VARIABLE]

        begin = 0
        dtypes = {'names': [], 'formats': []}
        rec_vars = []
        count = self._unpack_int()
        for var in range(count):
            name, dimensions, shape, attributes, typecode, size, dtype_, begin_, vsize = self._read_var()
            if shape and shape[0] is None:
                rec_vars.append(name)
                self.__dict__['_recsize'] += vsize
                if begin == 0: begin = begin_
                dtypes['names'].append(name)
                dtypes['formats'].append(str(shape[1:]) + dtype_)

                # Handle padding with a virtual variable.
                if typecode in 'bch':
                    actual_size = reduce(mul, (1,) + shape[1:]) * size
                    padding = -actual_size % 4
                    if padding:
                        dtypes['names'].append('_padding_%d' % var)
                        dtypes['formats'].append('(%d,)>b' % padding)

                # Data will be set later.
                data = None
            else:
                if self.use_mmap:
                    mm = mmap(self.fp.fileno(), begin_+vsize, access=ACCESS_READ)
                    data = ndarray.__new__(ndarray, shape, dtype=dtype_,
                            buffer=mm, offset=begin_, order=0)
                else:
                    pos = self.fp.tell()
                    self.fp.seek(begin_)
                    data = fromstring(self.fp.read(vsize), dtype=dtype_)
                    data.shape = shape
                    self.fp.seek(pos)

            # Add variable.
            self.variables[name] = netcdf_variable(
                    data, typecode, shape, dimensions, attributes,
                    maskandscale=self._maskandscale)

        if rec_vars:
            # Remove padding when only one record variable.
            if len(rec_vars) == 1:
                dtypes['names'] = dtypes['names'][:1]
                dtypes['formats'] = dtypes['formats'][:1]

            # Build rec array.
            if self.use_mmap:
                mm = mmap(self.fp.fileno(), begin+self._recs*self._recsize, access=ACCESS_READ)
                rec_array = ndarray.__new__(ndarray, (self._recs,), dtype=dtypes,
                        buffer=mm, offset=begin, order=0)
            else:
                pos = self.fp.tell()
                self.fp.seek(begin)
                rec_array = fromstring(self.fp.read(self._recs*self._recsize), dtype=dtypes)
                rec_array.shape = (self._recs,)
                self.fp.seek(pos)

            for var in rec_vars:
                self.variables[var].__dict__['data'] = rec_array[var]

    def _read_var(self):
        name = self._unpack_string()
        dimensions = []
        shape = []
        dims = self._unpack_int()
        
        for i in range(dims):
            dimid = self._unpack_int()
            dimname = self._dims[dimid]
            dimensions.append(dimname)
            dim = self.dimensions[dimname]
            shape.append(dim)
        dimensions = tuple(dimensions)
        shape = tuple(shape)

        attributes = self._read_att_array()
        nc_type = self.fp.read(4)
        vsize = self._unpack_int()
        begin = [self._unpack_int, self._unpack_int64][self.version_byte-1]()

        typecode, size = TYPEMAP[nc_type]
        if typecode is 'c':
            dtype_ = '>c'
        else:
            dtype_ = '>%s' % typecode
            if size > 1: dtype_ += str(size)

        return name, dimensions, shape, attributes, typecode, size, dtype_, begin, vsize

    def _read_values(self):
        nc_type = self.fp.read(4)
        n = self._unpack_int()

        typecode, size = TYPEMAP[nc_type]

        count = n*size
        values = self.fp.read(count)
        self.fp.read(-count % 4)  # read padding

        if typecode is not 'c':
            values = fromstring(values, dtype='>%s%d' % (typecode, size))
            if values.shape == (1,): values = values[0]
        else:
            values = values.rstrip('\x00') 
        return values

    def _pack_begin(self, begin):
        if self.version_byte == 1:
            self._pack_int(begin)
        elif self.version_byte == 2:
            self._pack_int64(begin)

    def _pack_int(self, value):
        self.fp.write(array(value, '>i').tostring())
    _pack_int32 = _pack_int

    def _unpack_int(self):
        return fromstring(self.fp.read(4), '>i')[0]
    _unpack_int32 = _unpack_int

    def _pack_int64(self, value):
        self.fp.write(array(value, '>q').tostring())

    def _unpack_int64(self):
        return fromstring(self.fp.read(8), '>q')[0]

    def _pack_string(self, s):
        count = len(s)
        self._pack_int(count)
        self.fp.write(s)
        self.fp.write('0' * (-count % 4))  # pad

    def _unpack_string(self):
        count = self._unpack_int()
        s = self.fp.read(count).rstrip('\x00')
        self.fp.read(-count % 4)  # read padding
        return s


class netcdf_variable(object):
    """
    ``netcdf_variable`` objects are constructed by calling the method
    ``createVariable`` on the netcdf_file object.

    ``netcdf_variable`` objects behave much like array objects defined in
    Numpy, except that their data resides in a file. Data is read by
    indexing and written by assigning to an indexed subset; the entire
    array can be accessed by the index ``[:]`` or using the methods
    ``getValue`` and ``assignValue``. ``netcdf_variable`` objects also
    have attribute ``shape`` with the same meaning as for arrays, but
    the shape cannot be modified. There is another read-only attribute
    ``dimensions``, whose value is the tuple of dimension names.

    All other attributes correspond to variable attributes defined in
    the NetCDF file. Variable attributes are created by assigning to an
    attribute of the ``netcdf_variable`` object.

    """
    def __init__(self, data, typecode, shape, dimensions, attributes=None,\
                 maskandscale=False):
        self.data = data
        self._typecode = typecode
        self._shape = shape
        self.dimensions = dimensions
        self._maskandscale = maskandscale

        self._attributes = attributes or {}
        for k, v in self._attributes.items():
            self.__dict__[k] = v

    def __setattr__(self, attr, value):
        # Store user defined attributes in a separate dict,
        # so we can save them to file later.
        try:
            self._attributes[attr] = value
        except AttributeError:
            pass
        self.__dict__[attr] = value

    @property
    def isrec(self):
        return self.data.shape and not self._shape[0]

    @property
    def shape(self):
        return self.data.shape

    def getValue(self):
        return self.data.item()

    def assignValue(self, value):
        self.data.itemset(value)

    def typecode(self):
        return self._typecode

    def __getitem__(self, index):
        data = squeeze(self.data[index])
        if self._maskandscale:
            return _unmaskandscale(self,data)
        else:
            return data

    def __setitem__(self, index, data):
        if self._maskandscale:
            data = _maskandscale(self,data)
        # Expand data for record vars?
        if self.isrec:
            if isinstance(index, tuple):
                rec_index = index[0]
            else:
                rec_index = index
            if isinstance(rec_index, slice):
                recs = (rec_index.start or 0) + len(data)
            else:
                recs = rec_index + 1
            if recs > len(self.data):
                shape = (recs,) + self._shape[1:]
                self.data.resize(shape)
        self.data[index] = data


NetCDFFile = netcdf_file
NetCDFVariable = netcdf_variable

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
def _unmaskandscale(var,data):
    # if _maskandscale mode set to True, perform
    # automatic unpacking using scale_factor/add_offset
    # and automatic conversion to masked array using
    # missing_value/_Fill_Value.
    totalmask = zeros(data.shape, bool)
    fill_value = None
    if hasattr(var, 'missing_value') and (data == var.missing_value).any():
        mask=data==var.missing_value
        fill_value = var.missing_value
        totalmask += mask
    if hasattr(var, '_FillValue') and (data == var._FillValue).any():
        mask=data==var._FillValue
        if fill_value is None:
            fill_value = var._FillValue
        totalmask += mask
    else:
        fillval = _default_fillvals[var.typecode()]
        if (data == fillval).any():
            mask=data==fillval
            if fill_value is None:
                fill_value = fillval
            totalmask += mask
    # all values where data == missing_value or _FillValue are
    # masked.  fill_value set to missing_value if it exists,
    # otherwise _FillValue.
    if fill_value is not None:
        data = ma.masked_array(data,mask=totalmask,fill_value=fill_value)
    # if variable has scale_factor and add_offset attributes, rescale.
    if hasattr(var, 'scale_factor') and hasattr(var, 'add_offset'):
        data = var.scale_factor*data + var.add_offset
    return data

def _maskandscale(var,data):
    # if _maskandscale mode set to True, perform
    # automatic packing using scale_factor/add_offset
    # and automatic filling of masked arrays using
    # missing_value/_Fill_Value.
    # use missing_value as fill value.
    # if no missing value set, use _FillValue.
    if hasattr(data,'mask'):
        if hasattr(var, 'missing_value'):
            fillval = var.missing_value
        elif hasattr(var, '_FillValue'):
            fillval = var._FillValue
        else:
            fillval = _default_fillvals[var.typecode()]
        data = data.filled(fill_value=fillval)
    # pack using scale_factor and add_offset.
    if hasattr(var, 'scale_factor') and hasattr(var, 'add_offset'):
        data = (data - var.add_offset)/var.scale_factor
    return data
