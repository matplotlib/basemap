"""DAP variables.

This module is a Python implementation of the DAP data model.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import copy
import itertools

from dap.lib import quote, to_list, _quote
from dap.util.ordereddict import odict
from dap.util.filter import get_filters

__all__ = ['StructureType', 'SequenceType', 'DatasetType', 'GridType', 'ArrayType', 'BaseType',
           'Float', 'Float0', 'Float8', 'Float16', 'Float32', 'Float64', 'Int', 'Int0', 'Int8',
           'Int16', 'Int32', 'Int64', 'UInt16', 'UInt32', 'UInt64', 'Byte', 'String', 'Url']

_basetypes = ['Float32', 'Float64', 'Int16', 'Int32', 'UInt16', 'UInt32', 'Byte', 'String', 'Url']
_constructors = ['StructureType', 'SequenceType', 'DatasetType', 'GridType', 'ArrayType']

# Constants.
Float   = 'Float64'
Float0  = 'Float64'
Float8  = 'Float32'
Float16 = 'Float32'
Float32 = 'Float32'
Float64 = 'Float64'
Int     = 'Int32'
Int0    = 'Int32'
Int8    = 'Byte'
Int16   = 'Int16'
Int32   = 'Int32'
Int64   = 'Int32'
UInt16  = 'UInt16'
UInt32  = 'UInt32'
UInt64  = 'UInt32'
UInt8   = 'Byte'
Byte    = 'Byte'
String  = 'String'
Url     = 'Url'

typemap = {
           # numpy
           'd': Float64,
           'f': Float32,
           'l': Int32,
           'b': Byte,
           'h': Int16,
           'q': Int32,
           'H': UInt16,
           'L': UInt32,
           'Q': UInt32,
           'B': Byte,
           'S': String,
          }


class StructureType(odict):
    """Structure contructor.  

    A structure is a dict-like object, which can hold other DAP variables.
    Structures have a 'data' attribute that combines the data from the
    stored variables when read, and propagates the data to the variables
    when set.
    
    This behaviour can be bypassed by setting the '_data' attribute; in
    this case, no data is propagated, and further reads do not combine the
    data from the stored variables.
    """

    def __init__(self, name='', attributes=None):
        odict.__init__(self)
        self.name = quote(name)
        self.attributes = attributes or {}

        self._id = name
        self._filters = []
        self._data = None

    def __iter__(self):
        # Iterate over the variables contained in the structure.
        return self.itervalues()
    
    walk = __iter__

    def __getattr__(self, attr):
        # Try to return stored variable.
        try:
            return self[attr]
        except KeyError:
            # Try to return attribute from self.attributes.
            try: return self.attributes[attr]
            except KeyError: raise AttributeError

    def __setitem__(self, key, item):
        # Assign a new variable and apply the proper id.
        self._dict.__setitem__(key, item)
        if key not in self._keys: self._keys.append(key)

        # Ensure that stored objects have the proper id.
        item._set_id(self._id)

    def _get_data(self):
        if self._data is not None:
            return self._data
        else:
            return [var.data for var in self.values()]
        
    def _set_data(self, data):
        # Propagate the data to the stored variables.
        for data_, var in itertools.izip(data, self.values()):
            var.data = data_
            
    data = property(_get_data, _set_data)

    def _get_id(self):
        return self._id
        
    def _set_id(self, parent=None):
        if parent:
            self._id = '%s.%s' % (parent, self.name)
        else:
            self._id = self.name

        # Propagate id to stored variables.
        for var in self.values(): var._set_id(self._id)
        
    id = property(_get_id)  # Read-only.
    
    def _get_filters(self):
        return self._filters
        
    def _set_filters(self, f):
        self._filters.append(f)

        # Propagate filter to stored variables.
        for var in self.values(): var._set_filters(f)
        
    filters = property(_get_filters, _set_filters)

    def __copy__(self):
        out = self.__class__(name=self.name, attributes=self.attributes.copy())
        out._id = self._id
        out._filters = self._filters[:]
        out._data = self._data

        # Stored variables *are not* copied.
        for k, v in self.items():
            out[k] = v
        return out

    def __deepcopy__(self, memo=None, _nil=[]):
        out = self.__class__(name=self.name, attributes=self.attributes.copy())
        out._id = self._id
        out._filters = self._filters[:]
        out._data = self._data

        # Stored variables *are* (deep) copied.
        for k, v in self.items():
            out[k] = copy.deepcopy(v)
        return out


class DatasetType(StructureType):
    """Dataset constructor.

    A dataset is very similar to a structure -- the main difference is that
    its name is not used when composing the fully qualified name of stored
    variables.
    """
    
    def __setitem__(self, key, item):
        self._dict.__setitem__(key, item)
        if key not in self._keys: self._keys.append(key)

        # Set the id. Here the parent should be None, since the dataset
        # id is not part of the fully qualified name.
        item._set_id(None)

    def _set_id(self, parent=None):
        self._id = self.name

        # Propagate id.
        for var in self.values(): var._set_id(None)


class SequenceType(StructureType):
    """Sequence constructor.

    A sequence contains ordered data, corresponding to the records in
    a sequence of structures with the same stored variables.
    """
    # Nesting level. Sequences inside sequences have a level 2, and so on.
    level = 1

    def __setitem__(self, key, item):
        # Assign a new variable and apply the proper id.
        self._dict.__setitem__(key, item)
        if key not in self._keys: self._keys.append(key)

        # Ensure that stored objects have the proper id.
        item._set_id(self._id)

        # If the variable is a sequence, set the nesting level.
        def set_level(seq, level):
            if isinstance(seq, SequenceType):
                seq.level = level
                for child in seq.walk(): set_level(child, level+1)
        set_level(item, self.level+1)

    def walk(self):
        # Walk over the variables contained in the structure.
        return self.itervalues()
    
    def _get_data(self):
        # This is similar to the structure _get_data method, except that data
        # is combined from stored variables using zip(), i.e., grouped values 
        # from each variable.
        if self._data is not None:
            return self._data
        else: 
            return _build_data(self.level, *[var.data for var in self.values()])
        
    def _set_data(self, data):
        for data_, var in itertools.izip(_propagate_data(self.level, data), self.values()):
            var.data = data_

    data = property(_get_data, _set_data)

    def __iter__(self):
        """
        When iterating over a sequence, we yield structures containing the
        corresponding data (first record, second, etc.).
        """
        out = self.__deepcopy__()

        # Set server-side filters. When the sequence is iterated in a
        # listcomp/genexp, this function inspects the stack and tries to
        # build a server-side filter from the client-side filter. This
        # is voodoo black magic, take care.
        filters = get_filters(out)
        for filter_ in filters:
            out._set_filters(filter_)

        for values in out.data:
            # Yield a nice structure.
            struct_ = StructureType(name=out.name, attributes=out.attributes)
            for data, name in zip(values, out.keys()):
                var = struct_[name] = out[name].__deepcopy__()
                var.data = data
            # Set the id. This is necessary since the new structure is not
            # contained inside a dataset.
            parent = out._id[:-len(out.name)-1]
            struct_._set_id(parent)
            yield struct_

    def filter(self, *filters):
        # Make a copy of the sequence.
        out = self.__deepcopy__()

        # And filter it according to the selection expression.
        for filter_ in filters:
            out._set_filters(_quote(filter_))
        return out


class GridType(object):
    """Grid constructor.

    A grid is a constructor holding an 'array' variable. The array has its
    dimensions mapped to 'maps' stored in the grid (lat, lon, time, etc.).
    Most of the requests are simply passed onto the stored array.
    """

    def __init__(self, name='', array=None, maps=None, attributes=None):
        self.name = quote(name)
        self.array = array
        self.maps = maps or odict()
        self.attributes = attributes or {}

        self._id = name
        self._filters = []

    def __len__(self):
        return self.array.shape[0]

    def __iter__(self):
        # Iterate over the grid. Yield the array and then the maps.
        yield self.array
        for map_ in self.maps.values(): yield map_

    walk = __iter__
            
    def __getattr__(self, attr):
        # Try to return attribute from self.attributes.
        try:
            return self.attributes[attr]
        except KeyError:
            raise AttributeError

    def __getitem__(self, index):
        # Return data from the array.
        return self.array[index]

    def _get_data(self):
        return self.array.data

    def _set_data(self, data):
        self.array.data = data

    data = property(_get_data, _set_data)

    def _get_id(self):
        return self._id

    def _set_id(self, parent=None):
        if parent: self._id = '%s.%s' % (parent, self.name)
        else: self._id = self.name

        # Propagate id to array and maps.
        if self.array: self.array._set_id(self._id)
        for map_ in self.maps.values():
            map_._set_id(self._id)

    id = property(_get_id)

    def _get_filters(self):
        return self._filters

    def _set_filters(self, f):
        self.filters.append(f)

        # Propagate filter.
        self.array._set_filters(f)
        for map_ in self.maps.values():
            map_._set_filters(f)

    filters = property(_get_filters, _set_filters)

    def _get_dimensions(self):
        # Return dimensions from stored maps.
        return tuple(self.maps.keys())

    dimensions = property(_get_dimensions)

    def _get_shape(self):
        return self.array.shape

    def _set_shape(self, shape):
        self.array.shape = shape

    shape = property(_get_shape, _set_shape)

    def _get_type(self):
        return self.array.type

    def _set_type(self, type):
        self.array.type = type

    type = property(_get_type, _set_type)

    def __copy__(self):
        out = self.__class__(name=self.name, array=self.array, maps=self.maps, attributes=self.attributes.copy())
        out._id = self._id
        out._filters = self._filters[:]
        return out

    def __deepcopy__(self, memo=None, _nil=[]):
        out = self.__class__(name=self.name, attributes=self.attributes.copy())
        out.array = copy.deepcopy(self.array)
        out.maps = copy.deepcopy(self.maps)
        out._id = self._id
        out._filters = self._filters[:]
        return out


class BaseType(object):
    """DAP Base type.

    Variable holding a single value, or an iterable if it's stored inside
    a sequence. It's the fundamental DAP variable, which actually holds
    data (together with arrays).
    """

    def __init__(self, name='', data=None, type=None, attributes=None):
        self.name = quote(name)
        self.data = data
        self.attributes = attributes or {}
        
        if type in _basetypes: self.type = type
        else: self.type = typemap.get(type, Int32)

        self._id = name
        self._filters = []

    def __iter__(self):
        # Yield the stored value.
        # Perhaps we should raise StopIteration?
        yield self.data

    def __getattr__(self, attr):
        # Try to return attribute from self.attributes.
        try:
            return self.attributes[attr]
        except KeyError:
            raise AttributeError

    def __getitem__(self, key):
        # Return data from the array.
        return self.data[key]

    def __setitem__(self, key, item):
        # Assign a new variable and apply the proper id.
        self.data.__setitem__(key, item)

    def _get_id(self):
        return self._id

    def _set_id(self, parent=None):
        if parent: self._id = '%s.%s' % (parent, self.name)
        else: self._id = self.name

    id = property(_get_id)  # Read-only.

    def _get_filters(self):
        return self._filters

    def _set_filters(self, f):
        self._filters.append(f)

        # Propagate to data, if it's a Proxy object.
        if hasattr(self.data, 'filters'):
            self.data.filters = self._filters

    filters = property(_get_filters, _set_filters)

    def __copy__(self):
        out = self.__class__(name=self.name, data=self.data, type=self.type, attributes=self.attributes.copy())
        out._id = self._id
        out._filters = self._filters[:]
        return out

    def __deepcopy__(self, memo=None, _nil=[]):
        out = self.__class__(name=self.name, type=self.type, attributes=self.attributes.copy())

        try:
            out.data = copy.copy(self.data)
        except TypeError:
            self.data = to_list(self.data)
            out.data = copy.copy(self.data)
        
        out._id = self._id
        out._filters = self._filters[:]
        return out

    # This allows the variable to be compared to numbers.
    def __ge__(self, other): return self.data >= other
    def __gt__(self, other): return self.data >  other
    def __le__(self, other): return self.data <= other
    def __lt__(self, other): return self.data <  other
    def __eq__(self, other): return self.data == other


class ArrayType(BaseType):
    """An array of BaseType variables.
    
    Although the DAP supports arrays of any DAP variables, pydap can only
    handle arrays of base types. This makes the ArrayType class very
    similar to a BaseType, with the difference that it'll hold an array
    of data in its 'data' attribute.

    Array of constructors will not be supported until Python has a
    native multi-dimensional array type. 
    """

    def __init__(self, name='', data=None, shape=None, dimensions=None, type=None, attributes=None):
        self.name = quote(name)
        self.data = data
        self.shape = shape or ()
        self.dimensions = dimensions or ()
        self.attributes = attributes or {}

        if type in _basetypes: self.type = type
        else: self.type = typemap.get(type, Int32)
        
        self._id = name
        self._filters = []

    def __len__(self):
        return self.shape[0]

    def __copy__(self):
        out = self.__class__(name=self.name, data=self.data, shape=self.shape, dimensions=self.dimensions, type=self.type, attributes=self.attributes.copy())
        out._id = self._id
        out._filters = self._filters[:]
        return out

    def __deepcopy__(self, memo=None, _nil=[]):
        out = self.__class__(name=self.name, shape=self.shape, dimensions=self.dimensions, type=self.type, attributes=self.attributes.copy())
        
        try:
            out.data = copy.copy(self.data)
        except TypeError:
            self.data = to_list(self.data)
            out.data = copy.copy(self.data)
        
        out._id = self._id
        out._filters = self._filters[:]
        return out


# Functions for propagating data up and down in sequences.
# I'm not 100% sure how this works.
def _build_data(level, *vars_):
    if level > 0:
        out = [_build_data(level-1, *els) for els in itertools.izip(*vars_)]
    else:
        out = vars_

    return out

def _propagate_data(level, vars_):
    if level > 0:
        out = zip(*[_propagate_data(level-1, els) for els in vars_])
    else:
        out = vars_

    return out
