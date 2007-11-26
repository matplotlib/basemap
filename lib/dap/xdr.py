# This Python file uses the following encoding: utf-8

"""Fast(er) implementation of xdrlib for the DAP.

This module reimplements Python's xdrlib module for the DAP. It uses the
"array" module for speed, and encodes bytes according to the DAP spec.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import array
import struct
import operator

try:
    from numpy import array as array_
except ImportError:
    array_ = None

from dap.lib import isiterable


_big_endian = struct.pack('i',1)[0] != '\x01'


# Convert from DAP types to array types.
typeconvert = {'float64': ('d', 8),
               'float32': ('f', 4),
               'float'  : ('f', 4),
               'uint'   : ('I', 4),
               'uint16' : ('I', 4),
               'uint32' : ('I', 4),
               'int'    : ('i', 4),
               'int16'  : ('i', 4),
               'int32'  : ('i', 4),
               'byte'   : ('b', 1), 
              }

                     
class DapPacker(object):
    r"""Pack variable data into XDR.

    This is a faster reimplementation of xdrlib using the native module
    "array".

        >>> from dap.dtypes import *
        >>> dapvar = BaseType(data=1, type='Int32')
        >>> xdrdata = DapPacker(dapvar)
        >>> for i in xdrdata:
        ...     print repr(i)
        '\x00\x00\x00\x01'

        >>> dapvar = ArrayType(data=['one', 'two'], shape=[2], type='String')
        >>> xdrdata = DapPacker(dapvar)
        >>> for i in xdrdata:
        ...     print repr(i)
        '\x00\x00\x00\x02'
        '\x00\x00\x00\x03one\x00'
        '\x00\x00\x00\x03two\x00'

        >>> dapvar = ArrayType(data=range(2), shape=[2], type='Int32')
        >>> xdrdata = DapPacker(dapvar)
        >>> for i in xdrdata:
        ...     print repr(i)
        '\x00\x00\x00\x02\x00\x00\x00\x02'
        '\x00\x00\x00\x00\x00\x00\x00\x01'

        >>> dapvar = ArrayType(data=range(2), shape=[2], type='Float64')
        >>> xdrdata = DapPacker(dapvar)
        >>> for i in xdrdata:
        ...     print repr(i)
        '\x00\x00\x00\x02\x00\x00\x00\x02'
        '\x00\x00\x00\x00\x00\x00\x00\x00?\xf0\x00\x00\x00\x00\x00\x00'
    """
    def __init__(self, dapvar, data=None):
        self.var = dapvar

        if data is None:
            data = dapvar.data

        # Make data iterable.
        if not isiterable(data): data = [data]

        # Put data in blocks.
        if not hasattr(data, 'shape') or len(data.shape) <= 1:
            data = [data]

        self.data = data

    def __iter__(self):
        """Iterate over the XDR encoded data."""
        # Yield length (twice) if array.
        if hasattr(self.var, 'shape'):
            if self.var.type.lower() in ['url', 'string']:
                yield self._pack_length()
            else:
                yield self._pack_length() * 2

        # Bytes are sent differently.
        if self.var.type.lower() == 'byte':
            for b in self._yield_bytes(): yield b
        # String are zero padded to 4n.
        elif self.var.type.lower() in ['url', 'string']:
            for block in self.data:
                if hasattr(block, 'flat'): block = block.flat  # for numpy str arrays
                for word in block:
                    yield self._pack_string(word)
        else:
            type_, size = typeconvert[self.var.type.lower()]
            for block in self.data:
                # Flatten block.
                if hasattr(block, 'flat'): block = block.flat

                # Encode block in XDR (big-endian). First we
                # try to encode using numpy's ``.astype``,
                # otherwise be use ``array.array`` and ensure
                # data is in big endian format.
                try:
                    dtype = ">%s%s" % (type_, size)
                    data = block.base.astype(dtype).data
                    data = str(data)
                except:
                    data = array.array(type_, list(block)) 
                    if not _big_endian: data.byteswap()
                    data = data.tostring()

                yield data

    def _pack_length(self):
        """Yield array length."""
        shape = getattr(self.var, 'shape', [1])
        length = reduce(operator.mul, shape)
        return struct.pack('>L', length)

    def _yield_bytes(self):
        r"""Yield bytes.

        Bytes are encoded as is, padded to a four-byte boundary. An array
        of five bytes, eg, is encoded as eight bytes:

            >>> from dap.dtypes import *
            >>> dapvar = ArrayType(data=range(5), shape=[5], type='Byte')
            >>> xdrdata = DapPacker(dapvar)
            >>> for i in xdrdata:
            ...     print repr(i)
            '\x00\x00\x00\x05\x00\x00\x00\x05'
            '\x00'
            '\x01'
            '\x02'
            '\x03'
            '\x04'
            '\x00\x00\x00'

        Again, the first line correponds to the array size packed twice,
        followed by the bytes and the padding.
        """
        count = 0
        for block in self.data:
            data = array.array('B', block).tostring()
            for d in data:
                yield d
                count += 1

        padding = (4 - (count % 4)) % 4
        yield padding * '\0'

    def _pack_string(self, s):
        """Pack a string.
        
        We first pack the string length, followed by the string padded
        to size 4n.
        """
        # Pack length first.
        n = len(s)
        length = struct.pack('>L', n)

        n = ((n+3)/4)*4
        data = length + s + (n - len(s)) * '\0'
        return data


class DapUnpacker(object):
    r"""A XDR data unpacker.

    Unpacking data from a base type:

        >>> from dap.dtypes import *
        >>> dapvar = BaseType(data=1, name='dapvar', type='Byte')
        >>> print dapvar.data
        1
        >>> from dap.server import SimpleHandler
        >>> from dap.helper import escape_dods
        >>> dataset = DatasetType(name='test')
        >>> dataset['dapvar'] = dapvar
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> print escape_dods(''.join(output), pad='')
        Dataset {
            Byte dapvar;
        } test;
        Data:
        \x01\x00\x00\x00
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> xdrdata = ''.join(output)
        >>> start = xdrdata.index('Data:\n') + len('Data:\n')
        >>> xdrdata = xdrdata[start:]
        >>> data = DapUnpacker(xdrdata, (), dapvar.type)
        >>> print data.getvalue()
        1
    
    An array of bytes:

        >>> dapvar = ArrayType(name='dapvar', data=range(5), shape=[5], type='Byte')
        >>> dataset['dapvar'] = dapvar
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> xdrdata = ''.join(output)
        >>> start = xdrdata.index('Data:\n') + len('Data:\n')
        >>> xdrdata = xdrdata[start:]
        >>> data = DapUnpacker(xdrdata, dapvar.shape, dapvar.type, outshape=[5])
        >>> print data.getvalue()
        [0 1 2 3 4]

    Another array:

        >>> dapvar = ArrayType(name='dapvar', data=range(25), shape=[5,5], type='Float32')
        >>> dataset['dapvar'] = dapvar
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> xdrdata = ''.join(output)
        >>> start = xdrdata.index('Data:\n') + len('Data:\n')
        >>> xdrdata = xdrdata[start:]
        >>> data = DapUnpacker(xdrdata, dapvar.shape, dapvar.type, outshape=[5,5])
        >>> print data.getvalue()
        [[  0.   1.   2.   3.   4.]
         [  5.   6.   7.   8.   9.]
         [ 10.  11.  12.  13.  14.]
         [ 15.  16.  17.  18.  19.]
         [ 20.  21.  22.  23.  24.]]

    One more:
        
        >>> dapvar = ArrayType(name='dapvar', data=['um', 'dois', 'três'], shape=[3], type='String')
        >>> dataset['dapvar'] = dapvar
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> for line in output:
        ...     print repr(line)
        'Dataset {\n'
        '    String dapvar[dapvar = 3];\n'
        '} test;\n'
        'Data:\n'
        '\x00\x00\x00\x03'
        '\x00\x00\x00\x02um\x00\x00'
        '\x00\x00\x00\x04dois'
        '\x00\x00\x00\x05tr\xc3\xaas\x00\x00\x00'
        >>> headers, output = SimpleHandler(dataset).dods()
        >>> xdrdata = ''.join(output)
        >>> start = xdrdata.index('Data:\n') + len('Data:\n')
        >>> xdrdata = xdrdata[start:]
        >>> data = DapUnpacker(xdrdata, dapvar.shape, dapvar.type, outshape=[3])
        >>> print data.getvalue()
        [um dois três]
    """
    def __init__(self, data, shape, type, outshape=None):
        self.__buf = data
        self.shape = shape
        self.type = type
        self.outshape = outshape
        
        # Buffer position.
        self.__pos = 0

    def getvalue(self):
        """Unpack the XDR-encoded data."""
        # Get current position.
        i = self.__pos

        # Check for empty sequence.
        if self.__buf[i:i+4] == '\xa5\x00\x00\x00':
            return []

        # Check for sequence with data.
        elif self.__buf[i:i+4] == '\x5a\x00\x00\x00':
            # Unpack sequence start marker (uint 1509949440).
            mark = self._unpack_uint()
            out = []
            while mark != 2768240640L:
                tmp = self.getvalue()
                out.append(tmp)
                # Unpack marker.
                mark = self._unpack_uint()
            return out

        # Get data length.
        n = 1
        if self.shape:
            n = self._unpack_uint()
            # Strings pack the size only once?
            if self.type.lower() not in ['url', 'string']:
                self._unpack_uint()

        # Bytes are treated differently.
        if self.type.lower() == 'byte':
            out = self._unpack_bytes(n)
            type_, size = typeconvert[self.type.lower()]
        # As are strings...
        elif self.type.lower() in ['url', 'string']:
            out = self._unpack_string(n)
            type_ = 'S'
        else:
            i = self.__pos
            type_, size = typeconvert[self.type.lower()]
            out = array.array(type_, self.__buf[i:i+n*size])
            self.__pos = i+n*size

            # Ensure big-endianess.
            if not _big_endian: out.byteswap()

        if not self.shape:
            out = out[0]  # BaseType.
        elif array_:
            # Convert to array and reshape.
            out = array_(out, type_)
            if self.outshape: out.shape = tuple(self.outshape)

        return out

    def _unpack_uint(self):
        i = self.__pos
        self.__pos = j = i+4
        data = self.__buf[i:j]
        if len(data) < 4:
            raise EOFError
        x = struct.unpack('>L', data)[0]
        try:
            return int(x)
        except OverflowError:
            return x

    def _unpack_bytes(self, count):
        i = self.__pos
        out = array.array('b', self.__buf[i:i+count])
        padding = (4 - (count % 4)) % 4
        self.__pos = i + count + padding
        
        return out
        
    def _unpack_string(self, count):
        out = []
        for s in range(count):
            # Unpack string length.
            n = self._unpack_uint()

            i = self.__pos
            j = i+n
            data = self.__buf[i:j]

            # Fix cursor position.
            padding = (4 - (n % 4)) % 4
            self.__pos = j + padding
            
            out.append(data)
        return out


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

