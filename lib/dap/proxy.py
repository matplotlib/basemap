from __future__ import division

"""Proxy class to DAP data.

This module implements a proxy object that behaves like an array and 
downloads data transparently from a DAP server when sliced. It is used
when building a representation of the dataset, and should be not 
directly used.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import sys

from dap.xdr import DapUnpacker
from dap.helper import fix_slice
from dap.util.http import openurl

try:
    from numpy import ndarray as ndarray_
except ImportError:
    ndarray_ = None


class Proxy(object):
    """
    A proxy to data stored in a DODS server.
    """
    def __init__(self, url, id, shape, type, filters=None, cache=None, username=None, password=None):
        self.url = url
        self.id = id
        self.shape = shape
        self.type = type

        self.filters = filters or []

        self.cache = cache
        self.username = username
        self.password = password

    def __iter__(self):
        return iter(self[:])

    def __getitem__(self, index):
        """
        Download data from DAP server.
        
        When the proxy object is sliced, it build an URL from the slice
        and retrieves the data from the DAP server.
        """
        # Build the base URL.
        url = '%s.dods?%s' % (self.url, self.id)

        if self.shape:
            # Fix the index for incomplete slices or ellipsis.
            index = fix_slice(len(self.shape), index)

            # Force index to tuple, to iterate over the slices.
            if not isinstance(index, tuple):
                index = index,

            # Iterate over the sliced dimensions.
            i = 0
            outshape = []
            for dimension in index:
                # If dimension is a slice, get start, step and stop.
                if isinstance(dimension, slice):
                    start = dimension.start  or 0
                    step  = dimension.step   or 1
                    if dimension.stop:
                        stop = dimension.stop-1
                    else:
                        stop = self.shape[i]-1

                # Otherwise, retrieve a single value.
                else:
                    start = dimension
                    stop  = dimension
                    step  = 1

                # When stop is not specified, use the shape.
                if stop == sys.maxint or stop > self.shape[i]-1:
                    stop = self.shape[i]-1
                # Negative slices. 
                elif stop < 0:
                    stop = self.shape[i]+stop

                # Negative starting slices.
                if start < 0: start = self.shape[i]+start

                # Build the URL used to retrieve the data.
                url = '%s[%s:%s:%s]' % (url, str(start), str(step), str(stop))

                # outshape is a list of the slice dimensions.
                outshape.append(1+(stop-start)//step)
                
                # Update to next dimension.
                i += 1

        else:
            # No need to resize the data.
            outshape = None

        # Make the outshape consistent with the numpy and pytables conventions.
        if outshape is not None:
            outshape = self._reduce_outshape(outshape)

        # Check for filters.
        if self.filters:
            ce = '&'.join(self.filters)
            url = '%s&%s' % (url, ce)

        # Fetch data.
        resp, data = openurl(url, self.cache, self.username, self.password)

        # First lines are ASCII information that end with 'Data:\n'.
        start = data.index('Data:\n') + len('Data:\n')
        xdrdata = data[start:]

        # Unpack data.
        output = DapUnpacker(xdrdata, self.shape, self.type, outshape).getvalue()

        # Convert length 1 arrays to scalars
        if ndarray_:
            if outshape == () and isinstance(output, ndarray_):
                output = output[0]

        return output


    def _reduce_outshape(self, outshape):
        """Make the outshape consistent with the numpy and pytables conventions.

        (1, N) -> (N,)
        (N, 1) -> (N,)
        (N, M, 1) -> (N, M)

        Basically, all ones are removed from the shape.
        """
        return tuple([index for index in outshape if index != 1])

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
