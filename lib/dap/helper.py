"""Helper functions.

These are generic functions used mostly for writing plugins.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import sys
import re
import operator
import itertools
import copy
from urllib import quote, unquote

from dap.dtypes import *
from dap.dtypes import _basetypes
from dap.exceptions import ConstraintExpressionError
from dap.lib import isiterable
from dap.util.safeeval import expr_eval
from dap.util.ordereddict import odict


def constrain(dataset, constraints):
    """A simple example. We create a dataset holding three variables:

        >>> dataset = DatasetType(name='foo')
        >>> dataset['a'] = BaseType(name='a', type='Byte')
        >>> dataset['b'] = BaseType(name='b', type='Byte')
        >>> dataset['c'] = BaseType(name='c', type='Byte')

    Now we give it a CE requesting only the variables ``a`` and ``b``:
    
        >>> dataset2 = constrain(dataset, 'a,b')
        >>> print dataset2  #doctest: +ELLIPSIS
        {'a': <dap.dtypes.BaseType object at ...>, 'b': <dap.dtypes.BaseType object at ...>}

    We can also request the variables in a different order:

        >>> dataset2 = constrain(dataset, 'b,a')
        >>> print dataset2  #doctest: +ELLIPSIS
        {'b': <dap.dtypes.BaseType object at ...>, 'a': <dap.dtypes.BaseType object at ...>}

    Another example. A dataset with two structures ``a`` and ``b``:

        >>> dataset = DatasetType(name='foo')
        >>> dataset['a'] = StructureType(name='a')
        >>> dataset['a']['a1'] = BaseType(name='a1', type='Byte')
        >>> dataset['b'] = StructureType(name='b')
        >>> dataset['b']['b1'] = BaseType(name='b1', type='Byte')
        >>> dataset['b']['b2'] = BaseType(name='b2', type='Byte')

    If we request the structure ``b`` we should get it complete:

        >>> dataset2 = constrain(dataset, 'a.a1,b')
        >>> print dataset2  #doctest: +ELLIPSIS
        {'a': {'a1': <dap.dtypes.BaseType object at ...>}, 'b': {'b1': <dap.dtypes.BaseType object at ...>, 'b2': <dap.dtypes.BaseType object at ...>}}

        >>> dataset2 = constrain(dataset, 'b.b1')
        >>> print dataset2  #doctest: +ELLIPSIS
        {'b': {'b1': <dap.dtypes.BaseType object at ...>}}

    Arrays can be sliced. Here we have a ``(2,3)`` array:

        >>> dataset = DatasetType(name='foo')
        >>> from numpy import array
        >>> data = array([1,2,3,4,5,6])
        >>> data.shape = (2,3)
        >>> dataset['array'] = ArrayType(data=data, name='array', shape=(2,3), type='Int32')
        >>> dataset2 = constrain(dataset, 'array')
        >>> from dap.server import SimpleHandler
        >>> headers, output = SimpleHandler(dataset).dds()
        >>> print ''.join(output)
        Dataset {
            Int32 array[2][3];
        } foo;
        <BLANKLINE>
        >>> print dataset2['array'].data
        [[1 2 3]
         [4 5 6]]

    But we request only part of it:

        >>> dataset2 = constrain(dataset, 'array[0:1:1][0:1:1]')
        >>> headers, output = SimpleHandler(dataset2).dds()
        >>> print ''.join(output)
        Dataset {
            Int32 array[2][2];
        } foo;
        <BLANKLINE>
        >>> print dataset2['array'].data
        [[1 2]
         [4 5]]

    The same is valid for grids:

        >>> dataset['grid'] = GridType(name='grid') 
        >>> data = array([1,2,3,4,5,6])
        >>> data.shape = (2,3)
        >>> dataset['grid'].array  = ArrayType(name='grid', data=data, shape=(2,3), dimensions=('x', 'y'))
        >>> dataset['grid'].maps['x'] = ArrayType(name='x', data=array([1,2]), shape=(2,))
        >>> dataset['grid'].maps['y'] = ArrayType(name='y', data=array([1,2,3]), shape=(3,))
        >>> dataset._set_id()
        >>> headers, output = SimpleHandler(dataset).dds()
        >>> print ''.join(output)
        Dataset {
            Int32 array[2][3];
            Grid {
                Array:
                    Int32 grid[x = 2][y = 3];
                Maps:
                    Int32 x[x = 2];
                    Int32 y[y = 3];
            } grid;
        } foo;
        <BLANKLINE>
        >>> dataset2 = constrain(dataset, 'grid[0:1:0][0:1:0]')
        >>> headers, output = SimpleHandler(dataset2).dds()
        >>> print ''.join(output)
        Dataset {
            Grid {
                Array:
                    Int32 grid[x = 1][y = 1];
                Maps:
                    Int32 x[x = 1];
                    Int32 y[y = 1];
            } grid;
        } foo;
        <BLANKLINE>
        >>> headers, output = SimpleHandler(dataset2).ascii()
        >>> print ''.join(output)
        Dataset {
            Grid {
                Array:
                    Int32 grid[x = 1][y = 1];
                Maps:
                    Int32 x[x = 1];
                    Int32 y[y = 1];
            } grid;
        } foo;
        ---------------------------------------------
        grid.grid
        [0][0] 1
        <BLANKLINE>
        grid.x
        [0] 1
        <BLANKLINE>
        grid.y
        [0] 1
        <BLANKLINE>
        <BLANKLINE>
        <BLANKLINE>
        <BLANKLINE>

    Selecting a map from a Grid should return a structure:

        >>> dataset3 = constrain(dataset, 'grid.x')
        >>> headers, output = SimpleHandler(dataset3).dds()
        >>> print ''.join(output)
        Dataset {
            Structure {
                Int32 x[x = 2];
            } grid;
        } foo;
        <BLANKLINE>

    Short notation also works:

        >>> dataset3 = constrain(dataset, 'x')
        >>> headers, output = SimpleHandler(dataset3).dds()
        >>> print ''.join(output)
        Dataset {
            Structure {
                Int32 x[x = 2];
            } grid;
        } foo;
        <BLANKLINE>

    It also works with Sequences:

        >>> dataset = DatasetType(name='foo')
        >>> dataset['seq'] = SequenceType(name='seq')
        >>> dataset['seq']['a'] = BaseType(name='a')
        >>> dataset['seq']['b'] = BaseType(name='b')
        >>> dataset['seq']['a'].data = range(5)
        >>> dataset['seq']['b'].data = range(5,10)
        >>> for i in dataset['seq'].data:
        ...     print i
        (0, 5)
        (1, 6)
        (2, 7)
        (3, 8)
        (4, 9)
        >>> dataset2 = constrain(dataset, 'seq.a')
        >>> for i in dataset2['seq'].data:
        ...     print i
        (0,)
        (1,)
        (2,)
        (3,)
        (4,)
        >>> dataset2 = constrain(dataset, 'seq.b')
        >>> for i in dataset2['seq'].data:
        ...     print i
        (5,)
        (6,)
        (7,)
        (8,)
        (9,)
        >>> dataset2 = constrain(dataset, 'seq.b,seq.a')
        >>> for i in dataset2['seq'].data:
        ...     print i
        (5, 0)
        (6, 1)
        (7, 2)
        (8, 3)
        (9, 4)

    The function also parses selection expressions. Let's create a
    dataset with sequential data:

        >>> dataset = DatasetType(name='foo')
        >>> dataset['seq'] = SequenceType(name='seq')
        >>> dataset['seq']['index'] = BaseType(name='index', type='Int32')
        >>> dataset['seq']['index'].data = [10, 11, 12, 13]
        >>> dataset['seq']['temperature'] = BaseType(name='temperature', type='Float32')
        >>> dataset['seq']['temperature'].data = [17.2, 15.1, 15.3, 15.1]
        >>> dataset['seq']['site'] = BaseType(name='site', type='String')
        >>> dataset['seq']['site'].data = ['Diamond_St', 'Blacktail_Loop', 'Platinum_St', 'Kodiak_Trail']

    Here's the data:

        >>> for i in dataset['seq'].data:
        ...     print i
        (10, 17.199999999999999, 'Diamond_St')
        (11, 15.1, 'Blacktail_Loop')
        (12, 15.300000000000001, 'Platinum_St')
        (13, 15.1, 'Kodiak_Trail')

    Now suppose we only want data where ``index`` is greater than 11:

        >>> dataset2 = constrain(dataset, 'seq&seq.index>11')
        >>> for i in dataset2['seq'].data:
        ...     print i
        (12, 15.300000000000001, 'Platinum_St')
        (13, 15.1, 'Kodiak_Trail')

    We can request only a few variables:

        >>> dataset2 = constrain(dataset, 'seq.site&seq.index>11')
        >>> for i in dataset2['seq'].data:
        ...     print i
        ('Platinum_St',)
        ('Kodiak_Trail',)

    A few more tests:

        >>> dataset = DatasetType(name='foo')
        >>> dataset['a'] = StructureType(name='a')
        >>> dataset['a']['shn'] = BaseType(name='shn')
        >>> dataset['b'] = StructureType(name='b')
        >>> dataset['b']['shn'] = BaseType(name='shn')
        >>> dataset2 = constrain(dataset, 'a.shn')
        >>> print dataset2  #doctest: +ELLIPSIS
        {'a': {'shn': <dap.dtypes.BaseType object at ...>}}
        >>> dataset3 = constrain(dataset, 'shn')
        Traceback (most recent call last):
        ...
        ConstraintExpressionError: 'Ambiguous shorthand notation request: shn'
        >>> dataset['shn'] = BaseType(name='shn')
        >>> dataset3 = constrain(dataset, 'shn')
        >>> print dataset3  #doctest: +ELLIPSIS
        {'shn': <dap.dtypes.BaseType object at 0x1746290>}
    """
    # Parse constraints.
    fields, queries = parse_querystring(constraints)

    # Ids and names are used to check that requests made using the 
    # shorthand notation are not ambiguous. Used names are stored to
    # make sure that at most only a single variables is returned from
    # a given name.
    ids = [var.id for var in walk(dataset)]
    names = []
    
    new = DatasetType(name=dataset.name, attributes=dataset.attributes.copy())
    new = build(dataset, new, fields, queries, ids, names)

    return new


def build(dapvar, new, fields, queries, ids, names):
    vars_ = fields.keys()
    order = []
    for var in dapvar.walk():
        # Make a copy of the variable, so that later we can possibly add it 
        # to the dataset we're building (that's why it's a candidate).
        candidate = copy.deepcopy(var)

        # We first filter the data in sequences. This has to be done 
        # before variables are removed, since we can select values based
        # on conditions on *other* variables. Eg: seq.a where seq.b > 1
        if queries and isinstance(candidate, SequenceType):
            # Filter candidate on the server-side, since the data may be 
            # proxied using ``dap.proxy.Proxy``.
            candidate = candidate.filter(*queries)
            # And then filter on the client side.
            candidate = filter_(candidate, queries)
        
        # If the variable was requested, either by id or name, or if no
        # variables were requested, we simply add this candidate to the
        # dataset we're building.
        if not vars_ or candidate.id in vars_ or (candidate.name in vars_ and candidate.name not in ids):
            new[candidate.name] = candidate

            # Check if requests done using shn are not ambiguous.
            if vars_ and candidate.id not in vars_:  # request by shn
                if candidate.name in names:
                    raise ConstraintExpressionError("Ambiguous shorthand notation request: %s" % candidate.name)
                names.append(candidate.name)

            # We also need to store the order in which the variables were
            # requested. Later, we'll rearrange the variables in our built
            # dataset in the correct order.
            if vars_:
                if candidate.id in vars_: index = vars_.index(candidate.id)
                else: index = vars_.index(candidate.name)
                order.append((index, candidate.name))

        # If the variable was not requested, but it's a constructor, it's
        # possible that one of its children has been requested. We apply
        # the algorithm recursively on the variable.
        elif not isinstance(var, BaseType):
            # We clear the candidate after storing a copy with the filtered
            # data and children. We will then append the requested children
            # to the cleared candidate.
            ccopy = copy.deepcopy(candidate)
            if isinstance(candidate, StructureType): 
                candidate.clear()
            else: 
                # If the variable is a grid we should return it as a 
                # structure with the requested fields.
                parent = candidate._id[:-len(candidate.name)-1]
                candidate = StructureType(name=candidate.name, attributes=candidate.attributes.copy())
                candidate._set_id(parent)

            # Check for requested children.
            candidate = build(ccopy, candidate, fields, queries, ids, names)

            # If the candidate has any keys, ie, stored variables, we add
            # it to the dataset we are building.
            if candidate.keys(): new[candidate.name] = candidate

        # Check if we need to apply a slice in the variable.
        slice_ = fields.get(candidate.id) or fields.get(candidate.name)
        if slice_: candidate = slicevar(candidate, slice_)
                
    # Sort variables according to order of requested variables.
    if len(order) > 1:
        order.sort()
        new._keys = [item[1] for item in order]
    
    return new


def filter_(dapvar, queries):
    # Get only the queries related to this variable.
    queries_ = [q for q in queries if q.startswith(dapvar.id)]
    if queries_:
        # Build the filter and apply it to the data.
        ids = [var.id for var in dapvar.values()]
        f = buildfilter(queries_, ids)
        data = itertools.ifilter(f, dapvar.data)

        # Set the data in the stored variables.
        data = list(data)
        dapvar.data = data

    return dapvar 


def slicevar(dapvar, slice_):
    if slice_ != (slice(None),):
        dapvar.data = dapvar.data[slice_]

    try:
        dapvar.shape = getattr(dapvar.data, 'shape', (len(dapvar.data),))
    except TypeError:
        pass

    if isinstance(dapvar, GridType):
        if not isiterable(slice_): slice_ = (slice_,)

        # Slice the maps.
        for map_,mapslice in zip(dapvar.maps.values(), slice_):
            map_.data = map_.data[mapslice]
            map_.shape = map_.data.shape

    return dapvar


def order(dataset, fields):
    """
    Order a given dataset according to the requested order.

        >>> d = DatasetType(name='d')
        >>> d['a'] = BaseType(name='a')
        >>> d['b'] = BaseType(name='b')
        >>> d['c'] = SequenceType(name='c')
        >>> d['c']['d'] = BaseType(name='d')
        >>> d['c']['e'] = BaseType(name='e')
        >>> print order(d, 'b,c.e,c.d,a'.split(','))  #doctest: +ELLIPSIS
        {'b': <dap.dtypes.BaseType object at ...>, 'c': {'e': <dap.dtypes.BaseType object at ...>, 'd': <dap.dtypes.BaseType object at ...>}, 'a': <dap.dtypes.BaseType object at ...>}
        >>> print order(d, 'c.e,c.d,a'.split(','))  #doctest: +ELLIPSIS
        {'c': {'e': <dap.dtypes.BaseType object at ...>, 'd': <dap.dtypes.BaseType object at ...>}, 'a': <dap.dtypes.BaseType object at ...>, 'b': <dap.dtypes.BaseType object at ...>}
        >>> print order(d, 'b,c,a'.split(','))  #doctest: +ELLIPSIS
        {'b': <dap.dtypes.BaseType object at ...>, 'c': {'d': <dap.dtypes.BaseType object at ...>, 'e': <dap.dtypes.BaseType object at ...>}, 'a': <dap.dtypes.BaseType object at ...>}
    """
    # Order the dataset.
    dataset = copy.copy(dataset)
    orders = []
    n = len(dataset._keys)
    for var in dataset.walk():
        # Search for id first.
        fields_ = [field[:len(var.id)] for field in fields]
        if var.id in fields_: index = fields_.index(var.id)
        # Else search by name.
        elif var.name in fields: index = fields.index(var.name)
        # Else preserve original order.
        else: index = n + dataset._keys.index(var.name)
        orders.append((index, var.name))

        # Sort children.
        if isinstance(var, StructureType):
            dataset[var.name] = order(var, fields)

    # Sort dataset.
    if len(orders) > 1:
        orders.sort()
        dataset._keys = [item[1] for item in orders]

    return dataset


def walk(dapvar):
    """
    Iterate over all variables, including dapvar.
    """
    yield dapvar
    try:
        for child in dapvar.walk():
            for var in walk(child): yield var
    except:
        pass


def getslice(hyperslab, shape=None):
    """Parse a hyperslab.

    Parse a hyperslab to a slice according to variable shape. The hyperslab
    follows the DAP specification, and ommited dimensions are returned in 
    their entirety.

        >>> getslice('[0:1:2][0:1:2]')
        (slice(0, 3, 1), slice(0, 3, 1))
        >>> getslice('[0:2][0:2]')
        (slice(0, 3, 1), slice(0, 3, 1))
        >>> getslice('[0][2]')
        (slice(0, 1, 1), slice(2, 3, 1))
        >>> getslice('[0:1:1]')
        (slice(0, 2, 1),)
        >>> getslice('[0:2:1]')
        (slice(0, 2, 2),)
        >>> getslice('')
        (slice(None, None, None),)
    """
    # Backwards compatibility. In pydap <= 2.2.3 the ``fields`` dict from 
    # helper.parse_querystring returned the slices as strings (instead of
    # Python slices). These strings had to be passed to getslice to get a 
    # Python slice. Old plugins still do this, but with pydap >= 2.2.4 
    # they are already passing the slices, so we simply return them.
    if not isinstance(hyperslab, basestring): return hyperslab or slice(None)

    if hyperslab:
        output = []
        dimslices = hyperslab[1:-1].split('][')
        for dimslice in dimslices:
            start, size, step = _getsize(dimslice)
            output.append(slice(start, start+size, step))
        output = tuple(output)
    else:
        output = (slice(None),)

    return output


def _getsize(dimslice):
    """Parse a dimension from a hyperslab.

    Calculates the start, size and step from a DAP formatted hyperslab.

        >>> _getsize('0:1:9')
        (0, 10, 1)
        >>> _getsize('0:2:9')
        (0, 10, 2)
        >>> _getsize('0')
        (0, 1, 1)
        >>> _getsize('0:9')
        (0, 10, 1)
    """
    size = dimslice.split(':')

    start = int(size[0])
    if len(size) == 1:
        stop = start
        step = 1
    elif len(size) == 2:
        stop = int(size[1])
        step = 1
    elif len(size) == 3:
        step = int(size[1])
        stop = int(size[2])
    else:
        raise ConstraintExpressionError('Invalid hyperslab: %s.' % dimslice)
    size = (stop-start) + 1
    
    return start, size, step


def buildfilter(queries, vars_):
    """This function is a filter builder.

    Given a list of DAP formatted queries and a list of variable names,
    this function returns a dynamic filter function to filter rows.

    From the example in the DAP specification:

        >>> vars_ = ['index', 'temperature', 'site']
        >>> data = []
        >>> data.append([10, 17.2, 'Diamond_St'])
        >>> data.append([11, 15.1, 'Blacktail_Loop'])
        >>> data.append([12, 15.3, 'Platinum_St'])
        >>> data.append([13, 15.1, 'Kodiak_Trail'])

    Rows where index is greater-than-or-equal 11:

        >>> f = buildfilter(['index>=11'], vars_)
        >>> for line in itertools.ifilter(f, data):
        ...     print line
        [11, 15.1, 'Blacktail_Loop']
        [12, 15.300000000000001, 'Platinum_St']
        [13, 15.1, 'Kodiak_Trail']

    Rows where site ends with '_St':

        >>> f = buildfilter(['site=~".*_St"'], vars_)
        >>> for line in itertools.ifilter(f, data):
        ...     print line
        [10, 17.199999999999999, 'Diamond_St']
        [12, 15.300000000000001, 'Platinum_St']

    Index greater-or-equal-than 11 AND site ends with '_St':

        >>> f = buildfilter(['site=~".*_St"', 'index>=11'], vars_)
        >>> for line in itertools.ifilter(f, data):
        ...     print line
        [12, 15.300000000000001, 'Platinum_St']

    Site is either 'Diamond_St' OR 'Blacktail_Loop':
    
        >>> f = buildfilter(['site={"Diamond_St", "Blacktail_Loop"}'], vars_)
        >>> for line in itertools.ifilter(f, data):
        ...     print line
        [10, 17.199999999999999, 'Diamond_St']
        [11, 15.1, 'Blacktail_Loop']

    Index is either 10 OR 12:

        >>> f = buildfilter(['index={10, 12}'], vars_)
        >>> for line in itertools.ifilter(f, data):
        ...     print line
        [10, 17.199999999999999, 'Diamond_St']
        [12, 15.300000000000001, 'Platinum_St']

    Python is great, isn't it? :)
    """
    filters = []
    p = re.compile(r'''^                          # Start of selection
                       {?                         # Optional { for multi-valued constants
                       (?P<var1>.*?)              # Anything
                       }?                         # Closing }
                       (?P<op><=|>=|!=|=~|>|<|=)  # Operators
                       {?                         # {
                       (?P<var2>.*?)              # Anything
                       }?                         # }
                       $                          # EOL
                    ''', re.VERBOSE)
    for query in queries:
        m = p.match(query)
        if not m: raise ConstraintExpressionError('Invalid constraint expression: %s.' % query)

        # Functions associated with each operator.
        op = {'<' : operator.lt,
              '>' : operator.gt,
              '!=': operator.ne,
              '=' : operator.eq,
              '>=': operator.ge,
              '<=': operator.le,
              '=~': lambda a,b: re.match(b,a),
             }[m.group('op')]
        # Allow multiple comparisons in one line. Python rulez!
        op = multicomp(op)

        # Build the filter for the first variable.
        if m.group('var1') in vars_:
            i = vars_.index(m.group('var1'))
            var1 = lambda L, i=i: operator.getitem(L, i)

            # Build the filter for the second variable. It could be either
            # a name or a constant.
            if m.group('var2') in vars_:
                i = vars_.index(m.group('var2'))
                var2 = lambda L, i=i: operator.getitem(L, i)
            else:
                var2 = lambda x, m=m: expr_eval(m.group('var2'))

            # This is the filter. We apply the function (op) to the variable
            # filters (var1 and var2).
            filter0 = lambda x, op=op, var1=var1, var2=var2: op(var1(x), var2(x))
            filters.append(filter0)

    if filters:
        # We have to join all the filters that were built, using the AND 
        # operator. Believe me, this line does exactly that.
        #
        # You are not expected to understand this.
        filter0 = lambda i: reduce(lambda x,y: x and y, [f(i) for f in filters])
    else:
        filter0 = bool

    return filter0


def multicomp(function):
    """Multiple OR comparisons.

    Given f(a,b), this function returns a new function g(a,b) which
    performs multiple OR comparisons if b is a tuple.

        >>> a = 1
        >>> b = (0, 1, 2)
        >>> operator.lt = multicomp(operator.lt)
        >>> operator.lt(a, b)
        True
    """
    def f(a, b):
        if isinstance(b, tuple):
            for i in b:
                # Return True if any comparison is True.
                if function(a, i): return True
            return False
        else:
            return function(a, b)

    return f


def fix_slice(dims, index):
    """Fix incomplete slices or slices with ellipsis.

    The behaviour of this function was reversed-engineered from numpy.

        >>> fix_slice(3, (0, Ellipsis, 0))
        (0, slice(None, None, None), 0)
        >>> fix_slice(4, (0, Ellipsis, 0))
        (0, slice(None, None, None), slice(None, None, None), 0)
        >>> fix_slice(4, (0, 0, Ellipsis, 0))
        (0, 0, slice(None, None, None), 0)
        >>> fix_slice(5, (0, Ellipsis, 0))
        (0, slice(None, None, None), slice(None, None, None), slice(None, None, None), 0)
        >>> fix_slice(5, (0, 0, Ellipsis, 0))
        (0, 0, slice(None, None, None), slice(None, None, None), 0)
        >>> fix_slice(5, (0, Ellipsis, 0, Ellipsis))
        (0, slice(None, None, None), slice(None, None, None), 0, slice(None, None, None))
        >>> fix_slice(4, slice(None, None, None))
        (slice(None, None, None), slice(None, None, None), slice(None, None, None), slice(None, None, None))
        >>> fix_slice(4, (slice(None, None, None), 0))
        (slice(None, None, None), 0, slice(None, None, None), slice(None, None, None))
    """
    if not isinstance(index, tuple): index = (index,)

    out = []
    length = len(index)
    for slice_ in index:
        if slice_ is Ellipsis:
            out.extend([slice(None)] * (dims - length + 1))
            length += (dims - length)
        else:
            out.append(slice_)
    index = tuple(out)

    if len(index) < dims:
        index += (slice(None),) * (dims - len(index))

    return index


def lenslice(slice_):
    """
    Return the number of values associated with a slice.
    
    By Bob Drach.
    """
    step = slice_.step
    if step is None: step = 1

    if step > 0:
        start = slice_.start
        stop = slice_.stop
    else:
        start = slice_.stop
        stop = slice_.start
        step = -step
    return ((stop-start-1)/step + 1)


def parse_querystring(query):
    """
    Parse a query_string returning the requested variables, dimensions, and CEs.

        >>> parse_querystring('a,b')
        ({'a': (slice(None, None, None),), 'b': (slice(None, None, None),)}, [])
        >>> parse_querystring('a[0],b[1]')
        ({'a': (slice(0, 1, 1),), 'b': (slice(1, 2, 1),)}, [])
        >>> parse_querystring('a[0],b[1]&foo.bar>1')
        ({'a': (slice(0, 1, 1),), 'b': (slice(1, 2, 1),)}, ['foo.bar>1'])
        >>> parse_querystring('a[0],b[1]&foo.bar>1&LAYERS=SST')
        ({'a': (slice(0, 1, 1),), 'b': (slice(1, 2, 1),)}, ['foo.bar>1', 'LAYERS=SST'])
        >>> parse_querystring('foo.bar>1&LAYERS=SST')
        ({}, ['foo.bar>1', 'LAYERS=SST'])

    """
    if query is None:  return {}, []

    query = unquote(query)
    constraints = query.split('&')

    # Check if the first item is either a list of variables (projection)
    # or a selection.
    relops = ['<', '<=', '>', '>=', '=', '!=',' =~']
    for relop in relops:
        if relop in constraints[0]:
            vars_ = []
            queries = constraints[:]
            break
    else:
        vars_ = constraints[0].split(',')
        queries = constraints[1:]
    
    fields = odict()
    p = re.compile(r'(?P<name>[^[]+)(?P<shape>(\[[^\]]+\])*)')
    for var in vars_:
        if var:
            # Check if the var has a slice.
            c = p.match(var).groupdict()
            id_ = quote(c['name'])
            fields[id_] = getslice(c['shape'])

    return fields, queries


def escape_dods(dods, pad=''):
    """
    Escape a DODS response.

    This is useful for debugging. You're probably spending too much time
    with pydap if you need to use this.
    """
    if 'Data:\n' in dods:
        index = dods.index('Data:\n') + len('Data:\n')
    else:
        index = 0

    dds = dods[:index]
    dods = dods[index:]

    out = []
    for i, char in enumerate(dods): 
        char = hex(ord(char))
        char = char.replace('0x', '\\x')
        if len(char) < 4: char = char.replace('\\x', '\\x0')
        out.append(char)
        if pad and (i%4 == 3): out.append(pad)
    out = ''.join(out)
    out = out.replace(r'\x5a\x00\x00\x00', '<start of sequence>')
    out = out.replace(r'\xa5\x00\x00\x00', '<end of sequence>\n')
    return dds + out


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
