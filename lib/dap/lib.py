from __future__ import division

"""Basic functions concerning the DAP.

These functions are mostly related to encoding data according to the DAP.
"""

from urllib import quote as _quote

__author__ = 'Roberto De Almeida <rob@pydap.org>'
__version__ = (2,2,6,2)  # module version
__dap__ = (2,0)        # protocol version

# Constants that used to live in __init__.py but had to be moved
# because we share the namespace with plugins and responses.
USER_AGENT = 'pydap/%s' % '.'.join([str(_) for _ in __version__])
INDENT = ' ' * 4
VERBOSE = False
CACHE = None
TIMEOUT = None


def isiterable(o):
    """Tests if an object is iterable.

        >>> print isiterable(range(10))
        True
        >>> print isiterable({})
        True
        >>> def a():
        ...     for i in range(10): yield i
        >>> print isiterable(a())
        True
        >>> print isiterable('string')
        False
        >>> print isiterable(1)
        False
    """
    # We DON'T want to iterate over strings.
    if isinstance(o, basestring): return False

    try:
        iter(o)
        return True
    except TypeError:
        return False


def to_list(L):
    if hasattr(L, 'tolist'): return L.tolist()  # shortcut for numpy arrays
    elif isiterable(L): return [to_list(item) for item in L]
    else: return L


def quote(name):
    """Extended quote for the DAP spec.

    The period MUST be escaped in names (DAP spec, item 5.1):

        >>> quote("White space")
        'White%20space'
        >>> _quote("Period.")
        'Period.'
        >>> quote("Period.")
        'Period%2E'
    """
    return _quote(name).replace('.', '%2E')


def encode_atom(atom):
    r"""Atomic types encoding.

    Encoding atomic types for the DAS. Integers should be printed using the
    base 10 ASCII representation of its value:

        >>> encode_atom(42)
        '42'

    Floating point values are printed using the base 10 ASCII
    representation of its value, conforming to ANSI C's description of
    printf using the %g format and precision 6.

        >>> encode_atom(1/3)
        '0.333333'

    String and URLs are printed in ASCII, escaped according to escape().

        >>> encode_atom('This is a "string".')
        '"This is a \\"string\\"."'
    """
    return {basestring: lambda s: escape(s),
            unicode   : lambda s: escape(s),
            str       : lambda s: escape(s),
            float     : lambda f: '%.6g' % f,
            long      : lambda f: '%.6g' % f,
            int       : lambda i: repr(i),
           }.get(type(atom), lambda obj: '%s' % obj)(atom)


def escape(s):
    r"""Escape a string.

    Enclose strings with double quotes, escape double quotes and backslashes:

        >>> escape('String')
        '"String"'
        >>> escape('This is a "test".')
        '"This is a \\"test\\"."'
    """
    s = s.replace(r'\\', r'\\\\')
    s = s.replace(r'"', r'\"')
    s = '"%s"' % s

    return s


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
