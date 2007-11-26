"""Ordered dictionary.

This is a dictionary class that preserves the order in which the keys are
stored. This is necessary to build Structures and Sequences that follow the
requested variable order.
"""

__author__ = "Roberto De Almeida <rob@pydap.org>"

import copy


class odict(dict):

    """Ordered dictionary."""
    
    def __init__(self, odict=None):
        self._keys = []
        self._dict = {}

        if odict is not None:
            self.update(odict)

    def __iter__(self):
        return iter(self._keys[:])

    def __setitem__(self, key, item):
        self._dict.__setitem__(key, item)
        if key not in self._keys: self._keys.append(key)

    def __getitem__(self, key):
        return self._dict.__getitem__(key)

    def __delitem__(self, key):
        self._dict.__delitem__(key)
        self._keys.remove(key)

    def keys(self):
        return self._keys[:]

    def items(self):
        return [(key, self._dict.__getitem__(key)) for key in self._keys]

    def values(self):
        return [self._dict.__getitem__(key) for key in self._keys]

    def iterkeys(self):
        for key in self._keys: yield key

    def iteritems(self):
        for key in self._keys: yield (key, self._dict.__getitem__(key))

    def itervalues(self):
        for key in self._keys: yield self._dict.__getitem__(key)

    def clear(self):
        self._dict.clear()
        self._keys = []

    def copy(self):
        new = odict(self)
        return new

    def update(self, odict):
        for k, v in odict.items():
            self.__setitem__(k, v)

    def setdefault(self, key, d=None):
        if key not in self._keys: self._keys.append(key)
        return self._dict.setdefault(key, d)

    def get(self, key, d=None):
        if key in self._keys: return self._dict.__getitem__(key)
        else: return d

    def has_key(self, key):
        return self._dict.has_key(key)

    def popitem(self):
        try: key = self._keys[-1]
        except IndexError: raise KeyError('dictionary is empty')

        self._keys.remove(key)
        return self._dict.pop(key)

    def pop(self, key, d=None):
        value = self._dict.pop(key, d)

        try: self._keys.remove(key)
        except ValueError: pass
        
        return value

    def fromkeys(keys, d=None):
        new = odict()
        for key in keys: new.__setitem__(key, d)
        return new

    def __contains__(self, key):
        return self._dict.has_key(key)

    def __len__(self):
        return self._dict.__len__()

    def __repr__(self):
        return '{%s}' % ', '.join(['%s: %s' % (k.__repr__(), v.__repr__()) for (k, v) in self.items()])

    def __str__(self):
        return '{%s}' % ', '.join(['%s: %s' % (k.__repr__(), v.__repr__()) for (k, v) in self.items()])

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo=None, _nil=[]):
        new = odict()
        for k, v in self.items():
            new.__setitem__(k, copy.deepcopy(v))
        return new 


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
