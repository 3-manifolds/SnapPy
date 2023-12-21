from ..verify.interval_tree import IntervalTree
from .floor import floor_as_integers

# For experimenting how the performance of the IntervalTree compares
# to using a python dict.
#
# Sadly, the IntervalTree is implemented in slow python and not fast
# C(++). This is because we do not have access to the mpfi number since
# we cannot include .pyd files from SageMath :(
#
_use_interval_tree = True

class RealHashDict:
    """
    A dictionary (following python dict API) where we provide for the keys:
    - a predicate that decides whether two keys are the same
    - a hash function returning an approximation of a real number

    The use case is, for example, to use as keys points in a metric space
    that are given by approximate numbers. Two close by approximations are
    probably the same point.

    The requirements of the hash function has the following requirements:
    - If it is an interval, it has to contain the true hash value, that is
      the true value of the hash function evaluated for the true value of
      the key.
    - If it is an interval and _use_interval_tree is False, the length of
      the interval cannot be much larger than 1/epsilon_inverse.
    - If it is a floating point number, the approximation cannot differ
      from the true value by more than 1/epsilon_inverse/16.

    Example Z as a discrete subset of R::

    >>> from snappy.number import *
    >>> def to_number(py_float):
    ...     return number_to_native_number(SnapPyNumbers(100)(py_float))

    >>> d = RealHashDict(
    ...        # If two floating approximations of integers are close by,
    ...        # they must be the same
    ...        equality_predicate = lambda x, y: abs(x - y) < to_number(0.5),
    ...        # For space with more than 1-dimension, use some function
    ...        # involving the other coordinates
    ...        hash_function = lambda x : x,
    ...        epsilon_inverse = 1024,
    ...        verified = False)

    # Keys near 1 should be treated as the same.
    >>> d[to_number(1.0)] = 'A'
    >>> d[to_number(1.0000001)] = 'B'
    >>> d.get(to_number(0.9999999))
    'B'

    >>> d.setdefault(to_number(10 + 0.1 / 1024), []).append('A')
    >>> d.get(to_number(10.0))
    ['A']
    >>> d.setdefault(to_number(10 - 0.1 / 1024), []).append('B')
    >>> d.get(to_number(10.0000001))
    ['A', 'B']

    sage: from sage.all import RIF
    sage: def equality_predicate(x, y):
    ...       d = abs(x - y)
    ...       if d < RIF(0.1):
    ...           return True
    ...       if d > 0:
    ...           return False
    ...       raise Exception("Could not determine whether points are the "
    ...                       "same or not")
    sage: d = RealHashDict(
    ...        equality_predicate = equality_predicate,
    ...        hash_function = lambda x : x,
    ...        epsilon_inverse = 1024,
    ...        verified = True)

    # Keys near 1 should be treated as the same.
    sage: d[RIF(0.9999999, 1.000001)] = 'A'
    sage: d[RIF(0.9999998, 1.0000001)] = 'B'
    sage: d.get(RIF(0.99999997, 1.0000002))
    'B'
    
    sage: d.setdefault(RIF(9.9999999, 10.00001), []).append('A')
    sage: d.get(RIF(10.0))
    ['A']
    sage: d.setdefault(RIF(9.9999993, 10.00002), []).append('B')
    sage: d.get(RIF(9.99999999, 10.0000001))
    ['A', 'B']

    sage: from snappy.tiling import real_hash_dict
    sage: original_use_interval_tree = real_hash_dict._use_interval_tree
    sage: real_hash_dict._use_interval_tree = True
    sage: d = RealHashDict(
    ...        equality_predicate = equality_predicate,
    ...        hash_function = lambda x : x,
    ...        epsilon_inverse = 1024,
    ...        verified = True)
    sage: d[RIF(1.0)] = 'A'
    sage: d[RIF(0.7, 1.3)] = 'B'
    Traceback (most recent call last):
    ...
    Exception: Could not determine whether points are the same or not
    sage: real_hash_dict._use_interval_tree = original_use_interval_tree
    
    
    """

    def __init__(self, equality_predicate, hash_function, epsilon_inverse, verified):

        self._equality_predicate = equality_predicate
        self._hash_function = hash_function

        if verified and _use_interval_tree:
            self._impl = IntervalTree()
        else:
            self._impl = _Dict(epsilon_inverse)

    def get(self, key, default = None):
        h = self._hash_function(key)

        for other_key, value in self._impl.find(h):
            if self._equality_predicate(key, other_key):
                return value

        return default

    # Not implemented because it was not needed yet:
    # __getitem__

    def __setitem__(self, key, value):
        h = self._hash_function(key)

        for pair in self._impl.find(h):
            if self._equality_predicate(key, pair[0]):
                pair[1] = value
                return
        
        self._impl.insert(h, [key, value])

    def setdefault(self, key, default):
        h = self._hash_function(key)

        for other_key, value in self._impl.find(h):
            if self._equality_predicate(key, other_key):
                return value

        self._impl.insert(h, [key, default])

        return default

class _Dict:
    def __init__(self, epsilon_inverse):
        self._epsilon_inverse = epsilon_inverse
        self._dict = {}

    def insert(self, key, value):
        for computed_key in self._computed_keys(key):
            self._dict.setdefault(computed_key, []).append(value)

    def find(self, key):
        for computed_key in self._computed_keys(key):
            yield from self._dict.get(computed_key, [])

    def _computed_keys(self, key):
        return floor_as_integers(key * self._epsilon_inverse)
