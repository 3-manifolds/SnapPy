from .line import R13LineWithMatrix
from .lifted_tetrahedron import LiftedTetrahedron
from .real_hash_dict import RealHashDict
from .canonical_key_dict import CanonicalKeyDict
from .canonical_keys import canonical_keys_function_for_line

from ..snap.t3mlite import Mcomplex # type: ignore
from ..hyperboloid import distance_unit_time_r13_points

from typing import Union

class ProductSet:
    """
    Behaves like a set of elements in the product AxB.
    It is implemented using a dictionary mapping A->B.

    Can be combined with RealHashDict or CanonicalKeyDict.
    """

    def __init__(self, dictionary):
        self._dictionary = dictionary

    def add(self, key, element):
        elements = self._get_elements_for_key(key)
        if element in elements:
            return False
        else:
            elements.add(element)
            return True

    def _get_elements_for_key(self, key):
        """
        Could be used for potential future optimization:
        When going to a neighboring tetrahedron paired through
        a trivial face-pairing matrix (so it is a tetrahedron within
        the same copy of the fundamental polyhedron), we can keep
        a reference to the set returned by _get_elements_for_key
        and just add the neighboring tetrahedron to this set rather
        than doing the look-up from scratch.
        """
        return self._dictionary.setdefault(key, set())

class LiftedTetrahedronSet:
    """
    A set of lifted tetrahedra in H^3 or a quotient space of
    H^3.

    Which space will be used is determined by the dictionary
    given when LiftedTetrahedronSet is constructed.
    """

    def __init__(self, dictionary, base_point):
        self._set = ProductSet(dictionary)
        self._base_point = base_point

    def add(self, lifted_tetrahedron : LiftedTetrahedron) -> bool:
        return self._set.add(
            lifted_tetrahedron.o13_matrix * self._base_point,
            lifted_tetrahedron.tet.Index)

def get_lifted_tetrahedron_set(mcomplex : Mcomplex,
                               geometric_object : Union[R13LineWithMatrix]
                               ) -> LiftedTetrahedronSet:
    """
    Returns a set to store lifted tetrahedra in H^3 or a quotient
    space of H^3.

    The type of geometric object determined which space will be used.

    """

    d = RealHashDict(
        _equality_predicate(mcomplex, geometric_object),
        _hash_function(mcomplex.RF),
        _epsilon_inverse,
        mcomplex.verified)

    if isinstance(geometric_object, R13LineWithMatrix):
        d = CanonicalKeyDict(
            d, canonical_keys_function_for_line(geometric_object))

    if isinstance(geometric_object, R13LineWithMatrix):
        base_point = mcomplex.R13_baseTetInCenter

    return LiftedTetrahedronSet(d, base_point)

def _equality_predicate(mcomplex, geometric_object):
    min_distance = mcomplex.baseTetInRadius

    if mcomplex.verified:
        right_dist = min_distance
        left_dist = 0
    else:
        RF = min_distance.parent()
        right_dist = min_distance * RF(0.125)
        left_dist = min_distance * RF(0.5)

    def result(point_0, point_1):
        d = distance_unit_time_r13_points(point_0, point_1)
        if d < right_dist:
            return True
        if d > left_dist:
            return False

        raise InsufficientPrecisionError(
            "Could neither verify that the two given tiles are "
            "the same nor that they are distinct. "
            "Distance between basepoint translates is: %r. "
            "Cut-offs are %r %r." % (
            d, left_dist, right_dist))

    return result

def _hash_function(RF):
    weights = [ RF(1.2003), RF(0.94553), RF(1.431112), RF(1.2342) ]

    def result(point):
        return (point[0] * weights[0] +
                point[1] * weights[1] +
                point[2] * weights[2] +
                point[3] * weights[3])

    return result

_epsilon_inverse = 1024

