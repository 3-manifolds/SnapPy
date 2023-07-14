from .lifted_tetrahedron import LiftedTetrahedron
from .real_hash_dict import RealHashDict
from .canonical_key_dict import CanonicalKeyDict

from ..snap.t3mlite import Mcomplex # type: ignore
from ..hyperboloid import r13_dot, o13_inverse

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

    def __init__(self, dictionary, base_point, act_on_base_point_by_inverse):
        self._set = ProductSet(dictionary)
        self._base_point = base_point
        self._act_on_base_point_by_inverse = act_on_base_point_by_inverse

    def add(self, lifted_tetrahedron : LiftedTetrahedron) -> bool:
        if self._act_on_base_point_by_inverse:
            m = o13_inverse(lifted_tetrahedron.o13_matrix)
        else:
            m = lifted_tetrahedron.o13_matrix
        
        return self._set.add(
            m * self._base_point,
            lifted_tetrahedron.tet.Index)

def get_lifted_tetrahedron_set(base_point,
                               canonical_keys_function,
                               act_on_base_point_by_inverse,
                               min_inner_product,
                               verified
                               ) -> LiftedTetrahedronSet:
    """
    Returns a set to store lifted tetrahedra in H^3 or a quotient
    space of H^3.

    The type of geometric object determined which space will be used.

    """

    d = RealHashDict(
        _equality_predicate(min_inner_product),
        _hash_function(base_point[0].parent()),
        _epsilon_inverse,
        verified)

    if canonical_keys_function:
        d = CanonicalKeyDict(d, canonical_keys_function)

    return LiftedTetrahedronSet(d, base_point, act_on_base_point_by_inverse)

def _equality_predicate(min_inner_product):
    def result(point_0, point_1):
        inner_product = r13_dot(point_0, point_1)
        if inner_product > min_inner_product:
            return True
        if inner_product < min_inner_product:
            return False

        raise InsufficientPrecisionError(
            "Could neither verify that the two given tiles are "
            "the same nor that they are distinct. "
            "Inner product is: %r, cut-off is: %s. " % (
                inner_product, min_inner_product))

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

