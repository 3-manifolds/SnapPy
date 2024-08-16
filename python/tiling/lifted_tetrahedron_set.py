from .lifted_tetrahedron import LiftedTetrahedron
from .canonical_key_dict import CanonicalKeyDict
from .dict_based_set import DictBasedProductSet
from .hyperboloid_dict import get_hyperboloid_dict

from ..hyperboloid import o13_inverse

class LiftedTetrahedronSet:
    """
    A set of lifted tetrahedra in H^3 or a quotient space of
    H^3.

    Which space will be used is determined by the dictionary
    given when LiftedTetrahedronSet is constructed.
    """

    def __init__(self, dictionary, base_point, act_on_base_point_by_inverse):
        self._set = DictBasedProductSet(dictionary)
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
                               max_neg_prod_equal, min_neg_prod_distinct,
                               verified
                               ) -> LiftedTetrahedronSet:
    """
    Returns a set to store lifted tetrahedra in H^3 or a quotient
    space of H^3.

    The type of geometric object determined which space will be used.

    """

    d = get_hyperboloid_dict(max_neg_prod_equal, min_neg_prod_distinct,
                             verified)

    if canonical_keys_function:
        d = CanonicalKeyDict(d, canonical_keys_function)

    return LiftedTetrahedronSet(d, base_point, act_on_base_point_by_inverse)


