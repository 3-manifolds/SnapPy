from .real_hash_dict import RealHashDict
from ..hyperboloid import r13_dot
from ..exceptions import InsufficientPrecisionError # type: ignore

_epsilon_inverse = 1024

def get_hyperboloid_dict(min_inner_product, verified):
    return RealHashDict(
        _equality_predicate(min_inner_product),
        _hash_function(min_inner_product.parent()),
        _epsilon_inverse,
        verified)

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
