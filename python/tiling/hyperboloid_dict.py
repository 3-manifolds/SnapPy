from .real_hash_dict import RealHashDict
from ..hyperboloid import r13_dot
from ..exceptions import InsufficientPrecisionError # type: ignore

def get_hyperboloid_dict(max_neg_prod_equal, min_neg_prod_distinct,
                         verified):
    RF = max_neg_prod_equal.parent()
    return RealHashDict(
        _equality_predicate(max_neg_prod_equal, min_neg_prod_distinct),
        _hash_function(RF),
        _compute_epsilon_inverse(RF),
        verified)

def _equality_predicate(max_neg_prod_equal, min_neg_prod_distinct):
    def result(point_0, point_1):
        neg_inner_product = -r13_dot(point_0, point_1)
        if neg_inner_product < max_neg_prod_equal:
            return True
        if neg_inner_product > min_neg_prod_distinct:
            return False

        raise InsufficientPrecisionError(
            "Could neither verify that the two given tiles are "
            "the same nor that they are distinct.\n"
            "Neg inner product is:     %r.\n"
            "Max to be same tile:      %r.\n"
            "Min to be distinct tiles: %r." % (
                neg_inner_product, max_neg_prod_equal, min_neg_prod_distinct))

    return result

def _hash_function(RF):
    weights = [ RF(1.2003), RF(0.94553), RF(1.431112), RF(1.2342) ]

    def result(point):
        return (point[0] * weights[0] +
                point[1] * weights[1] +
                point[2] * weights[2] +
                point[3] * weights[3])

    return result

_default_epsilon_inverse = 128

def _compute_epsilon_inverse(RF):
    result = _default_epsilon_inverse

    p = RF.precision()
    if p > 53:
        result *= 2 ** ((p - 53) // 3)

    return result


