from ..matrix import vector

__all__ = ['Infinity', 'ideal_point_to_r13']

Infinity = 'Infinity'


def ideal_point_to_r13(z, RF):
    """
    Takes a boundary point z of H^3, that is either a complex number
    or Infinite, and a real field type.

    Returns the corresponding unit light vector in the hyperboloid model.
    """
    if z == Infinity:
        return vector([RF(1), RF(1), RF(0), RF(0)])

    z_re = z.real()
    z_im = z.imag()
    z_abs_sqr = z_re**2 + z_im**2
    denom = z_abs_sqr + 1

    return vector([RF(1),
                   (z_abs_sqr - 1) / denom,
                   2 * z_re / denom,
                   2 * z_im / denom])
