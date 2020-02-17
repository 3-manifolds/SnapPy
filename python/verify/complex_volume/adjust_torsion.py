from ...sage_helper import _within_sage, sage_method
if _within_sage:
    from sage.all import pi
    import sage.all

from .extended_bloch import *
from ...snap import t3mlite as t3m

__all__ = ['verified_complex_volume_from_lifted_ptolemys',
           'normalize_by_pi_square_over_two']

_move_to_three = {
    t3m.F0 : t3m.Perm4((3,0,1,2)),
    t3m.F1 : t3m.Perm4((0,3,1,2)),
    t3m.F2 : t3m.Perm4((0,1,3,2)),
    t3m.F3 : t3m.Perm4((0,1,2,3))    
}

_move_from_three = {
    k : ~p for k, p in _move_to_three.items()
}

def _perm_for_q_tet(F, gluing):
    return _move_to_three[gluing.image(F)] * gluing * _move_from_three[F]

def _compute_adjustment_for_face(face):
    canonical_corner = face.Corners[0]
    tet     = canonical_corner.Tetrahedron
    F       = canonical_corner.Subsimplex
    gluing  = tet.Gluing[F]
    other_F = gluing.image(F)
    
    return -2 * _perm_for_q_tet(F, gluing)[0] * (-1) ** t3m.FaceIndex[other_F]

def _compute_adjustment(mcomplex):
    """
    Given an mcomplex, compute the adjustment term to account for the
    triangulation not being ordered.
    
    So far, only solves for the 3-torsion but 2-torsion remains
    """
    return sum([ _compute_adjustment_for_face(face)
                 for face in mcomplex.Faces ])

@sage_method
def verified_complex_volume_from_lifted_ptolemys(mcomplex, ptolemys):
    """
    Given lifted Ptolemy coordinates for a triangulation (as dictionary)
    and the number of tetrahedra, compute the complex volume (where
    the real part is the Chern-Simons and the imaginary part is the
    volume).

    The result is correct modulo pi^2/2.
    """

    # Simply add Neumann's dilog over all simplicies as if the triangulation
    # was ordered
    result = compute_complex_volume_from_lifted_ptolemys_no_torsion_adjustment(
        len(mcomplex.Tetrahedra), ptolemys)

    # Add suitable multiple of pi^2/6 to the result to account for the fact
    # that the triangulation was probably not ordered
    CIF = result.parent()
    return result + _compute_adjustment(mcomplex) * CIF(pi ** 2 / 6)


@sage_method
def normalize_by_pi_square_over_two(z):
    """
    Add multiples of pi^2/2 to the real part to try to bring the
    real part between -pi^2/4 and pi^2/4.
    """

    CIF = z.parent()
    RIF = CIF.real_field()

    pi_square_over_two = RIF(pi**2/2)

    # Round to integer
    q = (z.real().center() / pi_square_over_two.center()).round()
    
    # Subtract multiple of pi^2/6
    return z - q * pi_square_over_two

