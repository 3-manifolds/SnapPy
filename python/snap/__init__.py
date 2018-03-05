from .shapes import polished_tetrahedra_shapes
from ..sage_helper import _within_sage, sage_method
from .polished_reps import polished_holonomy
from . import nsagetools, interval_reps

if _within_sage:
    from .find_field import ListOfApproximateAlgebraicNumbers

@sage_method
def tetrahedra_field_gens(manifold):
    """
    The shapes of the tetrahedra as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m015')
        sage: tets = M.tetrahedra_field_gens()
        sage: tets.find_field(100, 10, optimize=True)    # doctest: +NORMALIZE_WHITESPACE
        (Number Field in z with defining polynomial x^3 - x - 1,
        <ApproxAN: -0.662358978622 - 0.562279512062*I>, [-z, -z, -z])
    """
    if manifold.is_orientable():
        def func(prec):
            return polished_tetrahedra_shapes(manifold, bits_prec=prec)
    else:
        double_cover = manifold.orientation_cover()
        def func(prec):
            return polished_tetrahedra_shapes(double_cover, bits_prec=prec)[::2]
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def trace_field_gens(manifold, fundamental_group_args = []):
    """
    The generators of the trace field as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m125')
        sage: traces = M.trace_field_gens()
        sage: traces.find_field(100, 10, optimize=True)    # doctest: +NORMALIZE_WHITESPACE
        (Number Field in z with defining polynomial x^2 + 1,
        <ApproxAN: -1.0*I>, [z + 1, z, z + 1])
    """
    def func(prec):
        return polished_holonomy(manifold, prec,
                                 fundamental_group_args).trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def invariant_trace_field_gens(manifold, fundamental_group_args = []):
    """
    The generators of the trace field as ApproximateAlgebraicNumbers. Can be
    used to compute the tetrahedra field, where the first two parameters
    are bits of precision and maximum degree of the field::

        sage: M = Manifold('m007(3,1)')
        sage: K = M.invariant_trace_field_gens().find_field(100, 10, optimize=True)[0]
        sage: L = M.trace_field_gens().find_field(100, 10, optimize=True)[0]
        sage: K.polynomial(), L.polynomial()
        (x^2 - x + 1, x^4 - 2*x^3 + x^2 + 6*x + 3)
    """
    def func(prec):
        return polished_holonomy(manifold, prec,
                                 fundamental_group_args).invariant_trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

@sage_method
def holonomy_matrix_entries(manifold, fundamental_group_args = []):
    def func(prec):
        G = polished_holonomy(manifold, prec, fundamental_group_args)
        return sum( [G.SL2C(g).list() for g in G.generators()], [])
    return ListOfApproximateAlgebraicNumbers(func)


def add_methods(mfld_class, hyperbolic=True):
    mfld_class.alexander_polynomial = nsagetools.alexander_polynomial
    mfld_class.homological_longitude = nsagetools.homological_longitude
    if hyperbolic:
        mfld_class.polished_holonomy = polished_holonomy
        mfld_class.tetrahedra_field_gens = tetrahedra_field_gens
        mfld_class.trace_field_gens = trace_field_gens
        mfld_class.invariant_trace_field_gens = invariant_trace_field_gens
        mfld_class.holonomy_matrix_entries = holonomy_matrix_entries
        mfld_class.hyperbolic_torsion = nsagetools.hyperbolic_torsion
        mfld_class.hyperbolic_adjoint_torsion = nsagetools.hyperbolic_adjoint_torsion
        mfld_class.hyperbolic_SLN_torsion = nsagetools.hyperbolic_SLN_torsion
