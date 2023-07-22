from .exceptions import *

from ..snap.t3mlite import simplex

def check_polynomial_edge_equations_exactly(mcomplex):
    """
    Check that the polynomial edge equations are fulfilled exactly.

    We use the conjugate inverse to support non-orientable manifolds.
    """

    # For each edge
    for edge in mcomplex.Edges:
        # The exact value when evaluating the edge equation
        val = 1

        # Iterate through edge embeddings
        for tet, perm in edge.embeddings():
            # Accumulate shapes of the edge exactly
            val *= _shape_for_edge_embedding(tet, perm)
            
        if not val == 1:
            raise EdgeEquationExactVerifyError(val)

def check_logarithmic_edge_equations_and_positivity(mcomplex, NumericalField):
    """
    Check that the shapes have positive imaginary part and that the
    logarithmic gluing equations have small error.

    The shapes are coerced into the field given as argument before the
    logarithm is computed. It can be, e.g., a ComplexIntervalField.
    """

    # For each edge
    for edge in mcomplex.Edges:

        # The complex interval arithmetic value of the logarithmic
        # version of the edge equation.
        log_sum = 0

        # Iterate through edge embeddings
        for tet, perm in edge.embeddings():

            shape = _shape_for_edge_embedding(tet, perm)

            numerical_shape = NumericalField(shape)

            log_shape = numerical_shape.log()

            # Note that this is true for z in R, R < 0 as well,
            # but then it would fail for 1 - 1/z or 1 / (1-z)

            if not (log_shape.imag() > 0):
                raise ShapePositiveImaginaryPartNumericalVerifyError(
                    numerical_shape)

            # Take logarithm and accumulate
            log_sum += log_shape

        twoPiI = NumericalField.pi() * NumericalField(2j)
        
        if not abs(log_sum - twoPiI) < NumericalField(1e-7):
            raise EdgeEquationLogLiftNumericalVerifyError(log_sum)

def _shape_for_edge_embedding(tet, perm):
    """
    Given an edge embedding, find the shape assignment for it.
    If the edge embedding flips orientation, apply conjugate inverse.
    """

    # Get the shape for this edge embedding
    subsimplex = perm.image(simplex.E01)

    # Figure out the orientation of this tetrahedron
    # with respect to the edge, apply conjugate inverse
    # if differ
    if perm.sign():
        return 1 / tet.ShapeParameters[subsimplex].conjugate()
    else:
        return tet.ShapeParameters[subsimplex]
