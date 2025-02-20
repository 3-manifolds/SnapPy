"""This module provides the classes NewtonShapeCertifer and KrawczyShapeCertifer.
These classes are instantiated with a Manifold and a precision specification,
given in binary or decimal digits with the options bit_prec or dec_prec.  They
produce a set of "certified" tetrahedra shapes.  This means that it can be
proven that there exists a solution to the Manifold's gluing equations which
is contained in an interval vector centered at the certified shapes and having
radius 2^-p + 2^-pj where p is the specified binary precision.

The certification involves first using Newton's method to compute a set Z0 of
tetrahedron shapes for M to a "working precision" higher than the precision
being certified.  The working precision may be provided explicitly with the
working_prec option.  If not provided it will be determined automatically.

We then consider the interval vector Z centered at Z0 with radius
2^-p + 2^-pj. We verify that Z contains a solution to the gluing equations.
Actually, we certify that Z contains a solution to the logarithmic gluing
equations, obtained by taking the log of both sides of the gluing equations.
This is because the logarithmic derivative of the rational function on the
left hand side of a gluing equation is a linear expression in log(z_i) and
log(z_i - 1) with small integer coefficients, making it very convenient to
work with.

The criteria used to show that Z contains a solution are given
by [XXX LOOKUP PUBLISHED REFERENCES].

"""

from ..snap.shapes import polished_tetrahedra_shapes, enough_gluing_equations
from ..number import Number, bit_precision
from flint import arb, acb, acb_mat
import math

__all__ = ['NewtonShapeCertifier', 'KrawczykShapeCertifier']

dec_to_bits = math.log(10) / math.log(2)

# Can it be true that flint does not provide an identity acb_mat?
def acb_identity(n:int) -> acb_mat:
    ID = acb_mat(n, n, 0)
    for n in range(n):
        ID[n, n] = 1
    return ID

class ShapeCertifierBase:
    """Base class for Newton and Krawczyk Shape Certifiers.

    Defines methods for computing the quantities used by both of the Certifier
    classes.

    """

    def __init__(self, M, bits_prec:int=None, dec_prec:int=None,
                 working_prec:int=None)-> None:
        if not self.test_interval:
            raise ValueError(
                'The ShapeCertifierBase class cannot be instantiated '
                'because it has no test_interval method.  Subclasses '
                'should define that method.')

        # Convert decimal precision to binary precision if needed.
        if dec_prec:
            self.precision = round(dec_to_bits * dec_prec)
        elif bits_prec:
            self.precision = bits_prec
        else:
            raise ValueError("One of dec_prec or bits_prec must not be None.")

        # If the working precision is not provided do something reasonable.
        if working_prec:
            self.high_precision = working_prec
        else:
            self.high_precision = self.precision + 30
        
        # Verify that the manifold is orientable.
        if not M.is_orientable():
            raise ValueError("Manifold needs to be orientable")

        # Compute the interval Z.  The component intervals of Z will be high
        # precision acbs with center at one of the precise shapes and radius set
        # to 2^-p where p is the bit precision being certified.
        
        self.manifold = M
        self.equations = enough_gluing_equations(M)
        self.dim = len(self.equations)
        precise_shapes = M.tetrahedra_shapes('rect',
                            bits_prec=self.high_precision)
        with bit_precision(self.high_precision):
            epsilon = arb(2.0) ** -self.precision
            self.Z0 = acb_mat(self.dim, 1,
                              [z.flint_obj.mid() for z in precise_shapes])
            self.Z = acb_mat(self.dim, 1,
                             [acb(arb(z.real.mid(), epsilon),
                                  arb(z.imag.mid(), epsilon))
                              for z in self.Z0])
        # Shapes have not been certified yet.
        self.certified_shapes = None

    def fZ0(self)-> acb_mat:
        """Compute and return f(z0).

        The result is a column of intervals containing the log of the
        left hand side of each gluing equation.

        """
        logs = []
        with bit_precision(self.high_precision):
            for a, b, c in self.equations:
                prod = arb(c)
                for i, z in enumerate(self.Z0):
                    # Use int exponents so multiplication will be used.
                    prod *=  (z ** int(a[i])) * ((1 - z) ** int(b[i]))
                logs.append(prod.log())
            return acb_mat(self.dim, 1, logs)

    def Df_Z0(self)-> acb_mat:
        """Compute and return Df_Z0.

        The result is a matrix of intervals containing the derivative matrix
        of f at Z0.

        """
        with bit_precision(self.high_precision):
            return acb_mat([[eqn[0][i] / z - eqn[1][i] / (1 - z)
                 for i, z in enumerate(self.Z0)] for eqn in self.equations])

    def Df_Z(self)-> acb_mat:
        """Compute and return the interval [Df_Z].

        [Df_Z] is a matrix of intervals containing the derivative matrix Df_z
        for every point z of the interval Z.
        """
        with bit_precision(self.high_precision):
            return acb_mat([[a[i] / z - b[i] / (1 - z)
                             for i, z in enumerate(self.Z)]
                             for a, b, _ in self.equations])

    def certify(self):
        # Construct the interval to be used for the test.
        T = self.test_interval()
        # Verify that the test interval is contained in Z
        with bit_precision(self.precision):
            is_contained = [z.contains(w) for z, w in zip(self.Z, T)]
            result = all(is_contained)
            if result:
                self.certified_shapes = [
                    Number(z, precision=self.high_precision)
                    for z in self.Z]
                for z in self.certified_shapes:
                    z._certified = True
        return result

class KrawczykShapeCertifier(ShapeCertifierBase):
    """
    >>> M = Manifold('v1234')
    >>> from snappy.verify.shape_certifiers import KrawczykShapeCertifier
    # Instantiate a shape certifier with a manifold and a precision.
    # An internal working precision will be selected automatically.
    >>> K = KrawczykShapeCertifier(M, bits_prec=53)
    >>> K.certify()
    True
    # The certify method finds and stores the certified shapes.
    >>> z = K.certified_shapes[1]; z
    1.06629329123300... + 0.24909676637578...j
    # The certified interval for each shape is available.
    >>> z.interval()
'[10662932912330067784e-19, 10662932912330070008e-19] + [24909676637578189449e-20, 24909676637578211657e-20]j'
    # The internal working precision can be specified explicitly.
    # If it is too small, certification will fail.
    >>> K = KrawczykShapeCertifier(M, bits_prec=53, working_prec=64)
    >>> K.certify()
    False
    >>> K = KrawczykShapeCertifier(M, bits_prec=53, working_prec=70)
    >>> K.certify()
    True

    """

    def krawczyk_interval(self):
        """
        Computeand return the interval K for the Krawczyk test.

        K is defined as

            K(Z0, [Z], f) := Z0 - c * f(Z0) + (Id - c * Df_Z)) * (Z - Z0)

        where
        
           - f evaluates the lhs of the logarithmic gluing equations;
           - [Z] is a column of intervals of shapes;
           - Z0 is the column of midpoints of the intervals in Z;
           - c is an approximate inverse of Df_Z0;
           - Df_Z is the interval matrix containing the derivative of f at
             each point of Z.

        Theorem XXX of [XXX] says that f has a root in Z if K ⊂ Z.
        
        """

        with bit_precision(self.high_precision):
            ID = acb_identity(self.dim)
            Df_Z = self.Df_Z()
            Df_Z0 = self.Df_Z0()
            fZ0 = self.fZ0()
            c = Df_Z0.inv()
            delta = self.Z - self.Z0
            return self.Z0 - c * fZ0 + (ID - c * Df_Z) * delta

    test_interval = krawczyk_interval

class NewtonShapeCertifier(ShapeCertifierBase):
    """
    >>> M = Manifold('v1234')
    >>> from snappy.verify.shape_certifiers import NewtonShapeCertifier
    # Instantiate a shape certifier with a manifold and a precision.
    # An internal working precision will be selected automatically.
    >>> N = NewtonShapeCertifier(M, bits_prec=53)
    # The certify method finds and stores the certified shapes.
    >>> N.certify()
    True
    >>> z = N.certified_shapes[1]; z
    1.06629329123300... + 0.24909676637578...j
    # The certified interval for each shape is available.
    >>> z.interval()
    '[10662932912330067784e-19, 10662932912330070008e-19] + [24909676637578189449e-20, 24909676637578211657e-20]j'
    # The internal working precision can be specified explicitly.
    # If it is too small, certification will fail.
    >>> N = NewtonShapeCertifier(M, bits_prec=53, working_prec=60)
    >>> N.certify()
    False
    >>> N = NewtonShapeCertifier(M, bits_prec=53, working_prec=70)
    >>> N.certify()
    True

    """

    def newton_interval(self)-> acb_mat: 
        """ Computes the interval N for the Newton test.

        N is defined as:
        
        N(Z0, Z, f) := Z0 - [Df_Z]^-1 * f(Z0)

        where
        
           - f evaluates the lhs of the logarithmic gluing equations;
           - [Z] is a column of intervals of shapes;
           - Z0 is the column of midpoints of the intervals in Z;
           - Df_Z is the interval matrix containing the derivative of f at
             each point of Z.

        Rather than compute the inverse of Df_Z and multiply by f(z0)
        we compute [Df_Z]^-1 * f(Z0) by solving for an interval vector
        [X] in the equation:

        [Df_Z] * [X]  =  f(Z0)
         
        Theorem XXX of [XXX] says that f has a root in Z if N ⊂ Z.

        """
        with bit_precision(self.high_precision):
            fZ0 = self.fZ0()
            Df_Z = self.Df_Z()
            Df_inv_fZ0 = Df_Z.solve(fZ0)
            image = self.Z0 - Df_inv_fZ0
            return image

    test_interval = newton_interval
