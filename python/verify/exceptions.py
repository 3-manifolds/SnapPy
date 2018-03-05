"""
All final excepions are deriving from two base classes:

- a subclass of VerifyErrorBase to indicate whether a numerical or exact
  verification failed
- a subclass of EquationType to indicate the type of equation of
  inequality for which the verification failed.

Intermediate subclasses (those without __init__) are not supposed to be
raised.

The hierarchy is as follows:

- VerifyErrorBase(RuntimeError)

  - NumericalVerifyError

    - InequalityNumericalVerifyError
    - LogLiftNumericalVerifyError

  - ExactVerifyError

    - IsZeroExactVerifyError

- EquationType

  - EdgeEquationType

    - EdgeEquationExactVerifyError
    - EdgeEquationLogLiftNumericalVerifyError
  - CuspConsistencyType

    - CuspEquationType

      - CuspEquationExactVerifyError
      - CuspEquationLogLiftNumericalVerifyError
    - CuspDevelopmentType

      - CuspDevelopmentTypeExactVerifyError
  - TiltType

    - TiltInequalityNumericalVerifyError

      - TiltProvenPositiveNumericalVerifyError
    - TiltIsZeroExactVerifyError
  - ShapeType

    - ShapePositiveImaginaryPartNumericalVerifyError
  - ConsistencyWithSnapPeaType

    - ConsistencyWithSnapPeaNumericalVerifyError
"""

class VerifyErrorBase(RuntimeError):
    """
    The base for all exceptions related to verification.
    """

class NumericalVerifyError(VerifyErrorBase):
    """
    The base for all exceptions resulting from a failed numerical
    verification of an equality (using some epsilon) or inequality
    (typically by interval arithmetics).
    """

class InequalityNumericalVerifyError(NumericalVerifyError):
    """
    The base for all exceptions resulting from a failed numerical
    verification of an inequality (typically by interval arithmetics).
    """

class LogLiftNumericalVerifyError(NumericalVerifyError):
    """
    To verify a logarithmic gluing equation, the verify module will usually
    first verify the corresponding polynomial gluing equation.
    This means that the logarithmic gluing equation will be fulfilled up
    to a multiple of 2 Pi I.
    It then computes the logarithms and numerically checks that the result
    is close (by some epsilon) to the right value. Because we already know
    that the difference is a multiple of 2 Pi I, checking closeness is enough.

    This exception is supposed to be raised if the polynomial gluing equations
    have passed but checking the logarithmic equation is epsilon-close has
    failed.
    """

class ExactVerifyError(VerifyErrorBase):
    """
    The base for all exceptions resulting from a failed verification of an
    equation using exact arithmetics.
    """

class IsZeroExactVerifyError(VerifyErrorBase):
    """
    The base for all exceptions resulting from verifying that a desired
    quantity is zero using exact arithmetics.
    """

class EquationType(object):
    """
    A base class to derive subclasses which indicate what kind of
    equation failed to be verified.
    """

class EdgeEquationType(EquationType):
    """
    A base class indicating that an edge equation could not be verified.
    """

class EdgeEquationExactVerifyError(ExactVerifyError,
                                   EdgeEquationType):
    """
    Exception for failed verification of a polynomial edge equation
    using exact arithmetics.
    """

    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return ('Verification of a polynomial edge equation using exact '
                'arithmetic failed: %r == 1' % self.value)

class EdgeEquationLogLiftNumericalVerifyError(LogLiftNumericalVerifyError,
                                              EdgeEquationType):
    """
    Exception for failed numerical verification that a logarithmic edge
    equation has error bound by epsilon.
    """
    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return ('Numerical verification that logarthmic edge equation has '
                'small error failed: %r == 2 Pi I' % self.value)

class CuspConsistencyType(EquationType):
    """
    A base class indicating that verificatin of an equation involving a cusp
    failed.
    """

class CuspEquationType(CuspConsistencyType):
    """
    A base class indicating that a cusp gluing equation (involving the
    shapes) failed.
    """

class CuspEquationExactVerifyError(ExactVerifyError,
                                   CuspEquationType):
    """
    Exception for failed verification of a polynomial cusp gluing equation
    using exact arithmetics.
    """

    def __init__(self, value, expected_value):
        self.value = value
        
    def __str__(self):
        return ('Verification of a polynomial cusp equation using exact '
                'arithmetic failed: %r == 1' % self.value)

class CuspEquationLogLiftNumericalVerifyError(LogLiftNumericalVerifyError,
                                              CuspEquationType):
    """
    Exception for failed numerical verification that a logarithmic cusp
    equation has error bound by epsilon.
    """
    def __init__(self, value, expected_value):
        self.value = value
        self.expected_value = expected_value
        
    def __str__(self):
        return ('Numerical verification that logarthmic cusp equation has '
                'small error failed: '
                '%r == %r' % (self.value, self.expected_value))


class CuspDevelopmentType(CuspConsistencyType):
    """
    A base class indicating that there was a failure to find a consistent
    assignment of side lengths to the Euclidean Horotriangles to form a
    Euclidean Horotorus for a cusp.
    """

class CuspDevelopmentExactVerifyError(ExactVerifyError,
                                      CuspDevelopmentType):
    """
    Raised when finding a consistent assignment of side lengths to the
    Euclidean Horotriangles to form a Euclidean Horotorus for a cusp failed
    using exact arithmetic.
    """
    
    def __init__(self, value1, value2):
        self.value1 = value1
        self.value2 = value2

    def __str__(self):
        return ('Inconsistency in the side lengths of the Euclidean '
                'Horotriangles for a cusp: '
                '%r = %r' % (self.value1, self.value2))

class TiltType(EquationType):
    """
    A base class relating to tilts.
    """

class TiltInequalityNumericalVerifyError(InequalityNumericalVerifyError,
                                         TiltType):
    """
    Numerically verifying that a tilt is negative has failed.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return ('Numerical verification that tilt is negative has '
                'failed: %r < 0' % self.value)


class TiltProvenPositiveNumericalVerifyError(
                                   TiltInequalityNumericalVerifyError):
    """
    Numerically verifying that a tilt is negative has not only failed, we
    proved that the tilt is positive and thus that this cannot be a
    proto-canonical triangulation.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return ('Numerical verification that tilt is negative has '
                'failed, tilt is actually positive. This is provably '
                'not the proto-canonical triangulation: %r <= 0' % self.value)

class TiltIsZeroExactVerifyError(IsZeroExactVerifyError,
                                 TiltType):
    """
    Verifying that a tilt is zero has failed using exact arithmetic.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return ('Verification that tilt is zero has failed using exact '
                'arithmetic: %r == 0' % self.value)

class ShapeType(EquationType):
    """
    Base class for failed verification of legal shapes.
    """

class ShapePositiveImaginaryPartNumericalVerifyError(
    InequalityNumericalVerifyError,
    ShapeType):
    """
    Failed numerical verification of a shape having positive imaginary part.
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return ('Numerical verification that shape has positive imaginary '
                'part has failed: Im(%r) > 0' % self.value)

class ConsistencyWithSnapPeaType(EquationType):
    """
    A base class for exceptions raised when there is a difference
    between the values computed by the SnapPea kernel and by this module
    for a given quantity.
    """

class ConsistencyWithSnapPeaNumericalVerifyError(
    NumericalVerifyError,
    ConsistencyWithSnapPeaType):
    """
    Exception raised when there is a significant numerical difference
    between the values computed by the SnapPea kernel and by this module
    for a given quantity.
    """
    def __init__(self, value, snappea_value):
        self.value = value
        self.snappea_value = snappea_value

    def __str__(self):
        return ('Inconsistency between SnapPea kernel and verify: '
                '%r == %r' % (self.snappea_value, self.value))
