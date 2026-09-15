class VerificationError(RuntimeError):
    pass

class BadVertexGramMatrixError(VerificationError):
    '''
    Failed to verify that vertex gram matrix is realized by a finite tetrahedron.
    '''
    pass

class AngleSumIntervalNotContainingTwoPiError(VerificationError):
    '''
    The interval for the sum of the dihedral angles adjacent to an edge is
    2 * pi.
    '''
    pass

class GimbalDerivativeNotInvertibleError(VerificationError):
    pass

class VertexHasNoApproxEdgeError(VerificationError):
    pass

class NewtonMethodError(VerificationError):
    pass

class NewtonStepError(NewtonMethodError):
    pass

class NewtonMethodConvergenceError(NewtonMethodError):
    pass

class BadDihedralAngleError(VerificationError):
    pass

class KrawczykFailedToFinishError(VerificationError):
    pass

class KrawczykFailedWithBadDihedralAngleError(VerificationError):
    pass

class OrbSolutionTypeError(VerificationError):
    pass

class UnknownOrbFailureError(VerificationError):
    pass

class OrbMissingError(VerificationError):
    pass

class OrbVertexGramMatrixError(VerificationError):
    pass

class PolishingError(VerificationError):
    pass

class PolishingFailedWithBadDihedralAngleError(VerificationError):
    pass

