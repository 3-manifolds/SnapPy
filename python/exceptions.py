class SnapPeaFatalError(Exception):
    """
    Exception raised by SnapPy when the SnapPea kernel encounters a fatal
    error.
    """


class InsufficientPrecisionError(Exception):
    """
    Exception raised when a computation fails and is likely to succeed if
    higher precision is used.
    """


class NonorientableManifoldError(ValueError):
    """
    Exception raised when a non-orientable manifold is given to a method
    only supporting orientable manifolds.
    """
    def __init__(self, method_name, manifold):
        self.method_name = method_name
        self.manifold = manifold

    def __str__(self):
        return ('%s only supports orientable manifolds but %s is '
                'non-orientable.') % (self.method_name, self.manifold)
