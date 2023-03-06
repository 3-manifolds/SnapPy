# Exceptions from the SnapPea kernel
class SnapPeaFatalError(Exception):
    """
    This exception is raised by SnapPy when the SnapPea kernel
    encounters a fatal error.
    """


class InsufficientPrecisionError(Exception):
    """
    This exception is raised when a computation fails and is likely
    to succeed if higher precision is used.
    """
