class DrillGeodesicError(RuntimeError):
    pass

class WordAppearsToBeParabolic(DrillGeodesicError):
    def __init__(self, word, trace):
        self.word = word
        self.trace = trace
        super().__init__(
            "Attempting to drill a geodesic corresponding to a matrix "
            "that could be parabolic. "
            "Word: %s, trace: %r." % (word, trace))

class GeodesicSystemNotSimpleError(DrillGeodesicError):
    def __init__(self, maximal_tube_radius):
        self.maximal_tube_radius = maximal_tube_radius
        super().__init__(
            "One of the given geodesics might not simple or two of the "
            "given geodesics might intersect. "
            "The maximal tube radius about the given system of geodesics "
            "was estimated to be: %r." % maximal_tube_radius)


class UnfinishedGraphTraceGeodesicError(DrillGeodesicError):
    def __init__(self, steps):
        self.steps = steps
        super().__init__(
            "The line fixed by the given word could not be moved (by "
            "Decktransformations) to intersect the fundamental domain after "
            "%d steps. "
            "This is probably due to a pathology, e.g., a bad conjugate "
            "was picked and the line is very far away from the fundamental "
            "domain or the given geodesic is very close to a core curve of "
            "a filled cusp." % steps)


class UnfinishedTraceGeodesicError(DrillGeodesicError):
    def __init__(self, steps):
        self.steps = steps
        super().__init__(
            "The geodesic seems to have more than %d pieces in the "
            "triangulation. This is probably due to a pathology, "
            "e.g., the geodesic is very close to a core curve of "
            "filled cusp." % steps)

class GeodesicHittingOneSkeletonError(DrillGeodesicError):
    """
    Base class for exceptions caused by the geodesic hitting the
    1-skeleton and that can be avoided by perturbing the geodesic.
    """


class GeodesicStartPointOnTwoSkeletonError(GeodesicHittingOneSkeletonError):
    """
    Raised when the start point given to GeodesicInfo appears not to be in the
    interior of a tetrahedron.
    """


class RayHittingOneSkeletonError(GeodesicHittingOneSkeletonError):
    """
    Raised when the geodesic appears to intersect the 1-skeleton of the
    original triangulation.
    """


class RetracingRayHittingOneSkeletonError(GeodesicHittingOneSkeletonError):
    """
    Raised when the geodesic appears to intersect the 1-skeleton of the
    subdivided triangulation.
    """
