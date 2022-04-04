class DrillGeodesicError(RuntimeError):
    pass

class WordAppearsToBeParabolic(DrillGeodesicError):
    def __init__(self, word, trace):
        self.word = word
        self.trace = trace
        super().__init__(
            "Attempting to drill a geodesic corresponding to parabolic "
            "matrix. Word: %s, trace: %r." % (word, trace))

class UnfinishedGraphTraceGeodesicError(DrillGeodesicError):
    def __init__(self, steps):
        self.steps = steps
        super().__init__(
            "??? after %d steps. This is probably due to a pathology, "
            "e.g., big conjugate" % steps)

class UnfinishedTraceGeodesicError(DrillGeodesicError):
    def __init__(self, steps):
        self.steps = steps
        super().__init__(
            "The geodesic seems to have more than %d pieces in the "
            "triangulation. This is probably due to a pathology, "
            "e.g., the geodesic is very close to a core curve of "
            "filled cusp." % steps)

class GeodesicStartPointOnTwoSkeletonError(DrillGeodesicError):
    pass

class RayHittingOneSkeletonError(DrillGeodesicError):
    pass

class RetracingRayHittingOneSkeletonError(DrillGeodesicError):
    pass

class GeodesicNotSimpleError(DrillGeodesicError):
    pass

class GeodesicCloseToCoreCurve(DrillGeodesicError):
    pass
