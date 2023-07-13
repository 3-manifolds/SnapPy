class WordAppearsToBeParabolic(RuntimeError):
    def __init__(self, word, trace):
        self.word = word
        self.trace = trace
        super().__init__(
            "Attempting to drill a geodesic corresponding to a matrix "
            "that could be parabolic. "
            "Word: %s, trace: %r." % (word, trace))

class UnfinishedGraphTraceGeodesicError(RuntimeError):
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
