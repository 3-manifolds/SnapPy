from ..math_basics import is_RealIntervalFieldElement # type: ignore
from ..exceptions import InsufficientPrecisionError # type: ignore

from typing import Sequence

__all__ = ['floor_as_intergers', 'SpatialDict']


def floor_as_integers(x) -> Sequence[int]:
    """
    Computes floor of a number or interval, returning a list of integers
    if evaluating floor is ambiguous.

        sage: floor_as_integers(RIF(1.8, 1.9))
        [1]
        sage: floor_as_integers(RIF(1.9, 2.1))
        [1, 2]

        >>> floor_as_integers(1.4)
        [1]
        >>> floor_as_integers(2.01)
        [1, 2]
        >>> floor_as_integers(1.99)
        [1, 2]

    """

    if is_RealIntervalFieldElement(x):
        f = x.floor()
        l = f.lower().round()
        u = f.upper().round() + 1

        n = u - l
        if n > 5:
            raise InsufficientPrecisionError(
                "Too many integers (%d) in given interval to compute floor. "
                "Increasing precision should fix it." % n)

        return list(range(l, u))
    else:
        f = x.floor()
        int_f = int(f)
        d = x - f
        if d < 0.125:
            return [ int_f - 1, int_f ]
        if d > 0.875:
            return [ int_f, int_f + 1 ]
        return [ int_f ]


class _Entry:
    """
    A helper for SpatialDict.

    The implementation of SpatialDict has the same instance of _Entry
    stored for multiple keys so that updating the value for all keys
    can be done by assigning the new value to _Entry.value only once.
    """
    def __init__(self, value):
        self.value = value


class SpatialDict:
    """
    A python dict-like object appropriate for using numerical points (e.g.,
    in the hyperboloid model) as keys. That is, look-ups return
    the same entry for points that are almost but not exactly the
    same due to rounding-errors.

    To achieve this, the points are assumed to be in some lattice
    and the minimal distance between any two points in the lattice
    must be given.
    """

    _scale = 1024

    def __init__(self, min_distance, verified):
        RF = min_distance.parent()

        self._min_distance = min_distance
        self._RF_scale = RF(self._scale)

        if verified:
            self._right_distance_value = min_distance
            self._left_distance_value = 0
        else:
            self._right_distance_value = min_distance * RF(0.125)
            self._left_distance_value = min_distance * RF(0.5)

        self._data = { }

    def setdefault(self, point, default):
        reps_and_ikeys = self._representatives_and_ikeys(point)

        for rep, ikey in reps_and_ikeys:
            for other_rep, entry in self._data.get(ikey, []):
                d = self.distance(rep, other_rep)
                if d < self._right_distance_value:
                    return entry.value
                if not (self._left_distance_value < d):
                    raise InsufficientPrecisionError(
                        "Could neither verify that the two given tiles are "
                        "the same nor that they are distinct. "
                        "Distance between basepoint translates is: %r. "
                        "Injectivty diameter about basepoint is: %r." % (
                            d, self._min_distance))

        entry = _Entry(default)
        for rep, ikey in reps_and_ikeys:
            self._data.setdefault(ikey, []).append((rep, entry))

        return default

    def distance(self, point_0, point_1):
        raise NotImplementedError()

    def representatives(self, point):
        # Applies, e.g., translation by geodesic matrix
        return [ point ]

    def float_hash(self, point):
        raise NotImplementedError()

    def _representatives_and_ikeys(self, point):
        return [
            (rep, ikey)
            for rep in self.representatives(point)
            for ikey in floor_as_integers(self._RF_scale * self.float_hash(rep)) ]
