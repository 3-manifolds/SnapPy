from ..sage_helper import _within_sage
from .mathHelpers import interval_aware_min

__all__ = ['unbiased_cusp_areas_from_cusp_area_matrix',
           'greedy_cusp_areas_from_cusp_area_matrix']

if _within_sage:
    from sage.rings.real_mpfi import is_RealIntervalFieldElement

    # python's sqrt only work for floats
    # They would fail or convert to float losing precision
    from sage.all import sqrt
else:
    import math

    # Otherwise, define our own sqrt which checks whether
    # the given type defines a sqrt method and fallsback
    # to python's log and sqrt which has the above drawback of
    # potentially losing precision.
    def sqrt(x):
        if hasattr(x, 'sqrt'):
            return x.sqrt()
        return math.sqrt(x)


def unbiased_cusp_areas_from_cusp_area_matrix(cusp_area_matrix):
    """

    Examples::

        sage: from sage.all import matrix, RIF
        sage: unbiased_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[RIF(9.0,9.0005),RIF(6.0, 6.001)],
        ...                     [RIF(6.0,6.001 ),RIF(4.0, 4.001)]]))
        [3.00?, 2.000?]

        >>> from snappy.SnapPy import matrix
        >>> unbiased_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[10.0, 40.0],
        ...                     [40.0, 20.0]]))
        [3.1622776601683795, 4.47213595499958]

    """

    if _within_sage:
        if is_RealIntervalFieldElement(cusp_area_matrix[0, 0]):
            return _verified_unbiased_cusp_areas_from_cusp_area_matrix(
                cusp_area_matrix)

    return _unverified_unbiased_cusp_areas_from_cusp_area_matrix(
                                                cusp_area_matrix)

def greedy_cusp_areas_from_cusp_area_matrix(cusp_area_matrix, first_cusps=[]):
    
    """

        sage: from sage.all import matrix, RIF
        sage: greedy_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[RIF(9.0,9.0005),RIF(6.0, 6.001)],
        ...                     [RIF(6.0,6.001 ),RIF(10.0, 10.001)]]))
        [3.0001?, 2.000?]

        >>> from snappy.SnapPy import matrix
        >>> greedy_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[10.0, 40.0],
        ...                     [40.0, 20.0]]))
        [3.1622776601683795, 4.47213595499958]
    
    """

    num_cusps = cusp_area_matrix.dimensions()[0]

    result = list(range(num_cusps)) # Make space for range; initial values irrelevant

    # Cusp permutation given in Cayley notation
    sigma = first_cusps + [ i for i in range(num_cusps) if i not in first_cusps ]
    
    for i in range(num_cusps):
        stoppers = [ cusp_area_matrix[sigma[i], sigma[j]] / result[sigma[j]]
                     for j in range(i) ]
        self_stopper = sqrt(cusp_area_matrix[sigma[i], sigma[i]])

        result[sigma[i]] = interval_aware_min(stoppers + [ self_stopper ])

    return result

def _interval_minimum_candidates(intervals_and_extras):
    
    result = [ intervals_and_extras[0] ]
    for this in intervals_and_extras[1:]:
        if not all([this[0] > other[0] for other in result]):
            if all([this[0] < other[0] for other in result]):
                result = [ this ]
            else:
                result.append(this)
            
    return result

def _find_potential_stoppers(cusp_area_matrix, assigned_areas):
    def stopper(i, j):
        if not assigned_areas[i] is None:
            return cusp_area_matrix[i,j] / assigned_areas[i]
        if not assigned_areas[j] is None:
            return cusp_area_matrix[i,j] / assigned_areas[j]
        return sqrt(cusp_area_matrix[i, j])

    num_cusps = cusp_area_matrix.dimensions()[0]

    return [ (stopper(i, j), (i, j))
             for i in range(num_cusps)
             for j in range(i, num_cusps)
             if (assigned_areas[j] is None) or (assigned_areas[i] is None) ]

def _find_stoppers(cusp_area_matrix, assigned_areas):
    return _interval_minimum_candidates(
        _find_potential_stoppers(cusp_area_matrix, assigned_areas))

def _union_intervals(intervals):
    result = intervals[0]
    for i in intervals[1:]:
        result = result.union(i)
    return result

def _get_cusps_from_stoppers(stoppers, assigned_areas):
    result = set()
    for stopper in stoppers:
        for cusp in stopper[1]:
            if assigned_areas[cusp] is None:
                result.add(cusp)
    return result

def _verified_unbiased_cusp_areas_from_cusp_area_matrix(
                                                cusp_area_matrix):

    num_cusps = cusp_area_matrix.dimensions()[0]

    result = num_cusps * [ None ]
    
    while None in result:
        stoppers = _find_stoppers(cusp_area_matrix, result)

        stoppers_union = _union_intervals([ stopper[0] for stopper in stoppers ])
        cusps = _get_cusps_from_stoppers(stoppers, result)
        stopper_pairs = set([stopper[1] for stopper in stoppers])

        stop_size = (stoppers_union * stoppers_union) / stoppers_union

        for area in result:
            if not area < stop_size:
                raise Exception("New area smaller than existing areas")

        for cusp in cusps:
            result[cusp] = stop_size

        for i in range(num_cusps):
            for j in range(i, num_cusps):
                if i in cusps or j in cusps:
                    if (not result[i] is None) and (not result[j] is None):
                        if not (i,j) in stopper_pairs:
                            if not result[i] * result[j] < cusp_area_matrix[i,j]:
                                raise Exception("Violated maximal cusp area", i, j)

    return result

def _find_minimal_stopper(cusp_area_matrix, assigned_areas):
    return min(_find_potential_stoppers(cusp_area_matrix, assigned_areas))

def _unverified_unbiased_cusp_areas_from_cusp_area_matrix(
                                                cusp_area_matrix):
    num_cusps = cusp_area_matrix.dimensions()[0]
    num_pending = num_cusps

    result = num_cusps * [ None ]

    while num_pending > 0:
        stop_size, (i, j) = _find_minimal_stopper(cusp_area_matrix, result)

        if result[i] is None:
            result[i] = stop_size
            num_pending -= 1
        if i != j:
            if result[j] is None:
                result[j] = stop_size
                num_pending -= 1

    return result

