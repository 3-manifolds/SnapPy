from ..math_basics import correct_min, is_RealIntervalFieldElement, lower, upper

__all__ = ['unbiased_cusp_areas_from_cusp_area_matrix',
           'greedy_cusp_areas_from_cusp_area_matrix']

def unbiased_cusp_areas_from_cusp_area_matrix(cusp_area_matrix):
    """

    Examples::

        sage: from snappy.sage_helper import matrix, RIF
        sage: unbiased_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[RIF(9.0,9.0005),RIF(6.0, 6.001)],
        ...                     [RIF(6.0,6.001 ),RIF(4.0, 4.001)]]))
        [3.000?, 2.001?]

        >>> from snappy.number import Number, number_to_native_number
        >>> def N(x): return number_to_native_number(Number(x))
        >>> from snappy.matrix import make_matrix
        >>> unbiased_cusp_areas_from_cusp_area_matrix(
        ...             make_matrix([[N(10.0), N(40.0)],
        ...                          [N(40.0), N(20.0)]])) # doctest: +NUMERIC9
        [3.1622776601683795, 4.47213595499958]

    """

    verified = is_RealIntervalFieldElement(cusp_area_matrix[0, 0])
    RIF = cusp_area_matrix[0,0].parent()

    n = cusp_area_matrix.dimensions()[0]

    sqrts = [ [ cusp_area_matrix[i, j].sqrt()
                for j in range(n) ]
              for i in range(n) ]

    result = n * [ None ]

    for _ in range(n):
        lower_t_next, i = min(
            ( lower(
                sqrts[i][j]
                if result[j] is None
                else cusp_area_matrix[i, j] / result[j]),
              i)
            for i in range(n) if result[i] is None
            for j in range(n) )
        if verified:
            lower_t_next_interval = RIF(lower_t_next)
            upper_t_next = min(
                upper(
                    sqrts[i][i]
                    if j == i
                    else cusp_area_matrix[i, j] / (
                            lower_t_next_interval
                            if result[j] is None
                            else result[j]))
                for j in range(n) )
            result[i] = RIF(lower_t_next, upper_t_next)
        else:
            result[i] = lower_t_next

    return result

def greedy_cusp_areas_from_cusp_area_matrix(cusp_area_matrix, first_cusps=[]):
    """

        sage: from snappy.sage_helper import matrix, RIF
        sage: greedy_cusp_areas_from_cusp_area_matrix(
        ...             matrix([[RIF(9.0,9.0005),RIF(6.0, 6.001)],
        ...                     [RIF(6.0,6.001 ),RIF(10.0, 10.001)]]))
        [3.0001?, 2.000?]

        >>> from snappy.number import Number, number_to_native_number
        >>> def N(x): return number_to_native_number(Number(x))
        >>> from snappy.matrix import make_matrix
        >>> greedy_cusp_areas_from_cusp_area_matrix(
        ...             make_matrix([[N(10.0), N(40.0)],
        ...                          [N(40.0), N(20.0)]])) # doctest: +NUMERIC9
        [3.1622776601683795, 4.47213595499958]

    """
    num_cusps = cusp_area_matrix.dimensions()[0]

    result = list(range(num_cusps)) # Make space for range; initial values irrelevant

    # Cusp permutation given in Cayley notation
    sigma = first_cusps + [ i for i in range(num_cusps) if i not in first_cusps ]

    for i in range(num_cusps):
        stoppers = [ cusp_area_matrix[sigma[i], sigma[j]] / result[sigma[j]]
                     for j in range(i) ]
        self_stopper = cusp_area_matrix[sigma[i], sigma[i]].sqrt()

        result[sigma[i]] = correct_min(stoppers + [ self_stopper ])

    return result
