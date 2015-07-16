def compute_geometric_solution(M, N = 2, numerical = False,
                               engine = None,
                               memory_limit = 750000000, directory = None,
                               prefer_rur = False, data_url = False,
                               verbose = None):

    """
    Given a manifold M, compute the exact or numerical solutions to the
    Ptolemy variety (pass numerical = True for numerical solutions). The additional
    options are the same as the ones passed to compute_solutions of a PtolemyVariety.

    >>> from snappy.ptolemy.geometricRep import compute_geometric_solutions
    >>> compute_geometric_solutions(Manifold("m004")) #doctest: +SKIP

    """

    if verbose is None:
        verbose = (engine == 'retrieve')

    obstruction_class = 'all' if N % 2 == 0 else None
    varities = M.ptolemy_variety(N, obstruction_class = obstruction_class)

    if engine == 'retrieve':
        sols = varities.retrieve_solutions(numerical = numerical,
                                           prefer_rur = prefer_rur,
                                           data_url = data_url,
                                           verbose = verbose)
    else:
        sols = varities.compute_solutions(engine = engine,
                                          numerical = numerical,
                                          memory_limit = memory_limit,
                                          directory = directory,
                                          verbose = verbose)

    for sol in sols.flatten(3):
        if sol.is_geometric():
            return sol

    return None

def retrieve_geometric_solution(M, N = 2,
                                numerical = False,
                                prefer_rur = False,
                                data_url = None,
                                verbose = True):

    """
    Given a manifold M, retrieve the exact or numerical solutions to the
    Ptolemy variety (pass numerical = True for numerical solutions). The additional
    options are the same as the ones passed to retrieve_solutions of a PtolemyVariety.

    >>> from snappy.ptolemy.geometricRep import retrieve_geometric_solutions
    >>> retrieve_geometric_solutions(Manifold("m004")) #doctest: +SKIP

    """

    
    return compute_geometric_solution(M, N, numerical = numerical,
                                      engine = 'retrieve',
                                      prefer_rur = prefer_rur,
                                      data_url = data_url, verbose = verbose)
