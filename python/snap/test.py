from __future__ import print_function
from ..sage_helper import _within_sage, doctest_modules
from ..pari import pari
import snappy
import snappy.snap as snap
import doctest, collections, getopt, sys

def _test_gluing_equations(manifold, shapes):
    """
    Given a manifold and exact shapes, test whether the rectangular gluing
    equations are fulfilled.
    """
    one_minus_shapes = [ 1 - shape for shape in shapes ]
    for A, B, c in manifold.gluing_equations('rect'):
        val = c
        for a, shape in zip(A, shapes):
            val *= shape ** a
        for b, one_minus_shape in zip(B, one_minus_shapes):
            val *= one_minus_shape ** b
        if not val == 1:
            return False
    return True

def test_polished(dec_prec=200):
    def test_manifold(manifold):
        eqns = manifold.gluing_equations('rect')
        shapes = manifold.tetrahedra_shapes('rect', dec_prec=dec_prec)
        return snap.shapes.gluing_equation_error(eqns, shapes)

    def test_census(name, census):
        manifolds = [M for M in census]
        print('Checking gluing equations for %d %s manifolds' % (len(manifolds), name))
        max_error = pari(0)
        for i, M in enumerate(manifolds):
            max_error = max(max_error, test_manifold(M))
            print('\r   ' + repr( (i, M) ).ljust(35) + '   Max error so far: %.2g' % float(max_error), end = '')
        print()
              
    test_census('cusped census', snappy.OrientableCuspedCensus(filter='cusps>1')[-100:])
    test_census('closed census', snappy.OrientableClosedCensus()[-100:])
    test_census('4-component links', [M for M in snappy.LinkExteriors(num_cusps=4) if M.solution_type() == 'all tetrahedra positively oriented'])

def test_holonomy(dec_prec=200):
    def test_manifold(manifold):
        # This has several internal checks which raise exceptions
        # if something is amiss
        G = snap.polished_holonomy(manifold, dec_prec=dec_prec)

    for census in [snappy.OrientableCuspedCensus, snappy.OrientableClosedCensus]:
        print('Testing holonomy of 100 manifolds in ', census)
        for manifold in census()[-100:]:
            test_manifold(manifold)

def test_fields(bits_prec=200, degree=20):
    for census in [snappy.OrientableCuspedCensus, snappy.OrientableClosedCensus]:
        print('Fields of 100 manifolds in ', census)
        for manifold in census()[:100]:
            S = snap.tetrahedra_field_gens(manifold)
            Tr = snap.trace_field_gens(manifold)
            InvTr = snap.trace_field_gens(manifold)
            hol = snap.holonomy_matrix_entries(manifold)
            for kind, X in [('shapes', S),
                            ('trace', Tr),
                            ('invtrace', InvTr),
                            ('hol', hol)]:
                K = X.find_field(bits_prec, degree)
                if K is None:
                    print('Problem with', manifold, kind)
                else:
                    if kind == 'shapes':
                        # Field is a sage number field, shapes are polynomials
                        field, numerical_root, shapes = K
                        # Turn the polynomials expressing the shapes in the
                        # root of the number field into expressions in the
                        # number field
                        shapes = [ field(shape) for shape in shapes ]
                        if not _test_gluing_equations(manifold, shapes):
                            print('Problem with', manifold,
                                  '(gluing equations violated)')

def test_ZHS(bits_prec=500, degree=20):
    for manifold in snappy.OrientableClosedCensus:
        if manifold.homology().order() == 1:
            T = snap.trace_field_gens(manifold)
            ans = T.find_field(bits_prec, degree, True)
            if ans:
                print(manifold, ans[0].polynomial())
            else:
                print(manifold, ans)

def big_test():
    test_polished()
    if _within_sage:
        test_holonomy()
        test_fields()


def run_doctests(verbose=False, print_info=True):
    from snappy.snap.t3mlite import linalg
    from snappy.snap.t3mlite import spun
    from snappy.snap.t3mlite import mcomplex
    from snappy.snap import character_varieties
    from snappy.snap import nsagetools
    from snappy.snap import polished_reps
    from snappy.snap import interval_reps
    from snappy.snap import fundamental_polyhedron

    modules = [
        mcomplex,
        linalg,
        spun,
        character_varieties,
        nsagetools,
        polished_reps,
        snap,
        interval_reps,
        fundamental_polyhedron]

    globs = {'Manifold':snappy.Manifold,
             'ManifoldHP':snappy.ManifoldHP,
             'Triangulation':snappy.Triangulation}

    return doctest_modules(modules, extraglobs=globs,
                           verbose=verbose, print_info=print_info)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    run_doctests(verbose)

