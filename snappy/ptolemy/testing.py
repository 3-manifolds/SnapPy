from __future__ import print_function
### Tests the ptolemy module
###
### Test in sage with precomputed results in testing_files_directory:
### sage -python testing.py 
###
### Test in python with precomputed results:
### python testing.py
###
### Test in sage computing the results:
### sage -python testing.py --compute
###
### Test in python computing the results using magma:
### python testing.py --compute

from snappy import Manifold, pari, ptolemy
from snappy.ptolemy import solutions_from_magma, Flattenings
from snappy.ptolemy.processMagmaFile import triangulation_from_magma
from snappy.ptolemy import __path__ as ptolemy_paths
from snappy.ptolemy.coordinates import PtolemyCannotBeCheckedException

import bz2
import sys

try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    from sage.misc.sage_eval import sage_eval
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

testing_files_directory = ptolemy_paths[0] + '/testing_files/'
testing_files_generalized_directory = ptolemy_paths[0] + '/testing_files_generalized/'

vol_tet = pari('1.014941606409653625021202554274520285941689307530299792017489106776597476258244022136470354228256695')

def check_volumes(complex_volumes, baseline_complex_volumes,
                  check_real_part_only = False,
                  torsion_imaginary_part = 6, epsilon = 1e-80):

    # add complex volumes from complex conjugation to baseline

    conjugates = [ -cvol.conj() for cvol in baseline_complex_volumes ]

    baseline_complex_volumes = baseline_complex_volumes + conjugates

    # real parts need to be equal
    # imaginary parts need to be equal up to pi^2/torstion
    p = pari('Pi * Pi') / torsion_imaginary_part

    def is_close(cvol1, cvol2):
        diff = cvol1 - cvol2

        if diff.real().abs() > epsilon:
            return False

        if check_real_part_only:
            return True

        return ((( diff.imag()) % p) < epsilon or
                ((-diff.imag()) % p) < epsilon)

    # check that every base line volume appears in computed volumes
    for cvol1 in baseline_complex_volumes:

        if not True in [
            is_close(cvol1, cvol2) for cvol2 in complex_volumes]:
            print("Missing base line volume:", cvol1)

            print("Volumes:")
            for i in complex_volumes:
                print("     ", i)

            raise Exception

    # check that every computed volume is a base line volume
    for cvol2 in complex_volumes:
        if not True in [
            is_close(cvol1, cvol2) for cvol1 in baseline_complex_volumes]:
            print("Extra volume:", cvol2)
            
            print("Volumes:")
            for i in complex_volumes:
                print("     ", i)
    
            raise Exception

def testSolutionsForManifold(M, N, solutions, baseline_cvolumes = None,
                             expect_non_zero_dimensional = None,
                             against_geometric = True):

    old_precision = pari.set_real_precision(100)

    # check solutions exactly

    found_non_zero_dimensional = False

    numerical_solutions = [ ]
    numerical_cross_ratios = [ ]

    numerical_cross_ratios_alt = [ ]

    for solution in solutions:

        if isinstance(solution,
                      ptolemy.ptolemyVariety.NonZeroDimensionalComponent):
            # encountered non-zero dimensional component

            found_non_zero_dimensional = True
        else:
            # check exact solutions
            solution.check_against_manifold(M)

            # compute numerical solutions and cross ratios
            for numerical_solution in solution.numerical():
                numerical_solutions.append(numerical_solution)
                numerical_cross_ratios.append(numerical_solution.cross_ratios())
            
            # check exact cross ratios
            cross_ratios = solution.cross_ratios()
            cross_ratios.check_against_manifold(M)
        
            # compute numerical cross ratios alternatively
            # (above we converted exact Ptolemy's to numerical and then
            # converted to cross ratio, here we first compute cross ratios
            # exactly and then convert to numerical)
            numerical_cross_ratios_alt += cross_ratios.numerical()

    # check we encountered non-zero dimensional component if that's
    # expected
    if not expect_non_zero_dimensional is None:
        assert expect_non_zero_dimensional == found_non_zero_dimensional

    # check the numerical solutions against the manifold
    for s in numerical_solutions:
        s.check_against_manifold(M, epsilon = 1e-80)

        # check that they make a flattening
        s.flattenings_numerical().check_against_manifold(M, epsilon = 1e-80)

    for s in numerical_cross_ratios:
        s.check_against_manifold(M, epsilon = 1e-80)

    for s in numerical_cross_ratios_alt:
        s.check_against_manifold(M, epsilon = 1e-80)

    # compute complex volumes and volumes
    complex_volumes = [s.complex_volume_numerical() for s in numerical_solutions]
    volumes = [s.volume_numerical() for s in numerical_cross_ratios]

    # there should be equally many
    assert len(complex_volumes) == len(volumes)

    # and real part of cvol should be equal to vol
    for vol, cvol in zip(volumes, complex_volumes):
        diff = vol - cvol.real()
        assert diff.abs() < 1e-80

    # volumes from cross ratios computed the different way
    volumes_alt = [s.volume_numerical() for s in numerical_cross_ratios_alt]

    # volumes should be equal
    assert len(volumes) == len(volumes_alt)
    volumes.sort(key = float)
    volumes_alt.sort(key = float)
    for vol1, vol2 in zip(volumes, volumes_alt):
        assert (vol1 - vol2).abs() < 1e-80

    if against_geometric:
        if M.solution_type() == 'all tetrahedra positively oriented':
            geom_vol = pari(M.volume()) * (N-1) * N * (N+1) / 6
            assert True in [
                (geom_vol - vol).abs() < 1e-12 for vol in volumes]

    # check that complex volumes match baseline volumes
    if not baseline_cvolumes is None:
        check_volumes(complex_volumes, baseline_cvolumes)

    pari.set_real_precision(old_precision)

def testComputeSolutionsForManifold(manifold, N, 
                                    compute_solutions = False, 
                                    baseline_cvolumes = None,
                                    expect_non_zero_dimensional = None):

    varieties = manifold.ptolemy_variety(N, obstruction_class = "all_original")

    if compute_solutions:
        def compute(variety):
            return variety.compute_solutions()

    else:
        def compute(variety):
            return compute_using_precomputed_magma(variety)

    solutions = sum([compute(variety) for variety in varieties], [])

    testSolutionsForManifold(manifold, N, solutions, 
                             baseline_cvolumes, expect_non_zero_dimensional)

def test_flattenings_from_tetrahedra_shapes_of_manifold():

    old_precision = pari.set_real_precision(100)

    # real parts need to be equal
    # imaginary parts need to be equal up to pi^2/6
    p = pari('Pi * Pi / 6')

    def is_close(cvol1, cvol2, epsilon):
        diff = cvol1 - cvol2
        return diff.real().abs() < epsilon and (
            ( diff.imag() % p) < epsilon or
            (-diff.imag() % p) < epsilon)

    from snappy import OrientableCuspedCensus

    for M in (list(OrientableCuspedCensus()[0:10])+
              list(OrientableCuspedCensus()[10000:10010])):
        
        flattening = Flattenings.from_tetrahedra_shapes_of_manifold(M)
        flattening.check_against_manifold(M, epsilon = 1e-80)

        is_close(flattening.complex_volume(),
                 M.complex_volume(), # returns only double precision
                 epsilon = 1e-13) 
        
    # test high precision

    M = Manifold("5_2")
    flattening = Flattenings.from_tetrahedra_shapes_of_manifold(M)
    flattening.check_against_manifold(M, epsilon = 1e-80)

    is_close(flattening.complex_volume(),
             pari('2.828122088330783162763898809276634942770981317300649477043520327258802548322471630936947017929999108 - 3.024128376509301659719951221694600993984450242270735312503300643508917708286193746506469158300472966*I'),
             epsilon = 1e-80)

    pari.set_real_precision(old_precision)

def checkSolutionsForManifoldGeneralizedObstructionClass(
    solutions_trivial, solutions_non_trivial,
    manifold, N, baseline_volumes, baseline_dimensions):

    torsionTrivial = pari('Pi^2/6 * I')
    torsionNonTrivial = pari('Pi^2/18 * I')
    
    solutions = (
        [ (s, False) for s in solutions_trivial ] +
        [ (s, True)  for s in solutions_non_trivial ])

    # Dimensions and volumes encounterd
    dimensions = set()
    volumes = []
    volumes_2 = []
    
    for solution, sol_is_non_trivial in solutions:
        # Add the dimension
        dimensions.add(solution.dimension)
        if solution.dimension == 0:

            if sol_is_non_trivial:
                got_exception = False
                try:
                    solution.check_against_manifold(manifold)
                except PtolemyCannotBeCheckedException:
                    got_exception = True
                    
                    assert got_exception, (
                        "check_against_manifold should not have passed")
            else:
                solution.check_against_manifold(manifold)
            
            fl = solution.flattenings_numerical()
            for f in fl:
                f.check_against_manifold(epsilon = 1e-80)

                cvol, modulo = f.complex_volume(with_modulo = True)

                if sol_is_non_trivial and N == 3:
                    assert (modulo - torsionNonTrivial).abs() < 1e-80, (
                        "Wrong modulo returned non-trivial case")
                else:
                    assert (modulo - torsionTrivial).abs() < 1e-80, (
                        "Wrong modulo returned trivial case")

                volumes_2.append(cvol.real())

            # Add the volumes
            volumes += solution.volume_numerical()
            
            # Check that the resulting cross ratios full fill
            # the gluing equations
            cross_ratios = solution.cross_ratios()
            cross_ratios.check_against_manifold(manifold)

    def is_close(a, b):
        return (a - b).abs() < 1e-80

    def make_unique(L):
        L.sort()
        result = L[:1]
        
        for i in L:
            if not is_close(result[-1], i):
                result.append(i)
        return result
    
    volumes = make_unique(volumes)
    volumes_2 = make_unique(volumes_2)
    all_expected_volumes = make_unique(baseline_volumes +
                                       [-vol for vol in baseline_volumes])

    assert len(all_expected_volumes) >= 2 * len(baseline_volumes) - 1

    for volume, expected_volume in zip(volumes, all_expected_volumes):
        assert is_close(volume, expected_volume)

    for volume, expected_volume in zip(volumes_2,all_expected_volumes):
        assert is_close(volume, expected_volume)

    assert dimensions == set(baseline_dimensions)

def testComputeSolutionsForManifoldGeneralizedObstructionClass(
    manifold, N, compute_solutions, baseline_volumes, baseline_dimensions):

    varieties = manifold.ptolemy_variety(N,
                                        obstruction_class = "all_generalized"
                                        #, simplify = False
                                        #, eliminate_fixed_ptolemys = True
                                        )
    
    assert len(varieties) == 2

    if compute_solutions:
        def compute(variety):
            return variety.compute_solutions()

    else:
        def compute(variety):
            return compute_using_precomputed_magma(
                variety, dir = testing_files_generalized_directory)

    # Solutions for the trivial obstruction class
    solutions_trivial = compute(varieties[0])
    
    # Solutions for the non-trivial obstruction class
    solutions_non_trivial = sum([compute(variety)
                                 for variety in varieties[1:]],
                                [])

    checkSolutionsForManifoldGeneralizedObstructionClass(
        solutions_trivial, solutions_non_trivial,
        manifold, N, baseline_volumes, baseline_dimensions)

def testGeneralizedObstructionClass(compute_solutions):

    vols = [
        pari('0'),
        2 * vol_tet
        ]
    test__m003__2 = (Manifold("m003"), # Manifold
                     2,                # N = 2
                     vols,             # expected volumes
                     [0])              # expected dimensions

    vols = [
        2 * vol_tet
        ]
    test__m004__2 = (Manifold("m004"), # Manifold
                     2,                # N = 2
                     vols,             # expected volumes
                     [0])              # expected dimensions

    vols = [
        pari('0'),
        pari('2.595387593686742138301993834077989475956329764530314161212797242812715071384508096863829303251915501'),
        2 * 4 * vol_tet,
        ]
    test__m003__3 = (Manifold("m003"), # Manifold
                     3,                # N = 3
                     vols,             # expected volumes
                     [0,1])            # expected dimensions

    test_cases = [ test__m003__2,
                   test__m004__2]

    if (not _within_sage) or (not compute_solutions):
        test_cases += [test__m003__3]

    for manifold, N, vols, dims in test_cases:

        print("Checking for", manifold.name(), "N = %d" % N)

        testComputeSolutionsForManifoldGeneralizedObstructionClass(
            manifold, N, compute_solutions, vols, dims)

def testNumericalSolutions():

    M = Manifold("m003")
    N = 3

    varities = M.ptolemy_variety(N, obstruction_class = 'all')

    solutions = [ 
        solutions_from_magma(
            get_precomputed_magma(variety,
                                  dir = testing_files_generalized_directory),
            numerical = True)
        for variety in varities ]

    for obstruction_index, obstruction in enumerate(solutions):
        for component in obstruction:
            for solution in component:
                flattenings = solution.flattenings_numerical()
                flattenings.check_against_manifold(epsilon = 1e-80)
                order = flattenings.get_order()

                if obstruction_index:
                    assert order == 6
                else:
                    assert order == 2

                cross_ratios = solution.cross_ratios()
                is_cr = cross_ratios.is_cr_structure(epsilon = 1e-80,
                                                     epsilon2 = 1e-10)

                if cross_ratios.volume_numerical().abs() < 1e-10:
                    # Every volume 0 representation of m003 happens to be
                    # cr
                    assert is_cr
                else:
                    assert not is_cr

    number_one_dimensional = 0

    allComponents = sum(solutions, [])

    dimension_dict = {}
    degree_dict = {}

    for component in allComponents:
        
        dim = component.dimension
        deg = len(component)

        dimension_dict[dim] = 1 + dimension_dict.get(dim, 0)
        degree_dict[deg] = 1 + degree_dict.get(deg, 0)

        assert (dim == 0) ^ (deg == 0)
        
    assert dimension_dict == {0: 4, 1: 2}
    assert degree_dict == {0 :2, 2: 2, 8: 2}
    
    allSolutions = sum(allComponents, [])

    allCVolumes = [s.complex_volume_numerical() for s in allSolutions]

    expected_cvolume = pari('2.595387593686742138301993834077989475956329764530314161212797242812715071384508096863829303251915501 + 0.1020524924166561605528051801006522147774827678324290664524996742369032819581086580383974219370194645*I')

    expected_cvolumes = [
        pari(0),
        expected_cvolume,
        expected_cvolume.conj(),
        2 * 4 * vol_tet
        ]

    check_volumes(allCVolumes, expected_cvolumes)

def testSageCommandLine():

    sage_eval('Manifold("m004").ptolemy_variety(3,1).compute_solutions().cross_ratios().check_against_manifold()',
              { 'Manifold' : Manifold })
    
def get_precomputed_magma(variety, dir):
    magma_file_name = ( dir + 
                        variety.filename_base() + 
                        '.magma_out.bz2')

    return bz2.BZ2File(magma_file_name,'r').read()

def compute_using_precomputed_magma(variety, dir = testing_files_directory):
    return solutions_from_magma(get_precomputed_magma(variety, dir))

def test_induced_representation():

    M = Manifold("m015")
    variety__sl2_c1 = M.ptolemy_variety(2, obstruction_class = 1)
    variety__sl3_c0 = M.ptolemy_variety(3, obstruction_class = 0)

    solutions__sl2_c1 = compute_using_precomputed_magma(
        variety__sl2_c1, dir = testing_files_generalized_directory)
    solutions__sl3_c0 = compute_using_precomputed_magma(
        variety__sl3_c0, dir = testing_files_generalized_directory)

    # Check that is_real 
    
    got_exception = False
    try:
        solutions__sl3_c0[0].cross_ratios().is_real(epsilon = 1e-80)
    except:
        got_exception = True
    assert got_exception, (
        "Expected error when calling is_real on exact solution")

    # For each component of sl3 c0, determine the number
    # of solutions which are real vs number of all solutions

    numbers_all_and_real = []
    for component in solutions__sl3_c0:
        number_real = 0
        number_all = 0
        for z in component.cross_ratios_numerical():
            if z.is_real(epsilon = 1e-80):
                number_real += 1
            number_all += 1
        numbers_all_and_real.append((number_all, number_real))
        
    # Bring into canonical form
    numbers_all_and_real.sort()

    expected_numbers = [
        (3, 1), # component has 3 solutions, 1 of them has real cross ratios
        (4, 2),
        (6, 0) ]
        
    assert numbers_all_and_real == expected_numbers, (
        "Order of components and number of real solutions is off")

    # Check is_induced_from_psl2

    number_psl2 = 0

    for component in solutions__sl3_c0:
        if component.cross_ratios().is_induced_from_psl2():
            number_psl2 += 1

    assert number_psl2 == 1, "Only one component can come from psl2"

    number_psl2 = 0

    for component in solutions__sl3_c0:
        is_induced_from_psl2 = [ z.is_induced_from_psl2(epsilon = 1e-80)
                    for z in component.cross_ratios_numerical() ]
        if True in is_induced_from_psl2:
            number_psl2 += 1
            assert not False in is_induced_from_psl2, (
                "Mixed up is_induced_from_psl2")

    assert number_psl2 == 1, (
        "Only one component can come from psl2 (numerically)")
            
    # Check that induced_representation for sl3 throws error
    got_exception = False
    try:
        solutions__sl3_c0[0].cross_ratios().induced_representation(3)
    except:
        got_exception = True
    assert got_exception, (
        "Expected error when calling induced_representation on sl3")

    # Check that induced_representation(3) works for m015

    m015_volume = pari("2.828122088330783162763898809276634942770981317300649477043520327258802548322471630936947017929999108")

    z = solutions__sl2_c1[0].cross_ratios().induced_representation(3)
    assert z.is_induced_from_psl2(), (
        "induced_representation failed to be detected as being induced")

    assert z.check_against_manifold, (
        "induced_representation fails to be valid")

    for v in z.volume_numerical():
        assert v.abs() < 1e-80 or (v.abs() - 4 * m015_volume).abs() < 1e-80, (
            "Did not get expected voluem for induced representation")

    for z in solutions__sl2_c1[0].cross_ratios_numerical():
        v = z.induced_representation(3).volume_numerical()
        assert v.abs() < 1e-80 or (v.abs() - 4 * m015_volume).abs() < 1e-80, (
            "Did not get expected voluem for induced representation")


def main():
    print("Testing in sage:", _within_sage)

    print("Running doctests...")

    import doctest
    doctest.testmod(ptolemy.component)
    doctest.testmod(ptolemy.coordinates)
    doctest.testmod(ptolemy.manifoldMethods)
    doctest.testmod(ptolemy.matrix)
    doctest.testmod(ptolemy.polynomial)
    doctest.testmod(ptolemy.processMagmaFile)
    doctest.testmod(ptolemy.ptolemyObstructionClass)
    doctest.testmod(ptolemy.ptolemyVariety)
    doctest.testmod(ptolemy.solutionsToGroebnerBasis)

    print("Testing Flattenings.from_tetrahedra_shapes_of_manifold...")

    test_flattenings_from_tetrahedra_shapes_of_manifold()

    compute_solutions = False

    if '--compute' in sys.argv:
        compute_solutions = True
        
    old_precision = pari.set_real_precision(100)

    print("Testing induced representation...")

    test_induced_representation()

    print("Running manifold tests for generalized obstruction class...")

    testGeneralizedObstructionClass(compute_solutions)

    print("Testing numerical solution retrieval method...")

    testNumericalSolutions()

    if _within_sage:
        print("Testing in sage command line...")
        testSageCommandLine()

    print("Running manifold tests...")

    ### Test for a non-hyperbolic manifold

    cvols = [ # Expected Complex volumes
        pari('0') ]
    test__3_1__2 = (Manifold("3_1"),  # Manifold
                    2,      # N = 2
                    cvols,  # expected complex volumes
                    False)  # No non-zero dimensional components

    ### Test for 4_1, amphichiral, expect zero CS

    cvols = [ # Expected Complex volumes
        2 * vol_tet
        ]
    test__4_1__2 = (Manifold("4_1"),  # Manifold
                    2,      # N = 2
                    cvols,  # expected complex volumes
                    False)   # expect non-zero dimensional components

    ### N = 3

    cvols = [ # Expected Complex volumes
        pari(0),
        2 * 4 * vol_tet
        ]
    test__4_1__3 = (Manifold("4_1"),  # Manifold
                    3,      # N = 3
                    cvols,  # expected complex volumes
                    False)   # expect non-zero dimensional components

    ### N = 4 

    cvols = [ # Expected Complex volumes
        2 * 10 * vol_tet,
        pari('8.355502146379565994303773768814775354386025704336131822255670190659042090899812770381061831431114740 - 0.7836375483069721973241186952472609466029570462055317453592429525294418326815666413123016980810744089*I'),
        pari('7.327724753417752120436828119459072886193194994253377074131984956974104158210038155834852103408920850'),
        pari('4.260549384199988638895626360783795200167603461491438308495654794999595513880915593832371886868625127 + 0.1361281655780504266700965891301086052404492485996004213292764928775628027210816440355745847599247091*I'),
        pari('3.276320849776062968038680855240250653840691843746753028019901828179320044237305241434903501877259688 + 0.03882948511714102091208888807575164800651790439786747350853616215556190252003379560451275222880494392*I'),
        pari('3.230859569867059515989136707781810556962637199799080875515780771251800351524537196576706995730457048')
        ]

    # amphicheiral, so complex volumes should appear in pairs

    cvols = cvols + [cvol.conj() for cvol in cvols]

    test__4_1__4 = (Manifold("4_1"),  # Manifold
                    4,      # N = 2
                    cvols,  # expected complex volumes
                    True)   # expect non-zero dimensional components


    ### Test for 5_2, expect non-trivial CS
    ### Number field has one real embedding with non-trival CS
    ### And one pair of complex embeddings with non-trivial CS

    cvols = [ # Expected Complex volumes
        pari('+ 1.113454552473924010022656943451126420312050780921075311799926598907813005362784871512051614741669817*I'),
        pari('2.828122088330783162763898809276634942770981317300649477043520327258802548322471630936947017929999108 - 3.024128376509301659719951221694600993984450242270735312503300643508917708286193746506469158300472966*I')
        ]
    test__5_2__2 = (Manifold("5_2"),  # Manifold
                    2,      # N = 2
                    cvols,  # expected complex volumes
                    False)   # expect no non-zero dimensional components

    ### m015 which is isometric to 5_2
    ### This example is also appearing on the website

    cvols = [ # Expected Complex volumes
        pari('+ 0.4033353624187061319128390413376061001062613536896471642214173008084760853375773296525240139108027276*I'),
        pari('- 0.4809839906489832693266177261335698864086465799360940660069682924787703897584631354526802428925968492*I'),
        pari('+ 0.6579736267392905745889660666584100756875799604827193750942232917480029881612803495334515602479034826*I'),
        pari('- 0.6579736267392905745889660666584100756875799604827193750942232917480029881612803495334515602479034826*I'),
        pari('6.332666642499251344664115407516762513592210378740462564820401485540596854320653136703448734707430999 + 0.6207993522147601522797880626542095445563442737585756367570704642807656925328117720905524433544779886*I'),
        pari('11.31248835332313265105559523710653977108392526920259790817408130903521019328988652374778807171999643 - 0.5819750380996215835728987202562276514051516606353521858642949684456185403223688691904743288635809278*I')
        ]

    test__m015__3 = (Manifold("m015"),   # Manifold
                     3,     # N = 3
                     cvols, # expected volumes
                     False)   # expect no non-zero dimensional components

    ### Test for m135
    ### Ptolemy Variety has one non-zero dimensional component

    cvols = [ # Expected Complex volumes
        pari('3.66386237670887606021841405972953644309659749712668853706599247848705207910501907791742605170446042499429769047678479831614359521330343623772637894992 + 4.93480220054467930941724549993807556765684970362039531320667468811002241120960262150088670185927611591201295688701157203888617406101502336380530883900*I')
        ]

    test__m135__2 = (Manifold("m135"),   # Manifold
                     2,           # N = 2
                     cvols,       # expected volumes
                     True)        # expect a non-zero dimensional component

    ### Test for s000
    ### Number field has one real embedding with non-trival CS
    ### And two complex embeddings with non-trivial CS

    cvols = [ # Expected Complex volumes
        pari('3.296902414326637335562593088559162089146991699269941951875989849869324250860299302482577730785960256 + 0.4908671850777469648812718224097607532197026477625119178645010117072332396428578681905186509136879130*I'),
        pari('2.270424126345803402614201085375714538056567876339790536089272218685353247084791707534069624462050547 - 0.4107712000441024931074403997367404299427214514732528671142771778838924102954454212055598588991237349*I'),
        pari('0.8061270618942289713279174295562193428686541442399728535727328789432842090598108984711933183445551654 - 0.08009598503364447177383142267302032327698119628925905075022383382334082934741244698495879201456417902*I')
        ]

    test__s000__2 = (Manifold("s000"),  # Manifold
                     2,       # N = 2
                     cvols,   # expected complex volumes
                     False)   # expect non-zero dimensional components

    ### Test for v0000
    ### This also tests the case of having more than one (here two)
    ### cusps

    cvols = [ # Expected Complex volumes
        pari('3.377597408231442496961257171798882829176020069714519460350380851901055794392493960110119513942153791 + 0.3441889979504813554136570264352067044809511501367282461849580661783699588789118851518262199077361427*I'),
        pari('1.380586178282342062845519278761524817753278777803238623999221385504661588501960339341800843520223671 + 0.06704451876157525370444676522629959282378632516497136324893149116413190872188833330658100524720353406*I')
        ]
    
    test__v0000__2 = (Manifold("v0000"),  # Manifold
                      2,       # N = 2
                      cvols,   # expected complex volumes
                      False)   # expect non-zero dimensional components


    ### Test for t00000
    ### This also tests the case of having more than one (here two)
    ### cusps

    cvols = [ # Expected Complex volumes
        pari('1.801231981344929412976326011175173926741195258454700704383597306674496002389749200442131507334117933 - 0.5490336239273931479260463113980794374780078424892352183284248305887346671477869107359071072271699046*I'),
        pari('3.434540885902563701250335473165274495380101665474171663779718033441916041515556280543077990039279103 - 0.2258335622748003900958464044389939976977579115267501160859900439335840118881246299028850924706341411*I'),
        pari('2.781430316400743786576311128801205534982324055570742403430395886754122643539081915302290955978507418 - 0.1135489194824466736456241002458946823709871812077727407261274931826622024508518397031626436088171650*I'),
        pari('0.6223401945262803853426860070209275705147832620602921139947907237662088851273537136119288228151713940 + 0.06594907226052699343130923275995552293727798462035885627276325301997714628516294342514039299674185799*I')
        ]
    
    test__t00000__2 = (Manifold("t00000"),  # Manifold
                      2,          # N = 2
                      cvols,      # expected complex volumes
                      False)      # expect non-zero dimensional components

    ### Check a big link
    ### Also an example from the website

    magma_file_name = (testing_files_directory + 
                       'DT[mcbbiceaibjklmdfgh]__sl2_c0.magma_out.bz2')
    magma_file = bz2.BZ2File(magma_file_name,'r').read()
    M = triangulation_from_magma(magma_file)

    cvols = [ # Expected complex volumes
        pari('+ 0.7322614121694386973039290771667413108310290182470824825342304051284154933661673470336385959164416421*I'),
        pari('+ 0.6559698855197901194421542799077526971594537987318704612418448259078988296037785866985260815280619566*I'),
        pari('+ 0.2247361398283578756297647616373802721002644994163825408181485469244129408530290519698535075998585820*I'),
        pari('- 0.2027422644881176424634704044404281277501368039784045126270573672316679046996639567741375438621581585*I'),
        pari('- 0.5165804080147074827389022895582362309867052040022748905555220755792330887283843134684241943604049933*I'),
        pari('- 0.6947370811282851128538476965898593866158843479399924653138487625940910866865914821776342801425993440*I'),
        pari('0.5540930202692745891960825388466975825977590423995169588735760161940190344506116205300980253997539517 + 0.4068295289298811130481298436200768517194764143946539028521976313233755795250413907401491036923026977*I'),
        pari('1.107573735539315638081527574351339252286206152514359620788778723352858972530297840033883810630372795 - 0.5969698437345460065934215233692845699554045087887375026603956452869723823599605032061849791769108861*I'),
        pari('2.930037687027569312623189945781681052278144410906236798293622026086427321638842320557907595726492995 + 0.3289208659101980133472847400704889056725461985759517831954390384444124055329869355087545294359968288*I'),
        pari('3.125178220789358723184584053682215536632004045535766511323850625359539915401419555941250425818672532 - 0.3045460696329709194663539414577718771554887499547088593863438059068830455709928110292357447645590516*I'),
        pari('4.421889017241566419690868074785716366646707154662264036268588279252613543336245268108457201245088608 + 0.6139898590804745843458673335891956572719167936270349811550161917334120397060146434004694625036091183*I'),
        pari('5.114841460302400721438106628468470523914735245629326811523719108779649554034504844499930511307337652 - 0.05922105399342951023220969900445652084832276940306885833938997476066822991673855494231393721309130829*I'),
        pari('6.222431566185758097679258799420253341029081314309961806975091338751806200485840320374430035873760649 - 0.4263260019368806638768168211894976779995931712761356294379973306148491137274105197747275282097587014*I'),
        pari('9.266495487232495996324305783134684559241087943387258046228284813849170601366096263130956023317432143 + 0.03172800820055749085258729263585301130211316529809616828980203075096987437740390637996179595789268661*I'),
        pari('11.33736073062520808632824769184286601589152410942391317127856679767752551409389094871647937428783709 - 0.5767482750608833159526763068340021498631801513808564565980389960163842898485087439099502988183940450*I'),
        pari('12.48005937399413830153052941177305767103579608241815594942757035041902652992672347868089986739774396 + 0.4828891402943609873677952178777231024869262986704386628808130740557195704279966401921665132533128189*I')
        ]

    test__morwen_link__2 = (M , # Manifold
                            2,       # N = 2
                            cvols,   # expected complex volumes
                            False)   # expect no non-zero dimensional components

    test_cases = [ test__3_1__2,
                   test__4_1__2,
                   test__5_2__2,
                   test__m135__2,
                   test__s000__2,
                   test__v0000__2,
                   test__t00000__2,
                   test__morwen_link__2,
                   test__4_1__3,
                   test__m015__3,
                   test__4_1__4]
    
    pari.set_real_precision(old_precision)

    for manifold, N, cvols, expect_non_zero_dim in test_cases:

        print("Checking for", manifold.name(), "N = %d" % N)

        testComputeSolutionsForManifold(
            manifold, N,
            compute_solutions = compute_solutions,
            baseline_cvolumes = cvols,
            expect_non_zero_dimensional = expect_non_zero_dim)

if __name__ == '__main__':
    main()
