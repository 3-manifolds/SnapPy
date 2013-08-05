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

import bz2
import sys

testing_files_directory = ptolemy_paths[0] + '/testing_files/'

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

        # add complex volumes from complex conjugation to baseline

        conjugates = [ -cvol.conj() for cvol in baseline_cvolumes ]

        baseline_cvolumes = baseline_cvolumes + conjugates

        # real parts need to be equal
        # imaginary parts need to be equal up to pi^2/6
        p = pari('Pi * Pi / 6')

        def is_close(cvol1, cvol2):
            diff = cvol1 - cvol2
            return diff.real().abs() < 1e-80 and (
                ( diff.imag() % p) < 1e-80 or
                (-diff.imag() % p) < 1e-80)

        # check that every base line volume appears in computed volumes
        for cvol1 in baseline_cvolumes:
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
                is_close(cvol1, cvol2) for cvol1 in baseline_cvolumes]:
                print("Extra complex volume:", cvol2)

                print("Volumes:")
                for i in complex_volumes:
                    print("     ", i)

                raise Exception

    pari.set_real_precision(old_precision)

def testComputeSolutionsForManifold(manifold, N, 
                                    compute_solutions = False, 
                                    baseline_cvolumes = None,
                                    expect_non_zero_dimensional = None):

    varities = manifold.ptolemy_variety(N, obstruction_class = "all")

    if compute_solutions:
        def compute(variety):
            return variety.compute_solutions()

    else:
        def compute(variety):
            magma_file_name = ( testing_files_directory + 
                                variety.filename_base() + 
                                '.magma_out.bz2')

            magma_file = bz2.BZ2File(magma_file_name,'r').read()
            
            return solutions_from_magma(magma_file)

    solutions = sum([compute(variety) for variety in varities], [])

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

def main():
    print("Running doctests...")

    import doctest
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

    print("Running manifold tests...")

    compute_solutions = False

    if '--compute' in sys.argv:
        compute_solutions = True
        
    old_precision = pari.set_real_precision(100)

    ### Test for a non-hyperbolic manifold

    cvols = [ # Expected Complex volumes
        pari('0') ]
    test__3_1__2 = (Manifold("3_1"),  # Manifold
                    2,      # N = 2
                    cvols,  # expected complex volumes
                    False)  # No non-zero dimensional components

    ### Test for 4_1, amphichiral, expect zero CS

    vol_tet = pari('1.014941606409653625021202554274520285941689307530299792017489106776597476258244022136470354228256695')

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
