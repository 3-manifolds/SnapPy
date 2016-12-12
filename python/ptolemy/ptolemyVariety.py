from __future__ import print_function
from . import matrix
from . import homology
from .polynomial import Polynomial
from .ptolemyObstructionClass import PtolemyObstructionClass
from .ptolemyGeneralizedObstructionClass import PtolemyGeneralizedObstructionClass
from .ptolemyVarietyPrimeIdealGroebnerBasis import PtolemyVarietyPrimeIdealGroebnerBasis
from . import processFileBase, processFileDispatch, processMagmaFile
from . import utilities
from string import Template
import signal
import re
import os

try:
    from sage.rings.rational_field import RationalField 
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.symbolic.ring import var as sageVariable
    from sage.rings.ideal import Ideal
    _within_sage = True
except ImportError:
    _within_sage = False

try:
    from urllib import urlopen
    from urllib import quote as urlquote
except ImportError: # Python 3
    from urllib.request import urlopen
    from urllib.request import quote as urlquote

class PtolemyFileMissingError(Exception):
    """
    An exception indicating that a requested solution file was missing.
    """

    def __init__(self, message):
        Exception.__init__(self, message)

class PtolemyVariety(object):
    """
    Holds a reduced Ptolemy variety.

    === Examples ===
    
    To generate such a variety, call:

    >>> from snappy import Manifold
    >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

    Show the equations and variables:

    >>> for e in p.equations: print(e)
    - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2
    c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
    - 1 + c_0011_0
    >>> p.variables
    ['c_0011_0', 'c_0101_0']

    Show as an ideal (sage object):

    >>> p.ideal    #doctest: +SKIP
    Ideal (-c_0011_0^2 + c_0011_0*c_0101_0 + c_0101_0^2, -c_0011_0^2 - c_0011_0*c_0101_0 + c_0101_0^2, c_0011_0 - 1) of Multivariate Polynomial Ring in c_0011_0, c_0101_0 over Rational Field                                                       
    (skip doctest because example only works in sage and not plain python)


    Produce magma input:

    >>> s = p.to_magma()
    >>> print(s.split('ring and ideal')[1].strip())    #doctest: +ELLIPSIS
    R<c_0011_0, c_0101_0> := PolynomialRing(RationalField(), 2, "grevlex");
    MyIdeal := ideal<R |
              - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2,
        ...

    Call ``p.compute_solutions()`` to automatically compute solutions!

    Show canonical representatives:

    (The Ptolemy coordinates c_0110_0 and c_0101_0 are identified, this 
    information is needed to recover all Ptolemy coordinates from the solutions
    of a simplified Ptolemy variety. The information is also packaged into a
    python section by py_eval_variable_dict().)

    >>> p.canonical_representative["c_0110_0"]
    (1, 0, 'c_0101_0')

    """

    
    def __init__(self, manifold, N, obstruction_class,
                 simplify, eliminate_fixed_ptolemys):
        self._manifold = manifold
        self._N = N
        self._obstruction_class = obstruction_class

        if obstruction_class:
            obstruction_class._checkManifoldAndN(manifold, N)

        if isinstance(obstruction_class, PtolemyObstructionClass):
            self._identified_variables_from_obstruction = (
                obstruction_class.identified_variables)
        else:
            self._identified_variables_from_obstruction = [ ]

        H2_class = None

        if isinstance(obstruction_class, PtolemyGeneralizedObstructionClass):
            H2_class = obstruction_class.H2_class

        self._identified_coordinates = (
            manifold._ptolemy_equations_identified_coordinates(N, H2_class))

        self._action_by_decoration_change = (
            manifold._ptolemy_equations_action_by_decoration_change(N))

        # find enough Ptolemy variables to set to one so that the
        # decoration is fixed
        self._fixed_ptolemy_coordinates = (
            _fix_decoration(N, self._action_by_decoration_change))
        
        self._identified_variables = (
            self._identified_coordinates +
            self._identified_variables_from_obstruction)

        self._ptolemy_relations = (
            _generate_ptolemy_relations(
                N, manifold.num_tetrahedra(),
                isinstance(obstruction_class, PtolemyObstructionClass)))

        self.equations = [eqn for eqn in self._ptolemy_relations]

        order_of_u = 1

        if isinstance(obstruction_class, PtolemyGeneralizedObstructionClass):
            order_of_u, equations = (
                obstruction_class._get_equation_for_u(N))
            self.equations += equations

        if eliminate_fixed_ptolemys:

            # each ptolemy set to 1 for fixing decoration is eliminated by
            # being replace by 1
            self._identified_variables += (
                [ (+1, 0, ptolemy_coord, 1)
                  for ptolemy_coord in self._fixed_ptolemy_coordinates])

        else:
            one = Polynomial.constant_polynomial(1)

            # we add an equation c_XXXX_X - 1 for enough ptolemy's
            # to fix the decoration
            self.equations += (
                [ Polynomial.from_variable_name(ptolemy_coord) - one
                  for ptolemy_coord in self._fixed_ptolemy_coordinates])
                  
        variables = _union([eqn.variables() for eqn in self.equations])

        if simplify:

            self.canonical_representative = _identified_variables_canonize(
                self._identified_variables)

            substitution = (
                _canonical_representative_to_polynomial_substituition(
                    self.canonical_representative, order_of_u))

            self.equations = [eqn.substitute(substitution)
                              for eqn in self.equations]

        else:

            self.canonical_representative = { }
            
            for sign, power, var1, var2 in self._identified_variables:

                self.canonical_representative[var1] = (+1, 0, var1)
                if var2 != 1:
                    self.canonical_representative[var2] = (+1, 0, var2)

                if order_of_u == 2:
                    u = Polynomial.constant_polynomial(-1)
                else:
                    u = Polynomial.from_variable_name('u')


                firstTerm = (
                    Polynomial.from_variable_name(var1) *
                    u ** (power % order_of_u))

                if var2 == 1:
                    secondTerm = (
                        Polynomial.constant_polynomial(sign))
                else:
                    secondTerm = (
                        Polynomial.constant_polynomial(sign) *
                        Polynomial.from_variable_name(var2))
                self.equations.append(firstTerm - secondTerm)

        self.variables = _union([ eqn.variables() 
                            for eqn in self.equations])
                
        # Process interior Ptolemy coordinates such as c_1111_x
        # Only invoked for N >= 4
        for var in self.variables:
            if var[0:2] == 'c_':
                if not var in self.canonical_representative:
                    self.canonical_representative[var] = (+1, 0, var)

        self.variables_with_non_zero_condition = [ "t" ] + self.variables

        # take out u, the root of unity
        vars_without_u = [ var
                           for var in self.variables_with_non_zero_condition
                           if not var == 'u']

        self._non_zero_condition = (
            _non_zero_condition(vars_without_u))

        self.equations_with_non_zero_condition = (
            self.equations + [ self._non_zero_condition ])

        if _within_sage:
            def sage_monomial(monomial):
                r = monomial.get_coefficient()
                for varName, expo in monomial.get_vars():
                    r = r * (sageVariable(varName) ** expo)
                return r

            def sage_eqn(eqn):
                return sum([sage_monomial(m) for m in eqn.get_monomials()])

            def sage_ideal(vars, eqns):
                
                polynomialRing = PolynomialRing(
                    RationalField(), vars, order = 'lex')

                return Ideal(
                    polynomialRing, [ sage_eqn(eqn) for eqn in eqns ])

            self.ideal = sage_ideal(
                self.variables,
                self.equations)

            self.ideal_with_non_zero_condition = sage_ideal(
                self.variables_with_non_zero_condition,
                 self.equations_with_non_zero_condition)
            
    def py_eval_variable_dict(self):

        def create_dict_entry(var1, val):
            sign, power, var2 = val

            assert sign in [+1, -1]
            
            p = ""
            if self._N == 2:
                sign *= (-1) ** power
            else:
                if power % self._N:
                    p = " * d['u'] ** %d" % (power % self._N)

            if sign == +1:
                return "'%s' : d['%s']%s" % (var1, var2, p)
            else:
                return "'%s' : - d['%s']%s" % (var1, var2, p)

        format_str = "(lambda d: {\n          %s})"

        return (
            format_str % ',\n          '.join(
                [create_dict_entry(key, val) 
                 for key, val 
                 in list(self.canonical_representative.items())
                 if not key == 1]))

    def py_eval_section(self):
        """
        Returns a string that can be evaluated in python and contains extra
        information needed to recover solutions from a simplified Ptolemy
        variety.

        >>> from snappy import Manifold, pari
        >>> M = Manifold('4_1')
        >>> p = M.ptolemy_variety(2, obstruction_class = 1)

        Get extra information

        >>> eval_section = p.py_eval_section()
        >>> print(eval_section)    #doctest: +ELLIPSIS
        {'variable_dict' : 
             (lambda d: {
            ...

        Turn it into a python object by evaluation.

        >>> obj = eval(eval_section)

        Access the function that expands a solution to the simplified
        Ptolemy variety to a full solution.

        >>> variable_dict = obj['variable_dict']

        Setup a solution and expand it to a full solution, '1' must map to 1

        >>> simplified_solution = {'c_0101_0' : pari('0.5 - 0.866025403784439*I'), '1' : pari(1), 'c_0011_0' : pari(1)}
        >>> full_solution = variable_dict(simplified_solution)

        Full solution is a dictionary with a key for every Ptolemy coordinate
        
        >>> full_solution['c_1010_1']
        1
        >>> for tet in range(2):
        ...     for i in utilities.quadruples_with_fixed_sum_iterator(2, skipVertices = True):
        ...         c = "c_%d%d%d%d" % i + "_%d" % tet
        ...         assert c in full_solution
        """

        result = "{"
        result += "'variable_dict' : \n     %s" % self.py_eval_variable_dict()

        # If we have a non-trivial generalized obstruction class,
        # add an extra key to the dictionary to mark it.

        # This will prevent PtolemyCoordinates to compute the ill-defined
        # flattenings and complex volume.

        if isinstance(self._obstruction_class,
                      PtolemyGeneralizedObstructionClass):
            if self._obstruction_class._is_non_trivial(self._N):
                result += (
                    ",\n "
                    "'non_trivial_generalized_obstruction_class' : True")

        result += "}"
    
        return result
                
    def to_magma_file(
            self, filename,
            template_path = "magma/default.magma_template"):
        
        """
        >>> from snappy import *
        >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

        >>> p.to_magma_file('/tmp/tmp_magma_file.magma')
        """
        open(filename,'w').write(self.to_magma(template_path = template_path))

    def to_magma(
            self,
            template_path = "magma/default.magma_template"):

        """
        Returns a string with the ideal that can be used as input for magma.

        The advanced user can provide a template string to write own magma
        code to process the equations.
        
        >>> from snappy import *
        >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

        Magma input to compute radical Decomposition
        >>> s = p.to_magma()
        >>> print(s.split('ring and ideal')[1].strip())    #doctest: +ELLIPSIS  +NORMALIZE_WHITESPACE
        R<c_0011_0, c_0101_0> := PolynomialRing(RationalField(), 2, "grevlex");
        MyIdeal := ideal<R | - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2,
        ...
        >>> "RadicalDecomposition" in p.to_magma()
        True
        """

        if os.path.isfile(template_path):
            template = open(template_path, 'r').read()
        else:
            from snappy.ptolemy import __path__ as base_paths
            abs_path = os.path.join(base_paths[0], template_path)
            if os.path.isfile(abs_path):
                template = open(abs_path, 'r').read()
            else:
                raise Exception("No file at template_path %s" % template_path)
            
        PREAMBEL = (
            "==TRIANGULATION=BEGINS==\n" +
            self._manifold._to_string() + "\n"
            "==TRIANGULATION=ENDS==\n" +
            "PY=EVAL=SECTION=BEGINS=HERE\n" +
            self.py_eval_section() + "\n"
            "PY=EVAL=SECTION=ENDS=HERE\n")

        # magma will wrap long lines and we potentially loose some information
        # that way when retrieving the ASCII text encoding the triangulation
        # and the py_eval dictionary.
        # To prevent this, we already break long lines here in such a way that
        # join_long_lines gives exactly back the same ASCII text.

        # Next, we quote the text so that we can give it to magma's print.

        QUOTED_PREAMBEL = utilities.quote_ascii_text(
            utilities.break_long_lines(PREAMBEL))

        return Template(template).safe_substitute(
            PREAMBEL = PREAMBEL,
            QUOTED_PREAMBEL = QUOTED_PREAMBEL,

            VARIABLES = (
                ", ".join(self.variables)),
            VARIABLES_QUOTED = (
                ", ".join(['"%s"' % v for v in self.variables])),
            VARIABLE_NUMBER = (
                len(self.variables)),

            VARIABLES_WITH_NON_ZERO_CONDITION = (
                ", ".join(self.variables_with_non_zero_condition)),
            VARIABLES_WITH_NON_ZERO_CONDITION_QUOTED = (
                ", ".join(['"%s"' % v
                           for v in self.variables_with_non_zero_condition])),
            VARIABLE_WITH_NON_ZERO_CONDITION_NUMBER = (
                len(self.variables_with_non_zero_condition)),

            EQUATIONS = (
                ',\n          '.join(
                    [str(eqn)
                     for eqn in self.equations])),
            EQUATIONS_WITH_NON_ZERO_CONDITION = (
                ',\n          '.join(
                    [str(eqn)
                     for eqn in self.equations_with_non_zero_condition])))
        
    def filename_base(self):
        """
        Preferred filename base for writing out this Ptolemy variety

        >>> from snappy import *
        >>> M = Manifold('4_1')
        >>> p = M.ptolemy_variety(2, obstruction_class = 1)
        >>> p.filename_base()
        '4_1__sl2_c1'

        >>> p = M.ptolemy_variety(2)
        >>> p.filename_base()
        '4_1__sl2_c0'
        """

        obstruction_class = "NoIndex"

        if self._obstruction_class is None:
            obstruction_class = "0"
        elif not self._obstruction_class._index is None:
            obstruction_class = "%d" % self._obstruction_class._index
            
        return '%s__sl%d_c%s' % (self._manifold.name(), self._N,
                                 obstruction_class)

    def path_to_file(self):

        name = self._manifold.name()

        if re.match('([msvt]|o9_)[0-9]+$', name):
            dir = 'OrientableCuspedCensus'
        elif re.match('[0-9]+([\^][0-9]+)?[_][0-9]+$', name):
            dir = 'LinkExteriors'
        elif re.match('[KL][0-9]+[an][0-9]+$', name):
            dir = 'HTLinkExteriors'
        else:
            raise Exception('No canonical path for manifold')

        tets = self._manifold.num_tetrahedra()

        return '/'.join(['data', 
                         'pgl%d' % self._N,
                         dir,
                         '%02d_tetrahedra' % tets])

    def _solution_file_url(self, data_url = None, rur = False):

        if data_url is None:
            from . import DATA_URL as data_url

        if not '://' in data_url:
            # No schema in url, assume file
            if not data_url[0] == '/':
                data_url = '/' + data_url
            data_url = 'file://' + data_url
 
        # Make it end in /
        if not data_url[-1] == '/':
            data_url = data_url + '/'

        if rur:
            ext = '.rur'
        else:
            ext = '.magma_out'

        filename = self.filename_base() + ext

        pathological_1dim = ["L14n24426__sl2_c3.magma_out"]

        if filename in pathological_1dim:
            filename = "truncated_" + filename

        return data_url + self.path_to_file() + '/' + urlquote(filename)

    def _retrieve_solution_file(self, data_url = None, prefer_rur = False,
                                verbose = False):

        # First try to retrieve solutions from the URL corresponding to
        # the prefered format (i.e., RUR vs mamga decomposition)

        url = self._solution_file_url(data_url = data_url,
                                      rur = prefer_rur)
        if verbose:
            print("Trying to retrieve solutions from %s ..." % url)

        try:
            return _retrieve_url(url)

        except PtolemyFileMissingError:

            # If that file wasn't there, try to retrieve solutions from URL
            # corresponding to the non-prefered format

            url = self._solution_file_url(data_url = data_url,
                                          rur = not prefer_rur)
            if verbose:
                print("Retrieving solutions instead from %s ...:" % url)
            return _retrieve_url(url)
            

    def retrieve_decomposition(self, data_url = None, verbose = True):
        
        url = self._solution_file_url(data_url = data_url, rur = False)
        if verbose:
            print("Retrieving decomposition from %s ..." % url)

        text = _retrieve_url(url)
        
        if verbose:
            print("Parsing...")
            
        M = processFileBase.get_manifold(text)
        assert M._to_bytes() == self._manifold._to_bytes(), (
            "Manifold does not match census manifold")

        return processMagmaFile.decomposition_from_magma(text)

    def retrieve_solutions(self, numerical = False,
                           prefer_rur = False,
                           data_url = None,
                           verbose = True):

        text = self._retrieve_solution_file(data_url = data_url,
                                            prefer_rur = prefer_rur,
                                            verbose = verbose)
        if verbose:
            print("Parsing...")

        M = processFileBase.get_manifold(text)
        assert M._to_bytes() == self._manifold._to_bytes(), (
            "Manifold does not match census manifold")

        return processFileDispatch.parse_solutions(text, numerical = numerical)

    def __repr__(self):
        
        res =  "Ptolemy Variety for %s, N = %d" % (self._manifold.name(), 
                                                   self._N)
        if not self._obstruction_class is None:
            res += ", obstruction_class = "
            if not self._obstruction_class._index is None:
                res += "%d" % self._obstruction_class._index
                if isinstance(self._obstruction_class,
                              PtolemyGeneralizedObstructionClass):
                    res += " (generalized)"
            elif isinstance(self._obstruction_class,
                            PtolemyGeneralizedObstructionClass):
                res += "%s" % self._obstruction_class.H2_class
            else:
                res += "..."

        res += "\n" + "\n".join(['    %s' % e for e in self.equations])

        return res

    def compute_decomposition(
        self,
        engine = None,
        memory_limit = 750000000,
        directory = None,
        verbose = False,
        template_path = "magma/default.magma_template"):

        """
        Starts an engine such as magma to compute the
        radical decomposition of the Ptolemy variety.

        If started in sage, uses sage, otherwise, uses magma.

        === Arguments ===

        engine --- engine to use, currently, only support magma and sage
        memory_limit --- maximal allowed memory in bytes
        directory --- location for input and output files, temporary directory used if not specified
        verbose --- print extra information
        """

        if engine is None:
            if _within_sage:
                engine = 'sage'
            else:
                engine = 'magma'

        if engine == 'magma':
            return processMagmaFile.run_magma(
                self.to_magma(template_path = template_path),
                filename_base = self.filename_base(),
                memory_limit = memory_limit,
                directory = directory,
                verbose = verbose)
            
        if engine == 'sage':

            M = self._manifold.copy()

            radical = self.ideal_with_non_zero_condition.radical()
            
            sage_radical_decomp = radical.primary_decomposition()

            def process_component(component):
                
                dimension = component.dimension()

                if dimension == 0:
                    sage_gb = component.groebner_basis()
                    polys = [ Polynomial.parse_string(str(p)) for p in sage_gb ]
                else:
                    polys = []

                return PtolemyVarietyPrimeIdealGroebnerBasis(
                    polys = polys,
                    term_order = 'lex',
                    size = None,
                    dimension = dimension,
                    is_prime = component.is_prime(),
                    free_variables = None,
                    py_eval = eval(self.py_eval_section()),
                    manifold_thunk = lambda :M)
                    
            return utilities.MethodMappingList(
                [ process_component(component)
                  for component in sage_radical_decomp 
                  if not component.is_one()])



    def compute_solutions(self,
                          engine = None,
                          numerical = False,
                          memory_limit = 750000000,
                          directory = None,
                          verbose = False):

        """
        Starts an engine such as magma to compute the
        radical decomposition of the ideal and computes exact solutions.

        If started in sage, uses sage, otherwise, uses magma.

        === Arguments ===

        engine --- engine to use, currently, only support magma and sage
        numerical --- get numerical solutions from magma instead of exact ones
        memory_limit --- maximal allowed memory in bytes
        directory --- location for input and output files, temporary directory used if not specified
        verbose --- print extra information
        """

        decomposition = self.compute_decomposition(
            engine = engine,
            memory_limit = memory_limit,
            directory = directory,
            verbose = verbose)

        
        return utilities.MethodMappingList(
                [ component.solutions(numerical = numerical)
                  for component in decomposition ])

    def degree_to_shapes(self):
        """
        In general, there can be d different solutions to the (reduced) Ptolemy
        variety for each solution to the gluing equations (with peripheral
        equations). This method computes d which is also the number of elements
        in H^1(Mhat; Z/N) where Mhat is the non-manifold cell comples obtained
        by gluing together the tetrahedra as non-ideal tetrahedra.


        For example, the Ptolemy variety for m009 has 4 points but there are
        only 2 distinct boundary-unipotent PSL(2,C) representations.
        Thus the following call returns 2 = 4 / 2

        >>> from snappy import Manifold
        >>> Manifold("m009").ptolemy_variety(2,1).degree_to_shapes()
        2

        >>> Manifold("m010").ptolemy_variety(2,1).degree_to_shapes()
        2
        >>> Manifold("m011").ptolemy_variety(2,1).degree_to_shapes()
        1

        >>> Manifold("m009").ptolemy_variety(3,1).degree_to_shapes()
        1
        >>> Manifold("m010").ptolemy_variety(3,1).degree_to_shapes()
        3
        >>> Manifold("m011").ptolemy_variety(3,1).degree_to_shapes()
        1
        
        """

        # Boundary maps for chain complex
        d2 = self._manifold._ptolemy_equations_boundary_map_2()[0]
        d1 = self._manifold._ptolemy_equations_boundary_map_1()[0]
        
        # Boundary maps for dual chain complex
        co_d1 = matrix.matrix_transpose(d2)
        co_d0 = matrix.matrix_transpose(d1)

        # All cohomology classes
        cohomology_classes = homology.homology_representatives(
            co_d1, co_d0, self._N)

        # Number of classes
        return len(list(cohomology_classes))

def _fix_decoration(N, action_by_decoration_change):
        
    action_matrix, ptolemy_coords, decorations_to_be_fixed = (
        action_by_decoration_change)

    return matrix.get_independent_rows(
        action_matrix, ptolemy_coords, desired_determinant = N)

def _generate_ptolemy_relations(N, num_tet,
                                has_obstruction_class):

    def generate_ptolemy_relation(tet, index):

        def generate_Ptolemy_coordinate(addl_index):
            total_index = matrix.vector_add(index, addl_index)
            return Polynomial.from_variable_name(
                "c_%d%d%d%d" % tuple(total_index) + "_%d" % tet)

        def generate_obstruction_variable(face):
            if has_obstruction_class:
                return Polynomial.from_variable_name(
                    "s_%d_%d" % (face, tet))
            else:
                return Polynomial.constant_polynomial(1)

        # implements equation 5.8 from paper
        
        return (
            generate_obstruction_variable(0) *
            generate_obstruction_variable(1) *
            generate_Ptolemy_coordinate((1,1,0,0)) *
            generate_Ptolemy_coordinate((0,0,1,1)) 
          - generate_obstruction_variable(0) *
            generate_obstruction_variable(2) *
            generate_Ptolemy_coordinate((1,0,1,0)) *
            generate_Ptolemy_coordinate((0,1,0,1))            
          + generate_obstruction_variable(0) *
            generate_obstruction_variable(3) *
            generate_Ptolemy_coordinate((1,0,0,1)) *
            generate_Ptolemy_coordinate((0,1,1,0)))

    return [generate_ptolemy_relation(tet, index)
            for tet in range(num_tet)
            for index in utilities.quadruples_with_fixed_sum_iterator(N-2)]

def _non_zero_condition(variables):
    one = Polynomial.constant_polynomial(1)
    
    polynomial = one
    
    for var in variables:
        polynomial = polynomial * Polynomial.from_variable_name(var)
        
    polynomial = polynomial - one

    return polynomial

def _union(lists):
    all = sum(lists, [])
    all = list(set(all))
    all.sort()
    return all

def _identified_variables_canonize(identified_variables):

    def merge_two_dicts(sign, power, var1, var2, dict1, dict2):

        sign1, power1 = dict1[var1]
        sign2, power2 = dict2[var2]

        new_sign  = sign1  * sign  * sign2
        new_power = power1 - power - power2 

        for v2, (s2, p2) in dict2.items():
            dict1[v2] = (s2 * new_sign, p2 + new_power)
        
        return dict1

    all_variables = { }

    for sign, power, var1, var2 in identified_variables:
        all_variables[var1] = { var1 : (+1, 0) }
        all_variables[var2] = { var2 : (+1, 0) }

    for sign, power, var1, var2 in identified_variables:
        if not all_variables[var1] is all_variables[var2]:
            new_dict = merge_two_dicts(sign, power, var1, var2,
                                       all_variables[var1],
                                       all_variables[var2])
            for var in new_dict.keys():
                all_variables[var] = new_dict
                
    result = { }

    for variable, variable_dict in all_variables.items():
        if variable not in result:
            vars = list(variable_dict.keys())
            if 1 in vars:
                canonical_rep = 1
            else:
                vars.sort()
                canonical_rep = vars[0]

            canonical_rep_sign, canonical_rep_power = (
                variable_dict[canonical_rep])

            for (var, (sign, power)) in variable_dict.items():
                result[var] = (canonical_rep_sign * sign,
                               canonical_rep_power - power,
                               canonical_rep)
    
    return result

def _canonical_representative_to_polynomial_substituition(
        canonical_representative, order_of_u):

    result = { }

    for var1, signed_var2 in canonical_representative.items():
        sign, power, var2 = signed_var2
        if not var1 == var2:

            if order_of_u == 2:
                u = Polynomial.constant_polynomial(-1)
            else:
                u = Polynomial.from_variable_name('u')

            sign_and_power = (
                Polynomial.constant_polynomial(sign) *
                 u ** (power % order_of_u))

            if var2 == 1:
                result[var1] =  sign_and_power
            else:
                result[var1] = (sign_and_power *
                                Polynomial.from_variable_name(var2))

    return result

def _retrieve_url(url):

    overview_url = "http://ptolemy.unhyperbolic.org/data/overview.html"
            
    try:
        # Remember SnapPy's SIGALRM handler (defined in app.py)
        # And temporarily disable it (except under windows which does not
        # have SIGALRM)
        sigalrm_handler = None
        if hasattr(signal, 'SIGALRM'):
            sigalrm_handler = signal.signal(signal.SIGALRM, signal.SIG_IGN)
        s = urlopen(url)
    except IOError as e:
        # IOError: this means the file wasn't there or we couldn't connect
        # to the server

        if url[:5] == 'http:': 
            # IOError for http means we could not connect to server
            raise RuntimeError(
                "Problem connecting to server while retrieving %s: "
                "%s" % (url, e))
        else:
            # Otherwise, it probably means the file wasn't there
            raise PtolemyFileMissingError(
                "The ptolemy variety has probably not been computed "
                "yet or the given data_url or environment variable "
                "PTOLEMY_DATA_URL "
                "is not configured correctly: %s. Also see %s" % (
                    e, overview_url))
    finally:
        # Always restore the original signal handler
        if sigalrm_handler: # Not supported under windows
            signal.signal(signal.SIGALRM, sigalrm_handler)
            
    # Read the text
    text = s.read().decode('ascii')
        
    if url[:5] != 'http:':
        # If this is a normal file, we are done
        return text

    # For http, we need to check that this isn't an error message
    # by checking the HTTP code
    code = s.getcode()
    if code == 200:
        # HTTP code 200 means alright
        return text

    # Otherwise we have an error
    httpErr = "(HTTP Error %d while accessing %s)" % (code, url)

    if code == 404:
        # 404 means file not found
        raise PtolemyFileMissingError(
            "The ptolemy variety has probably not been computed "
            "yet, see %s (%s)" % (overview_url, httpErr))

    # Everything else means something is wrong with the server
    # configuration or the connection to the server
    raise RuntimeError(
        "Problem retrieving file from server, please report to "
        "enischte@gmail.com: %s" % httpErr)


