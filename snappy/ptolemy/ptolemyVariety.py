import matrix
from polynomial import Polynomial
import coordinates
import solutionsToGroebnerBasis

try:
    from sage.rings.rational_field import RationalField 
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    from sage.symbolic.ring import var as sageVariable
    from sage.rings.ideal import Ideal
    _within_sage = True
except ImportError:
    _within_sage = False

class PtolemyVariety(object):
    """
    Represents a Ptolemy variety as described in
    Garoufalidis, Thurston, Zickert
    The Complex Volume of SL(n,C)-Representations of 3-Manifolds
    http://arxiv.org/abs/1111.2828
    
    Garoufalidis, Goerner, Zickert:
    Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
    http://arxiv.org/abs/1207.6711

    This variety can be used to find pSL(n,C) representations and associated
    invariants such as complex volume.

    === Examples ===
    
    To generate such a variety, call:

    >>> from snappy import Manifold
    >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

    Show the equations and variables:

    >>> for e in p.equations: print e
    1 - c_0101_0 + c_0101_0^2
    - 1 + c_0101_0 - c_0101_0^2
    >>> p.variables
    ['c_0101_0']

    Show as an ideal (sage object):

    >>> p.ideal    #doctest: +SKIP
    Ideal (c_0101_0^2 - c_0101_0 + 1, -c_0101_0^2 + c_0101_0 - 1) of Univariate Polynomial Ring in c_0101_0 over Rational Field
    (skip doctest because example is different in sage and plain python)


    Produce a magma file:

    >>> print p.to_magma()     #doctest: +ELLIPSIS
    P<t, c_0101_0> := PolynomialRing(RationalField(), 2);
    I := ideal<P |
    1 - c_0101_0 + c_0101_0^2,
        ...

    Call p.compute_solutions() to automatically compute solutions!

    Show canonical representatives:

    (The Ptolemy coordinates c_0110_0 and c_0101_0 are identified, this 
    information is needed to recover all Ptolemy coordinates from the solutions
    of a simplified Ptolemy variety. The information is also packaged into a
    python section by py_eval_variable_dict().)

    >>> p.canonical_representative["c_0110_0"]
    (1, 'c_0101_0')
    """

    
    def __init__(self, manifold, N, obstruction_class = None, simplify = True):
        self._manifold = manifold
        self._N = N
        self._obstruction_class = obstruction_class

        if obstruction_class:
            assert obstruction_class._manifold == manifold, (
                "PtolemyObstructionClass for wrong manifold")
            assert N % 2 == 0, (
                "PtolemyObstructionClass only makes sense for even N")
            self._identified_variables_from_obstruction = (
                obstruction_class.identified_variables)
        else:
            self._identified_variables_from_obstruction = [ ]

        self._identified_coordinates = (
            manifold._ptolemy_equations_identified_coordinates(N))

        self._action_by_decoration_change = (
            manifold._ptolemy_equations_action_by_decoration_change(N))

        self._identified_coordinates_fixing_decoration = (
            _fix_decoration(self._action_by_decoration_change))

        self._identified_variables = (
            self._identified_coordinates +
            self._identified_coordinates_fixing_decoration +
            self._identified_variables_from_obstruction)

        self._ptolemy_relations = (
            _generate_ptolemy_relations(N, manifold.num_tetrahedra(),
                                        obstruction_class))

        self.equations = [eqn for eqn in self._ptolemy_relations]

        variables = _union([eqn.variables() for eqn in self.equations])

        if simplify:

            self.canonical_representative = _identified_variables_canonize(
                self._identified_variables)

            substitution = (
                _canonical_representative_to_polynomial_substituition(
                    self.canonical_representative))

            self.equations = [eqn.substitute(substitution)
                              for eqn in self.equations]

        else:

            self.canonical_representative = { }
            
            for sign, var1, var2 in self._identified_variables:
                firstTerm = Polynomial.from_variable_name(var1)
                if var2 == 1:
                    secondTerm = (
                        Polynomial.constant_polynomial(sign))
                else:
                    secondTerm = (
                        Polynomial.constant_polynomial(sign) *
                        Polynomial.from_variable_name(var2))
                self.equations.append(firstTerm - secondTerm)

        for var in variables:
            if not self.canonical_representative.has_key(var):
                self.canonical_representative[var] = (+1, var)

        self.variables = _union([ eqn.variables() 
                                  for eqn in self.equations])

        self.variables_with_non_zero_condition = [ "t" ] + self.variables

        self._non_zero_condition = (
            _non_zero_condition(self.variables_with_non_zero_condition))

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
            sign, var2 = val

            assert sign in [+1, -1]
            
            if sign == +1:
                return "'%s' : d['%s']" % (var1, var2)
            else:
                return "'%s' : negation(d['%s'])" % (var1, var2)

        format_str = "(lambda d, negation = (lambda x:-x): {\n          %s})"

        return (
            format_str % ',\n          '.join(
                [create_dict_entry(key, val) 
                 for key, val 
                 in self.canonical_representative.items()
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
        >>> print eval_section    #doctest: +ELLIPSIS
        {'variable_dict' : 
             (lambda d, negation = (lambda x:-x): {
                  's_3_1' : d['1'],
                  's_3_0' : negation(d['1']),
            ...

        Turn it into a python object by evaluation.

        >>> obj = eval(eval_section)

        Access the function that expands a solution to the simplified
        Ptolemy variety to a full solution.

        >>> variable_dict = obj['variable_dict']

        Setup a solution and expand it to a full solution, '1' must map to 1

        >>> simplified_solution = {'c_0101_0' : pari('0.5 - 0.866025403784439*I'), '1' : pari(1)}
        >>> full_solution = variable_dict(simplified_solution)

        Full solution is a dictionary with a key for every Ptolemy coordinate
        
        >>> full_solution['c_1010_1']
        1
        >>> for tet in range(2):
        ...     for i in list_all_quadruples_with_fixed_sum(2, True):
        ...         assert full_solution.has_key("c_%d%d%d%d" % tuple(i) + "_%d" % tet)
        """
    
        return ("{'variable_dict' : \n     %s}" % 
                self.py_eval_variable_dict())
                
    def to_magma_file(self, filename, primary_decomposition = True):
        """
        >>> from snappy import *
        >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

        >>> p.to_magma_file('/tmp/tmp_magma_file.magma')
        """
        open(filename,'w').write(self.to_magma(primary_decomposition))

    def to_magma(self, primary_decomposition = True):

        """
        Returns a string with the ideal that can be used as input for magma.
        It will contain magma code to compute the Groebner basis or
        primary decomposition.
        
        >>> from snappy import *
        >>> p = Manifold("4_1").ptolemy_variety(2, obstruction_class = 1)

        Magma file to compute Primary Decomposition
        >>> print p.to_magma()          #doctest: +ELLIPSIS
        P<t, c_0101_0> := PolynomialRing(RationalField(), 2);
        I := ideal<P |
        1 - c_0101_0 + c_0101_0^2,
            ...
        
        >>> "PrimaryDecomposition" in p.to_magma()
        True

        Magma file just to compute the Groebner Basis
        >>> "PrimaryDecomposition" in p.to_magma(primary_decomposition = False)
        False
        >>> "GroebnerBasis" in p.to_magma(primary_decomposition = False)
        True
        """
        
        def quote_string(s):

            def quote_line(line):

                def split_chunks(line, chunk_length):
                    for i in range(0, len(line), chunk_length):
                        yield line[i:i+chunk_length]

                return r'\\n'.join(split_chunks(line, chunk_length = 60))
            
            return r'\n'.join([quote_line(line) for line in s.split('\n')])

        magma_format_str_begin = (
            'P<%s> := PolynomialRing(RationalField(), %d);\n'
            'I := ideal<P |\n'
            '%s>;\n'
            '\n'
            '\n'
            'print "==TRIANGULATION" cat "=BEGINS==";\n'
            'print "%s";\n'
            'print "==TRIANGULATION" cat "=ENDS==";\n'
            'print "PY=EVAL=SECTION" cat "=BEGINS=HERE";\n'
            'print "%s";\n'
            'print "PY=EVAL=SECTION=ENDS=HERE";\n'
            '\n'
            'cputime := Cputime();\n'
            '\n')

        magma_format_str_end = ( 
            '\n'
            'print "CPUTIME :", Cputime(cputime);\n')
        
        magma_format_str_primary_decomposition = (
            magma_format_str_begin +
            'print "PRIMARY=DECOMPOSITION" cat "=BEGINS=HERE";\n'
            'P,Q:=PrimaryDecomposition(I);\n'
            'P;\n'
            'print "PRIMARY=DECOMPOSITION=ENDS=HERE";\n' +
            magma_format_str_end)

        magma_format_str_groebner_basis = (
            magma_format_str_begin +            
            'print "GROEBNER=BASIS" cat "=BEGINS=HERE";\n'
            'GroebnerBasis(I);\n'
            'print "GROEBNER=BASIS=ENDS=HERE";\n' +
            magma_format_str_end)

        if primary_decomposition:
            magma_format_str = magma_format_str_primary_decomposition
            eqns = self.equations_with_non_zero_condition
            vars = self.variables_with_non_zero_condition
        else:
            magma_format_str = magma_format_str_groebner_basis
            eqns = self.equations
            vars = self.variables

        ideal_str = ',\n          '.join([str(eqn) for eqn in eqns])

        variables_str = ", ".join(vars)
        num_vars = len(vars)

        triangulation_str = quote_string(self._manifold._to_string())

        py_eval = self.py_eval_section()

        return magma_format_str % (
            variables_str, num_vars, ideal_str, triangulation_str, py_eval)
        
    def filename_base(self):
        """
        Preferred filename base for writing out this Ptolemy variety

        >>> from snappy import *
        >>> M = Manifold('4_1')
        >>> p = M.ptolemy_variety(2, obstruction_class = 1)
        >>> p.filename_base()
        '4_1__sl2_c1'
        """

        base = '%s__sl%d' % (self._manifold.name(), self._N)

        if self._N % 2 == 0 and not self._obstruction_class is None:
            base = base + '_c%d' % self._obstruction_class._index

        return base

    def __repr__(self):
        
        res =  "Ptolemy Variety for %s, N = %d" % (self._manifold.name(), 
                                                   self._N)
        if not self._obstruction_class is None:
            res += ", obstruction_class = %d" % self._obstruction_class._index

        return res

    def compute_solutions(self,
                          engine = None,
                          primary_decomposition = None,
                          memory_limit = 750000000,
                          directory = None,
                          cache_dir = None,
                          verbose = False):

        """
        Starts an engine such as magma to compute the Groebnes basis and
        primary decomposition of the ideal and computes exact solutions.

        If started in sage, uses sage, otherwise, uses magma.

        === Arguments ===

        engine --- engine to use, currently, only support magma and sage
        primary_decomposition --- use primary decomposition, slower but more reliable
        memory_limit --- maximal allowed memory in bytes
        directory --- location for input and output files, temporary directory used if not specified
        verbose --- print extra information
        """

        if engine is None:
            if _within_sage:
                engine = 'sage'
            else:
                engine = 'magma'

        if primary_decomposition is None:
            if engine == 'sage':
                primary_decomposition = False
            else:
                primary_decomposition = True

        if engine == 'magma':
            from . import processMagmaFile
            return processMagmaFile.run_magma(
                self.to_magma(primary_decomposition = primary_decomposition),
                filename_base = self.filename_base(),
                memory_limit = memory_limit,
                directory = directory,
                verbose = verbose)
            
        if engine == 'sage':

            if primary_decomposition:
                sage_prim_decomp = (
                    self.ideal_with_non_zero_condition.primary_decomposition())

                solutions = []
                for component in sage_prim_decomp:
                    sage_gb = component.groebner_basis()
                    gb = [Polynomial.parse_string(str(p)) for p in sage_gb]
                    new_solutions = (
                        solutionsToGroebnerBasis.exact_solutions_with_one(gb))
                    assert len(new_solutions) == 1
                    new_solution = new_solutions[0]
                    assert (
                        (new_solution is None) == 
                        (component.dimension() > 0))
                    solutions.append(new_solution)
                    
            else:
                if len(self.ideal.ring().variable_names()) == 1:
                    # sage doesn't do Groebner basis for an ideal over univariate
                    # polynomial ring
                    # 
                    # this is a principal ideal and we use the one generator
                    assert self.ideal.is_principal()
                    sage_gb = [ self.ideal.gen() ]
                else:
                    sage_gb = self.ideal.groebner_basis()

                gb = [Polynomial.parse_string(str(p)) for p in sage_gb]
                solutions = solutionsToGroebnerBasis.exact_solutions_with_one(gb)

            variable_dict = eval(self.py_eval_variable_dict())

            def process_solution(solution):
                if not solution is None:
                    return coordinates.PtolemyCoordinates(
                        variable_dict(solution), is_numerical = False)
                return None
            
            return [process_solution(solution) for solution in solutions]

        raise "No other engine supported"

def list_all_quadruples_with_fixed_sum(N, skipVerts):

    """
    All quadruples (a,b,c,d) of non-negative integers with a + b + c + d = N
    used to index cross ratios (use N - 2) or Ptolemy coordinates (use
    skipVerts = True).

    >>> list_all_quadruples_with_fixed_sum(2, skipVerts = True)
    [[0, 0, 1, 1], [0, 1, 0, 1], [0, 1, 1, 0], [1, 0, 0, 1], [1, 0, 1, 0], [1, 1, 0, 0]]
    """

    # all quadruples
    all_quads = _enumerate_all_tuples_with_fixed_sum(N, l = 4)
    
    if skipVerts:
        return [quad for quad in all_quads if not N in quad]
    else:
        return [quad for quad in all_quads]

def _fix_decoration(action_by_decoration_change):
        
    action_matrix, ptolemy_coords, decorations_to_be_fixed = (
        action_by_decoration_change)

    fixed_ptolemy_coords = matrix.get_independent_rows(
        action_matrix, ptolemy_coords, len(decorations_to_be_fixed))

    return [(+1, ptolemy_coord, 1) for ptolemy_coord in fixed_ptolemy_coords]

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

    indices = list_all_quadruples_with_fixed_sum(N - 2, skipVerts = False)

    return [generate_ptolemy_relation(tet, index)
            for tet in range(num_tet) for index in indices]

def _non_zero_condition(variables):
    one = Polynomial.constant_polynomial(1)
    
    polynomial = one
    
    for var in variables:
        polynomial = polynomial * Polynomial.from_variable_name(var)
        
    polynomial = polynomial - one

    return polynomial
    
def _enumerate_all_tuples_with_fixed_sum(N, l):
    if l == 1:
        yield [ N ]
    else:
        for i in range(N + 1):
            for j in _enumerate_all_tuples_with_fixed_sum(N-i, l-1):
                yield [i] + j

def _union(lists):
    all = sum(lists, [])
    all = list(set(all))
    all.sort()
    return all

def _identified_variables_canonize(identified_variables):

    def merge_two_dicts(sign, var1, var2, l1, l2):

        total_sign = sign * l1[var1] * l2[var2]

        new_dict = l1
        for key, val in l2.items():
            new_dict[key] = total_sign * val
        
        return new_dict

    all_variables = { }

    for sign, var1, var2 in identified_variables:
        all_variables[var1] = { var1 : + 1 }
        all_variables[var2] = { var2 : + 1 }

    for sign, var1, var2 in identified_variables:
        if not all_variables[var1] is all_variables[var2]:
            new_dict = merge_two_dicts(sign, var1, var2,
                                       all_variables[var1],
                                       all_variables[var2])
            for var in new_dict.keys():
                all_variables[var] = new_dict
                
    result = { }

    for variable, variable_dict in all_variables.items():
        if not result.has_key(variable):
            vars = variable_dict.keys()
            if 1 in vars:
                canonical_rep = 1
            else:
                vars.sort()
                canonical_rep = vars[0]

            for var in vars:
                result[var] = (variable_dict[var] *
                               variable_dict[canonical_rep],
                               canonical_rep)
    
    return result

def _canonical_representative_to_polynomial_substituition(
        canonical_representative):

    result = { }

    for var1, signed_var2 in canonical_representative.items():
        sign, var2 = signed_var2
        if not var1 == var2:

            if var2 == 1:
                result[var1] = (
                    Polynomial.constant_polynomial(sign))
            else:
                result[var1] = (
                    Polynomial.constant_polynomial(sign) *
                    Polynomial.from_variable_name(var2))

    return result
