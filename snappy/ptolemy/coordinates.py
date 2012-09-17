from ptolemyVariety import list_all_quadruples_with_fixed_sum
from solutionsToGroebnerBasis import AlgebraicNumber

try:
    from sage.libs.pari import gen 
    from sage.libs.pari.gen import pari
    from sage.rings.complex_field import ComplexField
    _within_sage = True
except ImportError:
    from cypari import gen
    from cypari.gen import pari
    _within_sage = False

import matrix
import re

class PtolemyCoordinates(dict):
    """
    Represents a solution of a Ptolemy variety as python dictionary.

    === Examples ===

    Construct solution from magma output:

    >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_Magma
    >>> from snappy import Manifold
    >>> solutions = solutions_from_Magma(_magma_output_for_4_1__sl3)
    >>> solution = solutions[2]

    Access a Ptolemy coordinate:

    >>> solution['c_2100_0']
    1

    >>> solution.number_field()
    x^2 - x + 1

    Check that it is really a solution, exactly:
    
    >>> M = Manifold("4_1")
    >>> solution.checkAgainstManifold(M)

    Turn into (Galois conjugate) numerical solutions:

    >>> old_precision = pari.set_real_precision(100) # with high precision
    >>> numerical_solutions = solution.numerical()
    
    Check that it is a solution, numerically:

    >>> numerical_solutions[0].checkAgainstManifold(M, 1e-80)
    >>> pari.set_real_precision(old_precision) # reset pari engine
    100

    Compute cross ratios from the Ptolemy coordinates:

    >>> cross = solution.CrossRatios()
    >>> cross['z_0001_0']
    Mod(-x + 1, x^2 - x + 1)

    Compute volumes:

    >>> volumes = cross.volume()

    Check that volume is 4 times the geometric one:

    >>> volume = volumes[0].abs()
    >>> diff = volume - 4 * M.volume()
    >>> diff < 1e-9
    True

    Compute flattenings:
    
    >>> flattenings = solution.Flattenings()

    Compute complex volumes:

    >>> cvols = [flattening.complexVolume() for flattening in flattenings]
    >>> volume = cvols[0].real().abs()
    >>> chernSimons = cvols[0].imag()
    >>> diff = volume - 4 * M.volume()
    >>> diff < 1e-9
    True

    >>> from snappy import pari
    >>> normalized = chernSimons * 6 / (pari('Pi')**2)

    Check that Chern Simons is zero up to 6 torsion:
    
    >>> normalized - normalized.round() < 1e-9
    True
    """
        
    def __init__(self, d, is_numerical = True):
        super(PtolemyCoordinates, self).__init__(d)
        self._is_numerical = is_numerical
        
    def number_field(self):
        """
        For an exact solution, return the number_field spanned by the
        Ptolemy coordinates. If number_field is Q, return None.
        """
        
        assert not self._is_numerical, "number_field for numerical solution"

        some_value = self.values()[0]
        if some_value.type() == 't_POLMOD':
            return some_value.mod()
        else:
            return None
        

    def numerical(self):
        """
        Turn exact solutions into numerical solutions using pari.

        Take an exact solution:

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_Magma
        >>> solutions = solutions_from_Magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Turn into a numerical solution:
        
        >>> old_precision = pari.set_real_precision(100) # with high precision
        >>> numerical_solutions = solution.numerical()
        >>> pari.set_real_precision(old_precision) # reset pari engine
        100

        Pick one of the Galois conjugates:

        >>> numerical_solution = numerical_solutions[0]
        >>> value = numerical_solution['c_1110_0']
        """
        
        if self._is_numerical:
            return self
        return [PtolemyCoordinates(d, is_numerical = True)
                for d in _to_numerical(self)]

    def CrossRatios(self):
        """
        Compute cross ratios from Ptolemy coordinates. Also see help(ptolemy.CrossRatios).

        Take an exact solution:

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_Magma
        >>> solutions = solutions_from_Magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Turn into cross Ratios:

        >>> crossRatios = solution.CrossRatios()
        
        Get a cross ratio:
        
        >>> crossRatios['zp_0010_0']
        Mod(-x + 1, x^2 - x + 1)

        Get information about what one can do with cross ratios
        """
        return CrossRatios(_ptolemy_to_cross_ratio(self, all_three = True),
                           is_numerical = self._is_numerical)

    def numericalCrossRatios(self):
        """
        Turn exact solutions into numerical and then compute cross ratios.
        See numerical() and CrossRatios().
        """
        
        if self._is_numerical:
            return self.CrossRatios()
        else:
            return [num.CrossRatios() for num in self.numerical()]

    def Flattenings(self):
        """
        Turn into numerical solutions and compute flattenings [z;p,q].
        Also see numerical()

        Get Ptolemy coordinates.

        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_Magma
        >>> solutions = solutions_from_Magma(_magma_output_for_4_1__sl3)
        >>> solution = solutions[2]

        Compute a numerical soluton

        >>> flattenings = solution.Flattenings()

        Get more information with help(flattenings[0])
        """
        if self._is_numerical:
            return Flattenings(
                    _ptolemy_to_cross_ratio(
                        self, with_flattenings = True))
        else:
            return [num.Flattenings() for num in self.numerical()]

    def volume(self, drop_negative_vols = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute volumes.
        If already numerical, only return the one volume.
        See numerical().

        If drop_negative_vols = True is given as optional argument,
        only return non-negative volumes.
        """
        if self._is_numerical:
            return self.CrossRatios().volume()
        else:
            vols = [num.volume() for num in self.numerical()]
            if drop_negative_vols:
                return [vol for vol in vols if vol > -1e-12]
            return vols

    def complexVolume(self, drop_negative_vols = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute complex
        volumes. If already numerical, return the volume.

        Complex volume is defined up to i*pi**2/6.

        See numerical(). If drop_negative_vols = True is given as optional
        argument, only return complex volumes with non-negative real part.
        """
        
        if self._is_numerical:
            return self.Flattenings().complexVolume()
        else:
            cvols = [num.Flattenings().complexVolume()
                     for num in self.numerical()]
            if drop_negative_vols:
                return [cvol for cvol in cvols if cvol.real() > -1e-12]
            return cvols

    def checkAgainstManifold(self, M, epsilon = None):
        """
        Checks that the given solution really is a solution to the Ptolemy
        variety of a manifold. See help(ptolemy.PtolemyCoordinates) for
        more example.

        === Arguments ===

        M --- manifold to check this for
        epsilon --- maximal allowed error when checking the relations, use
        None for exact comparision.
        """

        def get_obstruction_variable(face, tet):
            key = "s_%d_%d" % (face, tet)
            return self[key]

        def check(v, comment):
            if epsilon is None:
                assert v == 0, comment
            else:
                assert v.abs() < epsilon, comment + " error: %s" % v.abs()
        
        N, num_tets, has_obstruction_class = _find_N_tets_obstruction(
            self)

        assert M.num_tetrahedra() == num_tets, "Number tetrahedra not matching"

        if has_obstruction_class:
            # check cocycle condition
            for tet in range(num_tets):
                check( get_obstruction_variable(0, tet) *
                       get_obstruction_variable(1, tet) *
                       get_obstruction_variable(2, tet) *
                       get_obstruction_variable(3, tet) - 1,
                       "Obstruction cocycle condition violated")
            # check identified faces
            for dummy_sign, var1, var2 in (
                    M._ptolemy_equations_identified_face_classes()):
                check ( self[var1] - self[var2],
                        "Identified face classes violated")

        # Check identified Ptolemy coordinates
        for sign, var1, var2 in (
                M._ptolemy_equations_identified_coordinates(N)):
            check (self[var1] - sign * self[var2],
                   "Identified Ptolemy coordinates violated")

        # Check Ptolemy relationship
        indices = list_all_quadruples_with_fixed_sum(N - 2, skipVerts = False)

        for tet in range(num_tets):
            for index in indices:

                def get_Ptolemy_coordinate(addl_index):
                    total_index = matrix.vector_add(index, addl_index)
                    key = "c_%d%d%d%d" % tuple(total_index) + "_%d" % tet
                    return self[key]


                if has_obstruction_class:
                    s0 = get_obstruction_variable(0, tet)
                    s1 = get_obstruction_variable(1, tet)
                    s2 = get_obstruction_variable(2, tet)
                    s3 = get_obstruction_variable(3, tet)
                else:
                    s0 = 1
                    s1 = 1
                    s2 = 1
                    s3 = 1
                    
                rel = (  s0 * s1 * get_Ptolemy_coordinate((1,1,0,0))
                                 * get_Ptolemy_coordinate((0,0,1,1))
                       - s0 * s2 * get_Ptolemy_coordinate((1,0,1,0))
                                 * get_Ptolemy_coordinate((0,1,0,1))
                       + s0 * s3 * get_Ptolemy_coordinate((1,0,0,1))
                                 * get_Ptolemy_coordinate((0,1,1,0)))

                check(rel, "Ptolemy relation violated")

class Flattenings(dict):
    """
    Represents a flattening [z;p,q] assigned to each simplex as dictionary
    """
        
    def __init__(self, d):
        super(Flattenings, self).__init__(d)
        self._is_numerical = True

    def complexVolume(self):
        """
        Compute complex volume.

        Complex volume is defined up to i*pi**2/6.
        """
        p = pari('Pi^2/6')

        cvol = sum([ _L_function(flattening)
                     for flattening in self.values()]) / pari('I')
        vol  = cvol.real()
        cs   = cvol.imag() % p

        if cs > p/2 + pari('1e-12'):
            cs = cs - p

        return vol + cs * pari('I')

class CrossRatios(dict): 
    """
    Represents assigned shape parameters/cross ratios as
    dictionary.
    """
    
    def __init__(self, d, is_numerical = True):
        super(CrossRatios, self).__init__(d)
        self._is_numerical = is_numerical

    def numerical(self):
        """
        Turn exact solutions into numerical solutions using pari. Similar to
        numerical() of PtolemyCoordinates. See help(ptolemy.PtolemyCoordinates)
        for example.
        """        
        if self._is_numerical:
            return self
        return [CrossRatios(d, is_numerical = True) for d in _to_numerical(self)]

    def volume(self, drop_negative_vols = False):
        """
        Turn into (Galois conjugate) numerical solutions and compute volumes.
        If already numerical, only compute the one volume.
        See numerical().

        If drop_negative_vols = True is given as optional argument,
        only return non-negative volumes.
        """
        if self._is_numerical:
            return sum([_volume(z) for key, z in self.items() if 'z_' in key])
        else:
            vols = [num.volume() for num in self.numerical()]
            if drop_negative_vols:
                return [vol for vol in vols if vol > -1e-12]
            return vols

    def checkAgainstManifold(self, M, epsilon = None):
        """
        Checks that the given solution really is a solution to the PGL(N,C) gluing
        equations of a manifold. Usage similar to checkAgainstManifold of
        PtolemyCoordinates. See help(ptolemy.PtolemtyCoordinates) for example.

        === Arguments ===

        M --- manifold to check this for
        epsilon --- maximal allowed error when checking the relations, use
        None for exact comparision.
        """

        def check(v, comment):
            if epsilon is None:
                assert v == 0, comment
            else:
                assert v.abs() < epsilon, comment
        
        some_z = self.keys()[0]
        variable_name, index, tet_index = some_z.split('_')
        assert variable_name in ['z', 'zp', 'zpp']
        assert len(index) == 4
        N = sum([int(x) for x in index]) + 2
        
        matrix_with_explanations = M.gluing_equations_pgl(
            N, equation_type = 'all')

        matrix = matrix_with_explanations.matrix
        rows = matrix_with_explanations.explain_rows
        cols = matrix_with_explanations.explain_columns

        for row in range(len(rows)):
            product = 1
            for col in range(len(cols)):
                cross_ratio_variable = cols[col]
                cross_ratio_value = self[cross_ratio_variable]
                product = product * (cross_ratio_value ** matrix[row,col])
            check(product - 1, "Gluing equation %s violated" % rows[row])
        

def _ptolemy_to_cross_ratio(solution_dict,
                            all_three = False,
                            with_flattenings = False):

    assert not (all_three and with_flattenings)

    N, num_tets, has_obstruction_class = _find_N_tets_obstruction(
        solution_dict)

    def compute_cross_ratios(tet, index):
        def get_Ptolemy_coordinate(addl_index):
            total_index = matrix.vector_add(index, addl_index)
            key = "c_%d%d%d%d" % tuple(total_index) + "_%d" % tet
            return solution_dict[key]

        def get_obstruction_variable(face):
            key = "s_%d_%d" % (face, tet)
            return solution_dict[key]

        strIndicies = '_%d%d%d%d' % tuple(index) + '_%d' % tet
        
        z = ((get_Ptolemy_coordinate((1,0,0,1)) *
               get_Ptolemy_coordinate((0,1,1,0))) /
             (get_Ptolemy_coordinate((1,0,1,0)) *
               get_Ptolemy_coordinate((0,1,0,1))))
        if has_obstruction_class:
            z = z * (get_obstruction_variable(0) *
                     get_obstruction_variable(1))

        if with_flattenings:
            c01 = get_Ptolemy_coordinate((1,1,0,0))
            c02 = get_Ptolemy_coordinate((1,0,1,0))
            c03 = get_Ptolemy_coordinate((1,0,0,1))
            c12 = get_Ptolemy_coordinate((0,1,1,0))
            c13 = get_Ptolemy_coordinate((0,1,0,1))
            c23 = get_Ptolemy_coordinate((0,0,1,1))

            p, q = _compute_flattening(z, c01, c02, c03, c12, c13, c23)

            return [('z' + strIndicies, [z, p, q]),]

        if not all_three:
            return [('z' + strIndicies, z)]

        zp = - ((get_Ptolemy_coordinate((1,1,0,0)) *
                 get_Ptolemy_coordinate((0,0,1,1))) /
                (get_Ptolemy_coordinate((1,0,0,1)) *
                 get_Ptolemy_coordinate((0,1,1,0))))

        if has_obstruction_class:
            zp = zp * (get_obstruction_variable(0) *
                       get_obstruction_variable(2))

        # convention zp and zpp???
        zpp = ((get_Ptolemy_coordinate((0,1,0,1)) *
                get_Ptolemy_coordinate((1,0,1,0))) /
               (get_Ptolemy_coordinate((1,1,0,0)) *
                get_Ptolemy_coordinate((0,0,1,1))))

        if has_obstruction_class:
            zpp = zpp * (get_obstruction_variable(0) *
                         get_obstruction_variable(3))

        return [
            ('z'   + strIndicies, z),
            ('zp'  + strIndicies, zp),
            ('zpp' + strIndicies, zpp)]

    indices = list_all_quadruples_with_fixed_sum(N - 2, skipVerts = False)

    return dict(
        sum([compute_cross_ratios(tet, index)
          for tet in range(num_tets) for index in indices],[]))

    # looks in solution_dict keys to look for obstruction class and
    # ptolemy's
    # apply formula for cross_ratio

    # if with_flattenings, compute w0, w1, w2
    
    pass

def _find_N_tets_obstruction(solution_dict):
    N = None
    num_tets = 0
    has_obstruction_class = False

    for k in solution_dict.keys():
        variable_name, index, tet_index = k.split('_')
        assert variable_name in ['c', 's']
        num_tets = max(num_tets, int(tet_index)+1)
        if variable_name == 'c': # We are in the Ptolemy case
            assert len(index) == 4
            new_N = sum([int(x) for x in index])
            if N is None:
                N = new_N
            else:
                assert N == new_N
        else:
            has_obstruction_class = True
            
    return N, num_tets, has_obstruction_class

def _has_no_number_field(d):
    for key, value in d.items():
        if re.match('Mod\(.*,.*\)', str(value)):
            return False
    return True
            
def _to_numerical(d):
    if _has_no_number_field(d):
        return [d]
    else:
        return _to_numerical_iter(d)

def _to_numerical_iter(d):

    number_field = None
    new_dict = { }

    for key, value in d.items():
        if re.match('Mod\(.*,.*\)', str(value)):
            new_dict[key] = AlgebraicNumber.from_pari(value)
            number_field = new_dict[key].number_field
        else:
            new_dict[key] = value

    if number_field is None:
        roots = [ pari(0) ]
    else:
        # Bug in cypari: pari(str(number_field)).polroots()
        # gives less precision
        roots = pari('polroots(%s)' % number_field)

    for root in roots:

        def to_numerical(value):
            if re.match('Mod\(.*,.*\)', str(value)):
                return value.to_numerical(root)
            else:
                return value

        yield dict([ (key,to_numerical(value))
                     for key, value in new_dict.items()])

def _convert_to_pari_float(z):

    if z.type() in ['t_INT', 't_FRAC']:
        return z * pari('1.0')
    
    return z
 
def _compute_flattening(z, c01, c02, c03, c12, c13, c23):

    z = _convert_to_pari_float(z)
    c01 = _convert_to_pari_float(c01)
    c02 = _convert_to_pari_float(c02)
    c03 = _convert_to_pari_float(c03)
    c12 = _convert_to_pari_float(c12)
    c13 = _convert_to_pari_float(c13)
    c23 = _convert_to_pari_float(c23)

    log_c01 = (c01**2).log()/2
    log_c02 = (c02**2).log()/2
    log_c03 = (c03**2).log()/2
    log_c12 = (c12**2).log()/2
    log_c13 = (c13**2).log()/2
    log_c23 = (c23**2).log()/2

    w0 = log_c03 + log_c12 - log_c02 - log_c13
    w1 = log_c02 + log_c13 - log_c01 - log_c23

    PiI = pari('Pi * I')

    p = ((w0 -    z .log()) / PiI).round()
    q = ((w1 + (1-z).log()) / PiI).round()

    return p, q

# bug in pari

def _dilog(z):
    return pari("dilog(%s)" % z)

def _L_function(flattening):

    z, p, q = flattening

    z = _convert_to_pari_float(z)
    p = _convert_to_pari_float(p)
    q = _convert_to_pari_float(q)

    PiI = pari('Pi * I')
    Pi2 = pari('Pi * Pi')

    return (  _dilog(z)
            + (z.log() + p * PiI) * ((1-z).log() + q * PiI) / 2
            - Pi2 / 6)

def _volume(z):
    
    z = _convert_to_pari_float(z)
    
    return (1-z).arg() * z.abs().log() + _dilog(z).imag()
