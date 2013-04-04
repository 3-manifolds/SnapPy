import os, sys, operator, types, re, gzip, struct, tempfile
import tarfile, atexit, math, string
from signal import signal, SIGINT, default_int_handler
from manifolds import __path__ as manifold_paths
import database

import twister

include "SnapPy.pxi"

# This is part of the UCS2 hack.
cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    return string

# A stream for asynchronous messages
class MsgIO:
    def __init__(self):
        self.write = sys.stdout.write
        self.flush = sys.stdout.flush()

msg_stream = MsgIO()

# We need a matrix class
class SimpleMatrix:
    """
    A very simple matrix class that wraps a list of lists.  It has
    two indices and can print itself.  Nothing more.
    """
    def __init__(self, list_of_lists):
        self.data = list_of_lists
        self.type = type(self.data[0][0])
        self.shape = (len(list_of_lists), len(list_of_lists[0]))

    def __repr__(self):
        str_matrix = [[str(x) for x in row] for row in self.data]
        size = max([max([len(x) for x in row]) for row in str_matrix])
        str_rows = []
        for row in str_matrix:
            str_row = ['% *s'%(size, x) for x in row]
            str_rows.append('[' + ', '.join(str_row) + ']')
        result = 'matrix([' + ',\n        '.join(str_rows) + '])'
        return result

    def __str__(self):
        str_matrix = [[str(x) for x in row] for row in self.data]
        size = max([max([len(x) for x in row]) for row in str_matrix])
        str_rows = []
        for row in str_matrix:
            str_row = ['% *s'%(size, x) for x in row]
            str_rows.append(' [' + ' '.join(str_row) + ']')
        result = '\n'.join(str_rows)
        result = '[' + ('\n'.join(str_rows))[1:] + ']'
        return result

    def _check_indices(self, key):
        if type(key) == slice:
            raise TypeError("Simple matrices don't slice.")
        if type(key) == int:
            raise TypeError('Simple matrices need 2 indices.')
        i, j = key
        if i < 0 or j < 0:
            raise TypeError("Simple matrices don't have negative indices.") 
        return key

    def __getitem__(self, key):
        i, j = self._check_indices(key)
        return self.data[i][j]

    def __setitem__(self, key, value):
        i, j = self._check_indices(key)
        self.data[i][j] = value

    def _noalgebra(self, other):
        raise TypeError('To do matrix algebra, please install numpy '
                        'or run SnapPy in Sage.')
    def entries(self):
        return [x for row in self.data for x in row]

    __add__ = __sub__ = __mul__ = __div__ = __inv = _noalgebra

# Sage interaction
try:
    import sage.all
    import sage.structure.sage_object
    from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.interfaces.gap import gap
    from sage.interfaces.gap import is_GapElement
    from sage.interfaces.magma import magma
    from sage.interfaces.magma import is_MagmaElement
    from sage.matrix.constructor import matrix
    from sage.libs.pari.gen import pari as pari
    _within_sage = True
except ImportError:
    from cypari.gen import pari as pari
    matrix = SimpleMatrix
    _within_sage = False

# High precision stuff.

from . import snap

# Ptolemy variety stuff

from .ptolemy import manifoldMethods as ptolemyManifoldMethods
 
# The next two functions provide replacements for code in
# the SnapPea kernel module Dirichlet_precision.c which
# attempts to deal with round-off error when multiplying O31 matrices.

cdef public void precise_o31_product( O31Matrix a, O31Matrix b, O31Matrix product):
    # For now, just use SnapPea's built-in double-precision product
    o31_product(a, b, product);

cdef public void precise_generators( MatrixPairList* gen_list):
    # We don't need this at the moment.
    return

# Enable graphical link input
try:
    from plink import LinkEditor
except:
    LinkEditor = None

# Enable OpenGL display of DirichletDomains
try:
    from snappy.polyviewer import PolyhedronViewer
except ImportError:
    PolyhedronViewer = None

# Enable OpenGL display of CuspNeighborhoods
try:
    from horoviewer import HoroballViewer
except ImportError:
    HoroballViewer = None

# Enable Tk based save dialog

try:
    from snappy.filedialog import asksaveasfile
except ImportError:
    asksaveasfile = None

# String testing for Python 2 and 3.  
try:
    unicode
    def to_str(s):
        return s
    def bytearray_to_bytes(x):
        return str(x)
    
except NameError: # Python 3
    basestring = unicode = str
    def to_str(s):
        return s.decode()
    def bytearray_to_bytes(x):
        return bytes(x)

def to_byte_str(s):
    return s.encode('utf-8') if type(s) != bytes else s

# Paths
manifold_path = manifold_paths[0] + os.sep
# Obsolete
closed_census_directory = os.path.join(manifold_path, 'ClosedCensusData')
link_directory = os.path.join(manifold_path, 'ChristyLinks')
link_archive = os.path.join(manifold_path, 'ChristyLinks.tgz')
census_knot_archive = os.path.join(manifold_path, 'CensusKnots.tgz')
# Still in use
table_directory = os.path.join(manifold_path, 'HTWKnots')
morwen_link_directory = os.path.join(manifold_path, 'MTLinks')

# These are the gzipped files holding the knot tables.
Alternating_table = gzip.open(os.path.join(table_directory, 'alternating.gz') )
Nonalternating_table = gzip.open(os.path.join(table_directory, 'nonalternating.gz') )

# This flag can be used when testing algorithms in extension
# modules.  In the C code, declare:
# extern int SnapPy_test_flag

cdef public int SnapPy_test_flag = 0

def set_test_flag(int value):
    cdef int old
    global SnapPy_test_flag
    old = SnapPy_test_flag
    SnapPy_test_flag = value
    return old 

# Implementation of the SnapPea UI functions and their global variables.
cdef extern from *:
    ctypedef char* const_char_ptr "const char*"
    ctypedef int const_int "const int"

class SnapPeaFatalError(Exception):
    """
    This exception is raised by SnapPy when the SnapPea kernel
    encounters a fatal error.
    """

cdef public void uFatalError(char *function, char *file) except *:
    raise SnapPeaFatalError('SnapPea crashed in function %s(), '
                            'defined in %s.c.'%(function, file))

# global variables for interrupt processing 
cdef public Boolean gLongComputationInProgress
cdef public Boolean gLongComputationCancelled

# If not None, this will be called in gLongComputationContinues.
# This enables a GUI to do updates during long computations.
UI_callback = None

#cdef void SnapPea_sigint_handler(int signum, frame=None):
def SnapPea_sigint_handler(int signum, frame=None):
    """
    A SIGINT handler which cancels the SnapPea computation.
    """
    global gLongComputationCancelled
    gLongComputationCancelled = True
    sys.stderr.write('\nSnapPea computation aborted!\n')

def SnapPea_interrupt():
    """
    The UI can call this to stop SnapPea.  Returns True if SnapPea was busy.
    """
    global gLongComputationCancelled
    global gLongComputationInProgress
    gLongComputationCancelled = True
    return gLongComputationInProgress

cdef public void uLongComputationBegins(char *message, Boolean is_abortable):
    global gLongComputationCancelled
    global gLongComputationInProgress
    global old_sigint_handler, old_sigalrm_handler
    # Set SnapPea's flags
    gLongComputationCancelled = False
    gLongComputationInProgress = True
    # Install our sigint handler
    signal(SIGINT, SnapPea_sigint_handler)

cdef public c_FuncResult uLongComputationContinues() except *:
    global gLongComputationCancelled
    global gLongComputationInProgress
    if gLongComputationCancelled:
        gLongComputationCancelled = False
        gLongComputationInProgress = False
        return func_cancelled
    else:
        if UI_callback is not None:
            # Let the GUI update itself
            UI_callback()
        return func_OK

cdef public void uLongComputationEnds():
    global gLongComputationCancelled
    global gLongComputationInProgress
    # Restore the Python sigint handler
    signal(SIGINT, default_int_handler)
    # Reset SnapPea's flags
    gLongComputationCancelled = False
    gLongComputationInProgress = False

cdef public void uAcknowledge(const_char_ptr message):
    sys.stderr.write(<char *> message)
    sys.stderr.write('\n')

cdef public void uAbortMemoryFull():
    sys.stderr.write('Out of memory.\n')
    sys.exit(2)

cdef public int uQuery(const_char_ptr  message, 
                       const_int       num_responses,
                       const_char_ptr  responses[],
                       const_int       default_response):
    #  If desired you could write this function to obtain a response
    #  from the user, but for now it is set up to return the default
    #  response, to facilitate batch computations.
    cdef char *default = <char *> responses[<int> default_response]
    sys.stderr.write('Q: %s\nA:  %s\n'%(<char *> message, default))
    return <int> default_response

def cy_eval(s):
    return eval(s)

def smith_form(M):
    if _within_sage:
        if not hasattr(M, 'elementary_divisors'):
            M = matrix(M) 
        m, n = M.nrows(), M.ncols()
        result = M.elementary_divisors(algorithm='pari')
    else:
        if not isinstance(M, matrix):
            M = matrix(M)
        m, n = M.shape
        result = [int(x) for x in pari.matrix(m, n, M.entries()).matsnf()]

    # PARI views the input to matsnf0 as square.
    if m < n:
        result = result + [0]*(n-m)
    if m > n:
        for i in range(m - n):
            result.remove(0)

    # For consistency with SnapPea, need to switch the order of the factors.
    zeros = [x for x in result if x == 0]
    nonzeros = [x for x in result if x != 0]
    nonzeros.sort()
    return nonzeros + zeros

# end smith form code


# Enum conversions
CuspTopology = ['torus cusp', 'Klein bottle cusp', 'unknown']
MatrixParity = ['orientation-reversing', 'orientation-preserving']
Orientability = ['orientable', 'nonorientable', 'unknown']
Orbifold1 = ['unknown', 'circle', 'mirrored arc']
FuncResult = ['func_OK', 'func_cancelled', 'func_failed', 'func_bad_input']
SolutionType = ['not attempted', 'all tetrahedra positively oriented',
                'contains negatively oriented tetrahedra', 'contains flat tetrahedra',
                'contains degenerate tetrahedra', 'unrecognized solution type',
                'no solution found']

# global functions
def check_SnapPea_memory():
    verify_my_malloc_usage()

# convert and free an identification of variables structure

cdef convert_and_free_identification_of_variables(
                                        Identification_of_variables c_vars):
    var_list = []

    if c_vars.variables:
        for i in range(c_vars.num_identifications):
            var_list.append(
                  (c_vars.signs[i],
                   str(c_vars.variables[i][0]),
                   str(c_vars.variables[i][1])))
    return var_list

# convert and free an integer matrix from C

cdef convert_and_free_integer_matrix(Integer_matrix_with_explanations c_matrix):
    if not c_matrix.entries:
        return []

    python_matrix = [
        [ c_matrix.entries[i][j] for j in range(c_matrix.num_cols)]
                                 for i in range(c_matrix.num_rows)]

    explain_row = []

    for i in range(c_matrix.num_rows):
        if c_matrix.explain_row[i]:
            explain_row.append(str(c_matrix.explain_row[i]))
        else:
            explain_row.append(None)

    explain_column = []

    for i in range(c_matrix.num_cols):
        if c_matrix.explain_column[i]:
            explain_column.append(str(c_matrix.explain_column[i]))
        else:
            explain_column.append(None)	   

    free_integer_matrix_with_explanations(c_matrix)

    return python_matrix, explain_row, explain_column

class MatrixWithExplanations(object):

    def __init__(self, mat, explain_rows, explain_columns):

        self.matrix = mat
        self.explain_rows = explain_rows
        self.explain_columns = explain_columns

    def __add__(self, other):

        assert self.explain_columns == other.explain_columns, (
	    "matrices with different columns")

        if isinstance(self.matrix, SimpleMatrix):
            newMatrix = SimpleMatrix(self.matrix.data + other.matrix.data)
        elif _within_sage:
            newMatrix = self.matrix.stack(other.matrix)
        else:
            raise ValueError('Matrix type in MatrixWithExplanations '
	                     'not supported')

        return MatrixWithExplanations(
	    newMatrix,
	    self.explain_rows+other.explain_rows,
            self.explain_columns)

    def __repr__(self, _type_str = "MatrixWithExplanations"):

        def format_explain_list(l):
            if len(l) <= 6:
                return repr(l)
            
            return "[ %s, ..., %s]" % (
                ', '.join(map(repr,l[:3])),
                ', '.join(map(repr,l[-3:])))
                

        return (
            "%s(\n"
            "  %s,\n"
            "  explain_columns = %s,\n"
            "  explain_rows = %s)") % (
	                      _type_str,
                              '\n  '.join(repr(self.matrix).split('\n')),
                              format_explain_list(self.explain_columns),
                              format_explain_list(self.explain_rows))

class NeumannZagierTypeEquations(MatrixWithExplanations):

    def __init__(self, mat, explain_rows, explain_columns):
        MatrixWithExplanations.__init__(self, 
	                                mat, explain_rows, explain_columns)

    def __repr__(self):
        return MatrixWithExplanations.__repr__(self,
	                                       "NeumannZagierTypeEquations")
					 
    def __add__(self, other):
        mat = MatrixWithExplanations.__add__(self, other)

        return NeumannZagierTypeEquations(
            mat.matrix, mat.explain_rows, mat.explain_columns)

# Derivatives of basic classes; when these are called as a function
# they return self.  This means that they can be accessed either as
# attributes or as (fake) getter methods.
# These classes use __slots__ in order to save space by not creating
# a __dict__.  The only slot defined at the moment is accuracy.

class SnapPyStr(str):
    __slots__ = []
    def __call__(self):
        return self

class SnapPyBoolean(int):
    __slots__ = []
    def __repr__(self):
        return 'True' if self else 'False'
    def __call__(self):
        return self

class SnapPyInt(int):
    __slots__ = []
    def __call__(self):
        return self

class SnapPyFloat(float):
    __slots__ = ['accuracy']
    def __call__(self):
        return self
    def __repr__(self):
        try:
            return ('{0:.%sf}'%(self.accuracy+1)).format(self)
        except AttributeError:
            return float.__repr__(self)
    def __str__(self):
        return self.__repr__()
        
class SnapPyComplex(complex):
    __slots__ = ['accuracy']
    def __call__(self):
        return self
    def __repr__(self):
        try:
            D = self.accuracy+1
            return ('({0.real:.%sf}{0.imag:+.%sf}j)'%(D,D)).format(self)
        except AttributeError:
            return complex.__repr__(self)
    def __str__(self):
        return self.__repr__()


class SnapPyList(list):
    __slots__ = []
    def __call__(self):
        return self

# Immutable containers which hold information about SnapPy objects.
# The base class for these is Info. Subclasses should override
# __repr__ to appropriately display the information they contain.

class Info(dict):
    """
    Immutable dictionary whose data can be accessed either as
    attributes or by mapping.
    Intialize with keyword arguments, or **dict.
    """
    def __init__(self, **kwargs):
        # Hack to support obsolete keys
        content = dict(kwargs)
        for old, new in self._obsolete.items():
            try:
                content[old] = content[new]
            except KeyError:
                pass
        super(Info, self).__init__(content)
        self.__dict__.update(kwargs)
    def _immutable(self, *args):
        raise AttributeError('Info objects are immutable.')
    def keys(self):
        return self.__dict__.keys()
    __setattr__ = __delattr__ = __setitem__ = __delitem__ = _immutable
    pop = popitem = clear = update = _immutable
    _obsolete = {}

class CuspInfo(Info):
    def __repr__(self):
        if self.is_complete:
            if 'shape' in self:
                return ('Cusp %-2d: complete %s of shape %s' %
                        (self.index, self.topology, self.shape) )
            else:
                return ('Cusp %-2d: %s, not filled'%
                        (self.index, self.topology) )
        else:
            return ('Cusp %-2d: %s with Dehn filling coeffients (M, L) = %s'%
                    (self.index, self.topology, self.filling) )
    _obsolete = {'complete?'          : 'is_complete',
                 'holonomy precision' : 'holonomy_accuracy', 
                 'shape precision'    : 'shape_accuracy'}
    
class DualCurveInfo(Info):
    def __repr__(self):
        return ('%3d: %s curve of length %s'%
                (self.index, MatrixParity[self.parity], self.filled_length))
    _obsolete = {'complete length' : 'complete_length',
                 'filled length' : 'filled_length'}
    
class LengthSpectrumInfo(Info):
    def __repr__(self):
        return '%-4d %-32s %-14s%s'%(
            self.multiplicity, self.length, self.topology, self.parity )

class ShapeInfo(Info):
    _obsolete = {'precision' : 'accuracies'}
    def __repr__(self):
        return repr(self.__dict__)
    
class LengthSpectrum(list):
    def __repr__(self):
        base = ['%-4s %-32s %-12s  %s'%
                ('mult', 'length', 'topology', 'parity')]
        return '\n'.join(base + [repr(s) for s in self])

class ListOnePerLine(list):
    def __repr__(self):
        return '[' + ',\n '.join([repr(s) for s in self]) + ']'

# Abelian Groups

cdef class AbelianGroup:
    """
    An AbelianGroup object represents a finitely generated abelian group,
    usually the first homology group of a snappy Manifold.

    Instantiate as AbelianGroup(P) where P is a presentation matrix
    given as a list of lists of integers.  Alternatively, use
    AbelianGroup(elementary_divisors=[n_1, n_2, ... ]) where the n_i
    are the elementary divisors of the group.

    >>> AbelianGroup([[1,3,2],[2,0,6]])
    Z/2 + Z
    >>> A = AbelianGroup(elementary_divisors=[5,15,0,0])
    >>> A
    Z/5 + Z/15 + Z + Z
    >>> A[1]
    15
    >>> A.betti_number()
    2
    >>> A.order()
    'infinite'
    >>> len(A)
    4
    """
    cdef divisors
    # Backwards compatibility hack, part 1.
    cdef public coefficients
    
    def __init__(self, presentation=None, elementary_divisors=[]):
        if presentation is not None:
            self.divisors = smith_form(presentation)
        else:
            try:
                self.divisors = list(elementary_divisors)
            except:
                raise ValueError('Elementary divisors must be given '
                                 'as a sequence.')
        int_types = [int, long]
        if _within_sage:
            int_types += [sage.rings.integer.Integer]
        for c in self.divisors:
            assert type(c) in int_types and c >= 0,\
                   'Elementary divisors must be non-negative integers.\n'
        for i in range(len(elementary_divisors) - 1):
            n,m = elementary_divisors[i:i+2]
            assert (n == m == 0) or (m % n == 0),\
                   'The elementary divisors must form a divisibility chain\n'

        # So that the group is determined entirely by self.divisors
        # we don't allow '1' as a divisor.  
        self.divisors = [n for n in self.divisors if n != 1]
        # Backwards compatibility hack, part 2.
        self.coefficients = self.divisors

    def __repr__(self):
        if len(self.divisors) == 0:
            return '0'
        factors = ( ['Z/%d'%n for n in self.divisors if n != 0] +
                    ['Z' for n in self.divisors if n == 0] )
        return ' + '.join(factors)

    def __len__(self):
        return len(self.divisors)
    
    def __getitem__(self, i):
        return self.divisors[i]

    def __richcmp__(AbelianGroup self, AbelianGroup other, op):
        if op == 0:
            return self.divisors < other.elementary_divisors()
        elif op == 2:
            return self.divisors == other.elementary_divisors()
        elif op == 4:
            return self.divisors > other.elementary_divisors()
        elif op == 1:
            return self.divisors <= other.elementary_divisors()
        elif op == 3:
            return self.divisors != other.elementary_divisors()
        elif op == 5:
            return self.divisors >= other.elementary_divisors()
        else:
            return NotImplemented
        
    def __call__(self):
        return self

    def elementary_divisors(self):
        """
        The elementary_divisors of this finitely generated abelian group.
        """
        return SnapPyList(self.divisors)

    def rank(self):
        """
        The rank of the group.
        """
        return SnapPyInt(len(self.divisors))

    def betti_number(self):
        """
        The rank of the maximal free abelian subgroup.
        """
        return SnapPyInt(len([n for n in self.divisors if n == 0]))

    def order(self):
        """
        The order of the group.  Returns the string 'infinite' if the
        group is infinite.        
        """
        det = 1
        for c in self.divisors:
            det = det * c
        return SnapPyStr('infinite') if det == 0 else SnapPyInt(det)

# Isometry

def format_two_by_two(mat):
    a,b,c,d = ['%d' % x for x in [mat[0,0], mat[0,1], mat[1,0], mat[1,1]]]
    w0 = max(len(a), len(c))
    w1 = max(len(b), len(d))
    return ('[' + a.rjust(w0) + ' ' + b.rjust(w1) + ']',
            '[' + c.rjust(w0) + ' ' + d.rjust(w1) + ']')
    
class Isometry():
    """
    Represents an isometry from one manifold to another.
    """
    def __init__(self, cusp_images, cusp_maps, extends_to_link):
        self._cusp_images, self._cusp_maps = cusp_images, cusp_maps
        self._extends_to_link = extends_to_link

    def cusp_images(self):
        return self._cusp_images

    def cusp_maps(self):
        return self._cusp_maps

    def extends_to_link(self):
        return self._extends_to_link

    def num_cusps(self):
        return len(self.cusp_images())

    def __repr__(self):
        images = self.cusp_images()
        maps = self.cusp_maps()
        line0, line1, line2 = [], [], []
        for i in range(self.num_cusps()):
            l0 = '%d -> %d' % (i, images[i])
            l1, l2 = format_two_by_two(maps[i])
            L = max(len(l0), len(l1), len(l2))
            line0.append( l0.ljust(L) )
            line1.append( l1.ljust(L) )
            line2.append( l2.ljust(L) )
        line3 = 'Extends to link' if self.extends_to_link() else 'Does not extend to link'
        return '\n'.join(['  '.join(line0), '  '.join(line1), '  '.join(line2), line3])

cdef IsometryListToIsometries(IsometryList *isometries):
    cdef int n, c, i, j, c_cusp_image
    cdef MatrixInt22  c_cusp_map
    cdef Boolean extends
    n = isometry_list_size(isometries)
    c = isometry_list_num_cusps(isometries)

    ans = []
    for i in range(n):
        cusp_images = []
        cusp_maps = []
        for j in range(c):
            isometry_list_cusp_action(isometries, i, j, &c_cusp_image, c_cusp_map)
            cusp_images.append(c_cusp_image)
            cusp_maps.append(matrix(
                [ [c_cusp_map[0][0], c_cusp_map[0][1]],
                  [c_cusp_map[1][0], c_cusp_map[1][1]]] ))

        ans.append(Isometry(cusp_images, cusp_maps, B2B(isometry_extends_to_link(isometries, i))))

    return ans

# SnapPea Classes

# Triangulations

cdef class Triangulation(object):
    """
    A Triangulation object represents a compact 3-manifold with torus
    boundary components, given as an ideal triangulation of the
    manifold's interior.  A Dehn-filling can be specified for each
    boundary component, allowing the description of closed 3-manifolds
    and some orbifolds.  For non-orientable 3-manifolds, the boundary
    components can also be Klein bottles. Two Triangulations are equal
    ('==') if they represent combinatorially isomorphic
    triangulations.  A Triangulation does *not* have any geometric
    structure, and usually one works with the subclass Manifold which
    adds this.  Here's a quick example:

    >>> M = Triangulation('9_42')
    >>> M.num_tetrahedra()
    5
    >>> M.is_orientable()
    True

    A Triangulation can be specified in a number of ways, e.g.

    - Triangulation('9_42') : The complement of the knot 9_42 in S^3.
    - Triangulation('m125(1,2)(4,5)') : The SnapPea census manifold m125
      where the first cusp has Dehn filling (1,2) and the second cusp has
      filling (4,5).
    - Triangulation() : Opens a link editor window where can you
      specify a link complement.
    
    In general, the specification can be from among the below, with
    information on Dehn fillings added.

    - SnapPea cusped census manifolds: e.g. 'm123', 's123', 'v123'.

    - Link complements:
        + Rolfsen's table: e.g. '4_1', '04_1', '5^2_6', '6_4^7', 'L20935', 'l104001'.
        + Knots and links up to 14 crossings from tabulations by Hoste
          and Thistlethwaite: e.g. 'K12a456' or 'L13n579'.
        + Hoste-Thistlethwaite Knotscape table:  e.g. '11a17' or '12n345'
        + Dowker-Thistlethwaite code: e.g. 'DT[6,8,2,4]', 'DT[dadbcda]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a_0, B_1])' or 'Splitting(S_{1,0}, [b_1, A_0], [a_0,B_1])'

      See the help for the 'twister.twister' function for more.  

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """
    cdef c_Triangulation* c_triangulation
    cdef _DTcode
    cdef readonly _cache
    cdef readonly LE

    def __cinit__(self, spec=None, DTcode=None ):
        cdef c_Triangulation *c_triangulation = NULL
        self._DTcode = DTcode
        # Answers to potentially hard computations are cached
        self._cache = {}
        self.LE = None
        if spec is not None and spec != 'empty':
            if not isinstance(spec, basestring):
                raise TypeError(triangulation_help%
                                self.__class__.__name__)
            c_triangulation = get_triangulation(spec)
            if c_triangulation == NULL:
                raise RuntimeError, 'An empty triangulation was generated.'
        if spec is None:
            # Try to determine the name of the variable associated
            # to the manifold:
            if LinkEditor:
                try:
                    IP = cy_eval('get_ipython()')
                    fallback = 'Out[%d]'%IP.execution_count
                    cmd = IP._last_input_line
                    m = re.match('\s*([a-zA-Z_0-9]+)\s*=\s*Manifold\(\)', cmd)
                    link_title = m.group(1) if m else fallback
                except NameError:
                    link_title = None
                LE = LinkEditor(no_arcs=True,
                                callback=self._plink_callback,
                                cb_menu='Send to SnapPy')
                if link_title:
                    print('Starting the link editor.\n'
                          'Select PLink->Send to SnapPy to load the '
                          'link complement as the variable %s' % link_title)

                    LE.window.title('PLink Editor - %s' % link_title)

                else:
                    print('Starting the link editor.\n'\
                          'Select PLink->Send to SnapPy to load the link complement.')
                self.LE = LE

            else:
                raise RuntimeError, 'PLink was not imported.'

        if c_triangulation != NULL:    
            self.set_c_triangulation(c_triangulation)
            remove_hyperbolic_structures(c_triangulation)

    def clear_cache(self, key=None):
        if not key: 
            self._cache.clear()
        else:
            self._cache.pop(key)

    def _plink_callback(self):
        cdef c_Triangulation* c_triangulation = NULL
        if self.LE is not None:
            klp = self.LE.SnapPea_KLPProjection() 
            if klp is not None:
                c_triangulation = get_triangulation_from_PythonKLP(klp)
                self.set_c_triangulation(c_triangulation)
                self._cache = {}
                msg_stream.write('\nNew triangulation received from PLink!\n')
                return
        else:
            raise RuntimeError, 'Communication with PLink failed.' 

    def plink(self):
        """
        Brings the corresponding link editor window to the front, if
        there is one.
        """
        if self.LE is not None:
            self.LE.reopen()
        else:
            raise ValueError, 'This manifold does not have a PLink window.'

    cdef set_c_triangulation(self, c_Triangulation* c_triangulation):
        self.c_triangulation = c_triangulation

    def num_cusps(self, cusp_type='all'):
        """
        Return the total number of cusps.  By giving the optional argument
        'orientable' or 'nonorientable' it will only count cusps of that type.

        >>> M = Triangulation('m125')
        >>> M.num_cusps()
        2
        """
        if cusp_type == 'all':
            return get_num_cusps(self.c_triangulation)
        elif cusp_type == 'orientable':
            return get_num_or_cusps(self.c_triangulation)
        elif cusp_type == 'nonorientable':
            return get_num_nonor_cusps(self.c_triangulation)
        else:
            raise ValueError("Acceptable cusp types are "
                             "['all', 'orientable', 'nonorientable'].")

    def orientation_cover(self):
        """
        For a non-orientable manifold, returns the 2-fold cover which
        is orientable.

        >>> X = Triangulation('x123')
        >>> Y = X.orientation_cover()
        >>> (X.is_orientable(), Y.is_orientable())
        (False, True)
        >>> Y
        x123~(0,0)(0,0)
        """ 
        
        if self.is_orientable():
            raise ValueError, 'The Triangulation is already orientable.'

        cdef c_Triangulation* cover_c_triangulation = NULL
        cdef Triangulation new_tri
        
        cover_c_triangulation = double_cover(self.c_triangulation)
        new_tri = Triangulation('empty')
        new_tri.set_c_triangulation(cover_c_triangulation)
        new_tri.set_name(self.name() + '~')
        return new_tri

    def is_orientable(self):
        """
        Return whether the underlying 3-manifold is orientable.

        >>> M = Triangulation('x124')
        >>> M.is_orientable()
        False
        """
        orientability = Orientability[get_orientability(self.c_triangulation)]
        if orientability == 'orientable': return SnapPyBoolean(True)
        elif orientability == 'nonorientable': return SnapPyBoolean(False)
        else: return None

    def copy(self):
        """
        Returns a copy of the triangulation.

        >>> M = Triangulation('m125')
        >>> N = M.copy()
        """
        cdef c_Triangulation* copy_c_triangulation = NULL
        cdef Triangulation new_tri
        
        copy_triangulation(self.c_triangulation, &copy_c_triangulation)
        new_tri = Triangulation('empty')
        new_tri.set_c_triangulation(copy_c_triangulation)
        return new_tri

    def randomize(self):
        """
        Perform random Pachner moves on the underlying triangulation.

        >>> M = Triangulation('braid[1,2,-3,-3,1,2]')
        >>> M.randomize()
        """
        if self.c_triangulation is NULL: return
        randomize_triangulation(self.c_triangulation)
        self._cache = {}

    def simplify(self):
        """
        Try to simplify the triangulation by doing Pachner moves.

        >>> M = Triangulation('12n123')
        >>> M.simplify()
        """
        if self.c_triangulation is NULL: return
        basic_simplification(self.c_triangulation)
        self._cache = {}

    def _two_to_three(self, tet_num, face_index):
        cdef c_FuncResult result
        cdef c_Tetrahedron* tet

        n = range(self.num_tetrahedra())[tet_num]
        tet = self.c_triangulation.tet_list_begin.next
        for i in range(n):
            tet = tet.next
        result = two_to_three(tet, face_index, &self.c_triangulation.num_tetrahedra)
        return result
        
    def with_hyperbolic_structure(self):
        """
        Add a (possibly degenerate) hyperbolic structure, turning the
        Triangulation into a Manifold.

        >>> M = Triangulation('m004')
        >>> N = M.with_hyperbolic_structure()
        >>> N.volume()
        2.02988321282
        """
        return Manifold_from_Triangulation(self)

    def _empty_save(self):
        """
        Called by save when no file name is specified, so that in
        theory it can be overloaded by the UI.
        """
        if asksaveasfile:
            savefile = asksaveasfile(
                mode='w', title='Save Triangulation', defaultextension='.tri',
                filetypes = [
                ('Triangulation and text files', '*.tri *.txt', 'TEXT'),
                ('All text files', '', 'TEXT'),
                ('All files', '')])
            if savefile:
                filename = savefile.name
                savefile.close()
                self.save(filename)
                return 

        raise ValueError, 'Please specify a file name.'
        
    def save(self, file_name=None):
        """
        Save the triangulation as a SnapPea triangulation file.

        >>> M = Triangulation('m004')
        >>> M.save('fig-eight.tri')
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if file_name == None:
            self._empty_save()
        else:
            b_file_name = to_byte_str(file_name)
            write_triangulation(self.c_triangulation, b_file_name)

    def _to_string(self):
        """
        Return a string containing the contents of a SnapPea triangulation
        file.
	>>> M = Manifold('7_4')
	>>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N == M
        True
        """
        cdef char *c_string
        cdef result
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        else:
            try:
                c_string = string_triangulation(self.c_triangulation)
                result = c_string
            finally:
                free(c_string)
            return to_str(result)

    def _from_string(self, string):
        """
        Fill an empty Triangulation from a string generated by _to_string.
	>>> M = Manifold('7_4')
	>>> seed = M._to_string()
        >>> N = Manifold('empty')
        >>> N._from_string(seed)
        >>> N == M
        True
        """
        cdef c_Triangulation* c_triangulation = NULL
        if not self.c_triangulation is NULL:
            raise ValueError('The Triangulation must be empty.')
        b_string = to_byte_str(string)
        c_triangulation = read_triangulation_from_string(b_string)
        self.set_c_triangulation(c_triangulation)

    def _to_bytes(self):
        """
        Return a reasonably short byte sequence which encodes the
        combinatorics of this triangulation.
        >>> M = Manifold('m125')
        >>> seed = M._to_bytes()
        >>> len(seed)
        12
        >>> N = Manifold('empty')
        >>> N._from_bytes(seed)
        >>> N == M
        True
        """
        cdef TerseTriangulation* c_terse
        cdef int n, byte
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if False in [ c.is_complete for c in self.cusp_info()]:
            raise ValueError('Byte encoding requires complete cusps.')
        c_terse = tri_to_terse(self.c_triangulation)
        byteseq = bytearray([c_terse.num_tetrahedra])
        byte, bit = 0, 0
        for n from 0 <= n < 2*c_terse.num_tetrahedra:
            if c_terse.glues_to_old_tet[n]:
                byte |= 1 << bit
            bit += 1
            if bit%8 == 0:
                byteseq.append(byte)
                byte, bit = 0, 0
        if bit:
            byteseq.append(byte)
        for n from 0 <= n < 1 + c_terse.num_tetrahedra:
            byteseq.append(c_terse.which_old_tet[n])
        for n from 0 <= n < 1 + c_terse.num_tetrahedra:
            byteseq.append(c_terse.which_gluing[n])
        free_terse_triangulation(c_terse)    
        return bytearray_to_bytes(byteseq)

    def _from_bytes(self, bytestring):
        """
        Fill an empty triangulation from a byte sequence generated by
        _to_bytes.
        >>> M = Manifold('m125')
        >>> seed = M._to_bytes()
        >>> len(seed)
        12
        >>> N = Manifold('empty')
        >>> N._from_bytes(seed)
        >>> N == M
        True
        """
        cdef c_Triangulation* c_triangulation = NULL
        if not self.c_triangulation is NULL:
            raise ValueError('The Triangulation must be empty.')
        c_triangulation = triangulation_from_bytes(bytestring)
        self.set_c_triangulation(c_triangulation)

    def _reindex_cusps(self, permutation):
        cdef int* indices
        cdef int n, num = self.num_cusps()
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if ( len(permutation) != num or set(permutation) != set(range(num)) ):
            raise ValueError('Not a valid permutation')
        indices = <int *>malloc(num*sizeof(int))
        for n in range(num):
            indices[n] = permutation[n]
        reindex_cusps(self.c_triangulation, indices)
        free(indices)
            
    def _get_peripheral_curve_data(self):
        cdef int i, j, k, v, f
        cdef TriangulationData* data
        triangulation_to_data(self.c_triangulation, &data)
        result = []
        for i from 0 <= i < self.num_tetrahedra():
          for j from 0 <= j < 2:       # meridian, longitude 
            for k from 0 <= k < 2:     # righthanded, lefthanded
              row = []
              for v from 0 <= v < 4:
                for f from 0 <= f < 4:
                  row.append(
                     data.tetrahedron_data[i].curve[j][k][v][f]
                     )
              result.append(row)
        free_triangulation_data(data)
        return result

    def isomorphisms_to(self, Triangulation other):
        """
        Returns a complete list of combinatorial isomorphisms between
        the two triangulations:

        >>> M = Manifold('5^2_1')
        >>> N = Manifold('5^2_1')
        >>> N.set_peripheral_curves([[[2,3],[-1,-1]],[[1,1],[0,1]]])
        >>> isoms = M.isomorphisms_to(N)
        >>> isoms[6]
        0 -> 1  1 -> 0
        [ 1 0]  [-1 1]
        [-1 1]  [-3 2]
        Does not extend to link

        Each transformation between cusps is given by a matrix which
        acts on the left.  That is, the two *columns* of the matrix
        give the image of the meridian and longitude respectively.  In
        the above example, the meridian of cusp 0 is sent to the
        meridian of cusp 1.
        """
        cdef IsometryList *isometries = NULL

        if self.c_triangulation is NULL or other.c_triangulation is NULL:
            raise ValueError('Manifolds must be non-empty.')

        compute_cusped_isomorphisms(self.c_triangulation,
                                    other.c_triangulation, 
                                    &isometries,
                                    NULL)

        if isometry_list_size(isometries) == 0:
            result = []
        else:
            result = IsometryListToIsometries(isometries)
        free_isometry_list(isometries)
        return result 
    
    def __dealloc__(self):
        if self.c_triangulation is not NULL:
            free_triangulation(self.c_triangulation)

    def __richcmp__(Triangulation self, Triangulation other, op):
        """
        Two triangulations are equal if they are combinatorially
        isomorphic.  Currently we don't handle the case where there
        are non-trivial Dehn fillings.

        >>> M = Triangulation('m004')
        >>> N = M.copy()
        >>> N == M
        True
        >>> M.dehn_fill( (5,3), 0)
        >>> N == M
        ValueError: Can't compare triangulations of manifolds with Dehn fillings.
        """
        cdef c_Triangulation *c_triangulation1
        cdef c_Triangulation *c_triangulation2
        cdef Boolean answer
        if op != 2:
            return NotImplemented
        if type(self) != type(other):
            return False
        if False in ( self.cusp_info('is_complete') +
                      other.cusp_info('is_complete') ):
            raise ValueError("Can't compare triangulations of manifolds "
                             "with Dehn fillings.")
        if same_triangulation(self.c_triangulation, other.c_triangulation):
            return True
        else:
            return False

    def __repr__(self):
        if self.c_triangulation is NULL:
            return 'Empty Triangulation'
        else:
            repr = self.name()
            for i in range(self.num_cusps()):
                info = self.cusp_info(i)
                if info.is_complete:
                    repr += '(0,0)'
                else:
                    repr += '(%g,%g)'% info['filling']
            return repr
        
    def name(self):
        """
        Return the name of the triangulation.

        >>> M = Triangulation('4_1')
        >>> M.name()
        '4_1'
        """
        if self.c_triangulation is NULL: return
        return SnapPyStr(to_str(get_triangulation_name(self.c_triangulation)))

    def set_name(self, new_name):
        """
        Give the triangulation a new name.

        >>> M = Triangulation('4_1')
        >>> M.set_name('figure-eight-comp')
        >>> M
        figure-eight-comp(0,0)
        """
        b_new_name = to_byte_str(new_name)
        cdef char* c_new_name = b_new_name
        if self.c_triangulation is NULL:
            raise ValueError('The empty triangulation has no name.')
        set_triangulation_name(self.c_triangulation, c_new_name)
    
    def DT_code(self, alpha=False):
        """
        Return the Dowker-Thistlethwaite code supplied when the
        Manifold was instantiated.  This is an immutable attribute,
        intended for use with knot and link exteriors only.  By
        default this returns a list of tuples of even integers.  With
        the flag alpha=True it returns the compressed alphabetical
        form used in the tabulations by Hoste and Thistletwaite.
        """
        if self._DTcode:
            if alpha:
                return self._DTcode
            result = []
            ints = [(64-ord(c)) if ord(c)<96 else (ord(c)-96) 
                    for c in self._DTcode]
            components = ints[1]
            start = 2 + components
            sizes = ints[2:start]
            twox = [x<<1 for x in ints[start:]]
            n = 0
            for size in sizes:
                result.append(tuple(twox[n:n+size]))
                n += size
            return result

    def num_tetrahedra(self):
        """
        Return the number of tetrahedra in the triangulation.

        >>> M = Triangulation('m004')
        >>> M.num_tetrahedra()
        2
        """
        if self.c_triangulation is NULL: return 0
        return SnapPyInt(get_num_tetrahedra(self.c_triangulation))
    
    def dehn_fill(self, filling_data, which_cusp=None):
        """
        Set the Dehn filling coefficients of the cusps.  This can be
        specified in the following ways, where the cusps are numbered
        by 0,1,...,(num_cusps - 1).  

        - Fill cusp 2:

          >>> M = Triangulation('8^4_1')
          >>> M.dehn_fill((2,3), 2)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(1,5)
        
        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          8^4_1(3,0)(1,-4)(2,3)(1,5)

        - When there is only one cusp, there's a shortcut

          >>> N = Triangulation('m004')
          >>> N.dehn_fill( (-3,4) )
          >>> N
          m004(-3,4)
        
        Does not return a new Triangulation.
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        if which_cusp != None:
            try:
                which_cusp = range(self.num_cusps())[which_cusp]
            except IndexError:
                raise IndexError('The specified cusp (%s) does not '
                                 'exist.'%which_cusp)

            meridian, longitude = filling_data
            complete = ( meridian == 0 and longitude == 0)
            set_cusp_info(self.c_triangulation,
                          which_cusp, complete, meridian, longitude)
            self._cache = {}
        else:
            if self.num_cusps() > 1 and len(filling_data) == 2:
                if not hasattr(filling_data, '__getitem__') or not hasattr(filling_data[0], '__getitem__'):
                    raise IndexError('If there is more than one cusp '
                                     'you must specify which one you\n'
                                     'are filling, e.g. M.dehn_fill((2,3),1)')
            if self.num_cusps() == 1 and len(filling_data) == 2:
                self.dehn_fill(filling_data, 0)
                return 
            if len(filling_data) > self.num_cusps():
                raise IndexError('You provided filling data for too '
                                 'many cusps.  There are only %s.'%
                                 self.num_cusps())
            for i, fill in enumerate(filling_data):
                self.dehn_fill(fill, i)

    def cusp_info(self, data_spec=None):
        """
        Returns an info object containing information about the given
        cusp.   Usage:

        >>> M = Triangulation('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0)
        >>> c = M.cusp_info(1)
        >>> c.is_complete
        False
        >>> c.keys()
        ['index', 'filling', 'is_complete', 'topology']

        You can get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : torus cusp, not filled,
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('is_complete')
        [True, False, False]
        """
        cdef int cusp_index
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef double m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_accuracy, current_shape_accuracy,
        cdef Complex initial_modulus, current_modulus
        cdef int meridian_accuracy, longitude_accuracy, singularity_index, accuracy
        cdef Complex c_meridian, c_longitude, c_core_length

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        if data_spec == None:
            return ListOnePerLine([self.cusp_info(i)
                                   for i in range(self.num_cusps())])
        if type(data_spec) == type(''):
            return [c[data_spec] for c in self.cusp_info()]
        try:
            cusp_index = range(self.num_cusps())[data_spec]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%data_spec)

        get_cusp_info(self.c_triangulation, cusp_index,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_accuracy, &current_shape_accuracy,
                      &initial_modulus, &current_modulus)
        info = {
           'index' : cusp_index,
           'topology' : CuspTopology[topology],
           'is_complete' : B2B(is_complete),
           'filling' : (m, l)
           }

        #If there's a hyperbolic structure, there more information to
        #pass on.
        if hasattr(self, 'tetrahedra_shapes'):
            get_holonomy(self.c_triangulation, cusp_index,
                         &c_meridian, &c_longitude,
                         &meridian_accuracy, &longitude_accuracy)
            shape = SnapPyComplex(C2C(current_shape))
            shape.accuracy = current_shape_accuracy
            meridian = SnapPyComplex(C2C(c_meridian))
            meridian.accuracy = meridian_accuracy
            longitude = SnapPyComplex(C2C(c_longitude))
            longitude.accuracy = longitude_accuracy
            info.update({
                'shape':shape,
                'shape_accuracy':current_shape_accuracy,
                'modulus':SnapPyComplex(C2C(current_modulus)),
                'holonomies':(meridian, longitude),
                'holonomy_accuracy':min(meridian_accuracy,longitude_accuracy)
                })

            core_geodesic(self.c_triangulation, cusp_index,
                          &singularity_index, &c_core_length, &accuracy)
            
            if singularity_index != 0:
                core_length = SnapPyComplex(C2C(c_core_length))
                core_length.accuracy = accuracy
                info.update({
                    'core_length':core_length,
                    'singularity_index':singularity_index
                    })
                
        return CuspInfo(**info)

    def reverse_orientation(self):
        """
        Reverses the orientation of the Triangulation, presuming that
        it is orientable.

        >>> M = Manifold('m015')
        >>> cs = M.chern_simons()
        >>> M.reverse_orientation()
        >>> cs + M.chern_simons()
        0.0
        """
        if not self.is_orientable():
            raise ValueError("The Manifold is not orientable, so its "
                             "orientation can't be reversed.")
        reorient(self.c_triangulation)
        self._cache = {}
            
    def filled_triangulation(self, cusps_to_fill='all'):
        """
        Return a new manifold where the specified cusps have been
        permanently filled in.  Examples:

        Filling all the cusps:
        
        >>> M = Triangulation('m125(1,2)(3,4)')
        >>> N = M.filled_triangulation()
        >>> N.num_cusps()
        0

        Filling cusps 0 and 2 :

        >>> M = Triangulation('v3227(1,2)(3,4)(5,6)')
        >>> M.filled_triangulation([0,2])
        v3227_filled(3,4)
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        n = self.num_cusps()
        if cusps_to_fill == 'all':
            cusps_to_fill = [c for c in range(n) if cusp_is_fillable(self.c_triangulation, c)]
                
        if False in [(c in range(n)) for c in cusps_to_fill]:
            raise IndexError('The specified indices to be filled are beyond '
                             'the actual number of cusps.')
        if 0 in [cusp_is_fillable(self.c_triangulation, c)
                 for c in cusps_to_fill]:
            raise IndexError('To permanently fill a cusp, the Dehn '
                             'filling coefficients must be relatively '
                             'prime integers.')

        cdef c_Triangulation* c_filled_tri = NULL
        cdef Triangulation filled_tri
        cdef Boolean *fill_cusp_spec = NULL
        
        fill_cusp_spec = <Boolean*>malloc(n*sizeof(Boolean))
        for i in range(n):
            fill_cusp_spec[i] = 1 if i in cusps_to_fill else 0
        fill_all = 1 if not False in [i in cusps_to_fill
                                      for i in range(n)] else 0
        c_filled_tri = fill_cusps(self.c_triangulation,
                                  fill_cusp_spec, '', fill_all)
        free(fill_cusp_spec)
        filled_tri = Triangulation('empty')
        filled_tri.set_c_triangulation(c_filled_tri)
        filled_tri.set_name(self.name() + '_filled')
        return filled_tri
        
    def edge_valences(self):
        """
        Returns a dictionary whose keys are the valences of the edges
        in the triangulation, and the value associated to a key is the
        number of edges of that valence.

        >>> M = Triangulation('v3227')
        >>> M.edge_valences()
        {10: 1, 4: 1, 5: 2, 6: 3}
        """
        cdef int c, v = 1
        ans = {}
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        while get_num_edge_classes(self.c_triangulation, v, 1) > 0:
            c = get_num_edge_classes(self.c_triangulation, v, 0)
            if c > 0:
                ans[v] = c
            v += 1
        return ans
        
    def gluing_equations_pgl(self, N=2, equation_type='all'):

        """
        This method returns a matrix of exponents for gluing equations of
        cross ratios of PGL(N,C) representations where N (default 2) can be
        specified. See 

	Garoufalidis, Goerner, Zickert: 
	Gluing Equations for PGL(n,C)-Representations of 3-Manifolds.
        
        and http://www.unhyperbolic.org/ptolemy.html .
        
        In the default mode, the function returns an equation object containing
        a matrix with rows of the form 

        a b c d ... f ...
        
        which means

         z_0001_0^a * zp_0001_0^b * zpp_0001_0^c * z_0010_0^d * ... * z_0001_1^a = 1

        where z's denote cross ratios at the edge (0,1), zp's (usually denoted
	z') at (0,2) and zpp's (usually denoted z'') at (1,2). 
        See kernel_code/edge_classes.c for a detailed account
        of the convention used. The first index denotes the Ptolemy index 
        (integral point), the second index the tetrahedron.

	Note that the SnapPy conventions are slightly different from the above
	paper: 1. the paper used upper indicies instead of ' and ''. 2. the
        paper is refering to the z' and z'' notation, but switches z' and z''.
    
        The value of equation_type can be (default is 'all'):

        * 'all'               # list all gluing equations
        * 'non_peripheral'    # list non-peripheral equations

          * 'edge'            # list edge gluing equations
          * 'face'            # list face gluing equations
          * 'internal'        # list internal gluing equations

        * 'peripheral'        # list cusp gluing equations

          * 'meridian'        # list cusp gluing equations for meridians
          * 'longitude'       # list cusp gluing equations for longitudes
    
        >>> M = Triangulation('m004')
        >>> M.gluing_equations_pgl(N=2).explain_columns
        ['z_0000_0', 'zp_0000_0', 'zpp_0000_0', 'z_0000_1', 'zp_0000_1', 'zpp_0000_1']
        >>> M.gluing_equations_pgl(N=2).explain_rows
        ['edge_0_0', 'edge_0_1', 'meridian_0_0', 'longitude_0_0']
        >>> M.gluing_equations_pgl(N=2).matrix    
        matrix([[ 2,  1,  0,  1,  0,  2],
                [ 0,  1,  2,  1,  2,  0],
                [ 1,  0,  0,  0, -1,  0],
                [ 0,  0,  0,  0, -2,  2]])
        """


        cdef Integer_matrix_with_explanations c_matrix

        if N < 2 or N > 15:
            raise ValueError('N has to be 2...15')

        if not equation_type in ['all', 
                                     'non_peripheral', 
                                         'edge', 'face', 'internal',
                                     'peripheral',
                                         'longitude', 'meridian']:
            raise ValueError('Wrong equation_type')
        
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        equations = []
        explain_rows = []
        explain_cols = []

        if equation_type in [ 'all', 'non_peripheral', 'edge']:

            # Add edge equations
            get_edge_gluing_equations_pgl(self.c_triangulation,
                                          &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in [ 'all', 'non_peripheral', 'face']:

            # Add face equations
            get_face_gluing_equations_pgl(self.c_triangulation,
                                          &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in [ 'all', 'non_peripheral',  'internal']:

            # Add internal equations
            get_internal_gluing_equations_pgl(self.c_triangulation,
                                              &c_matrix, N)
            eqns, r, explain_cols = convert_and_free_integer_matrix(c_matrix)
            equations += eqns
            explain_rows += r

        if equation_type in ['all', 'peripheral', 'longitude', 'meridian']:

            # Add peripheral equations
            
            for i in range(self.num_cusps()):
                cusp_info = self.cusp_info(i)

                to_do = [] # keep a todo list where we add (meridian,longitude)
                           # pairs to process later

                if cusp_info.is_complete:

                    # Add meridians for complete cusps
                    if equation_type in [ 'meridian', 'peripheral', 'all']:
                        to_do += [ (1,0) ]
                        explain_rows += [
                            "meridian_%d_%d" % (j, i) for j in range(N-1) ]

                    # Add longitudes for complete cusps
                    if equation_type in [ 'longitude', 'peripheral', 'all']:
                        to_do += [ (0,1) ]
                        explain_rows += [
                            "longitude_%d_%d" % (j, i) for j in range(N-1) ]

                else:

                    # Add Dehn-filling for incomplete cusp
                    to_do += [ cusp_info.filling ]
                    explain_rows += [
                        "filling_%d_%d" % (j, i) for j in range(N-1) ]

                # process the todo list
                for (m, l) in to_do:

                    get_cusp_equations_pgl(self.c_triangulation,
                                           &c_matrix,
                                           N, i, m, l)

                    eqns, r, explain_cols = (
                        convert_and_free_integer_matrix(c_matrix))
                    equations += eqns

        if equations == []: # cover cases N = 2, 3 and equation_type = 'internal'
            return None

        return NeumannZagierTypeEquations(matrix(equations),
                                          explain_rows,
                                          explain_cols)

    def _ptolemy_equations_identified_face_classes(self):
        """
        This function returns an identification structure where s_f_t gets 
        identified with -s_g_u if face f of tetrahedron t is glued to face g of
        tetrahedron u. 

        We can represent a 2-cohomology class H^2(M,boundary M) by denoting by
        s_f_t the value the 2-cohomology class takes on the face f of 
        tetrahedron t with the orientation being the one induced from the
        orientation of the tetrahedron.
        Because a face class of the triangulation has two representatives
        (tet_index, face_index) and the gluing is orientation-reversing on the
        face, one s will be the negative of another s.
        """
 
        cdef Identification_of_variables c_vars

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_identified_face_classes(
            self.c_triangulation, &c_vars)

        return convert_and_free_identification_of_variables(c_vars)        
                
    def _ptolemy_equations_identified_coordinates(self, N):

        """
        Ptolemy coordinates that need to be identified for the given
        triangulation when computing pSL(N,C) representations. 
	"""
 
        cdef Identification_of_variables c_vars

        if N < 2 or N > 15:
            raise ValueError('N has to be 2...15')

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_identified_coordinates(
            self.c_triangulation, &c_vars, N)

        return convert_and_free_identification_of_variables(c_vars)

    def _ptolemy_equations_action_by_decoration_change(self, int N):
        """
        We can change a decoration by multiplying a coset of a cusp by a
        diagonal matrix. Let's let a diagonal matrix SL(n,C) with diagonal
        entries 1 1 ... z 1 ... 1 1/z (z at positon j) act on cusp i. It
        changes some Ptolemy coordinate c_p_t by some power z^n.
        This is expressed in the following matrix as the entry in row
        labeld c_p_t and the column labeled diagonal_entry_j_on_cusp_i.
        """
        
        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_action_by_decoration_change(
            self.c_triangulation,N, &c_matrix)
        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def _ptolemy_equations_boundary_map_3(self):
        """
        Boundary map C_3 -> C_2 in cellular homology represented as matrix

        The following map represents the boundary map in the cellular chain
        complex when representing a linear map as a matrix m acting on a column
        vector v by left-multiplication m * v. With right-multiplication acting
        on row vectors, the matrix represents the coboundary map in the cochain
        complex. 
        
        The basis for C_3 are just the oriented tetrahedra of the triangulation.
        The basis for C_2 are the face classes, see 
        _ptolemy_equations_identified_face_classes.
	"""
        
        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_boundary_map_3(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols
                                             
    def _ptolemy_equations_boundary_map_2(self):
        """
        Boundary map C_2 -> C_1 in cellular homology represented as matrix.

        Also see _ptolemy_equations_boundary_map_3.
	"""
        
        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_boundary_map_2(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def ptolemy_obstruction_classes(self):

        """
        Generates a list of obstruction cocycles representing each class in
        H^2(M,bd M; Z/2) suitable as argument for get_ptolemy_variety.
        The first element in the list is always the trivial obstruction class.

        See Definition 1.7 of
        Garoufalidis, Thurston, Zickert
        The Complex Volume of SL(n,C)-Representations of 3-Manifolds
        http://arxiv.org/abs/1111.2828

        s_f_t takes values +/-1 and is the value of evaluating the cocycle on
        face f of tetrahedron t.

        === Examples ===

        Get the obstruction classes for 4_1:
    
        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_obstruction_classes()
        
        There are two such clases for 4_1:
    
        >>> len(c)
        2
    
        Print the non-trivial obstruction class:

        >>> c[1]
        PtolemyObstructionClass(s_0_0 + 1, s_1_0 - 1, s_2_0 - 1, s_3_0 + 1, s_0_0 - s_0_1, s_1_0 - s_3_1, s_2_0 - s_2_1, s_3_0 - s_1_1)

        Construct Ptolemy variety for non-trivial obstruction class and N = 2:
    
        >>> p = M.ptolemy_variety(2, obstruction_class = c[1])

        Short cut for the above code:
    
        >>> p = M.ptolemy_variety(2, obstruction_class = 1)

        Obstruction class only makes sense for even N:
    
        >>> p = M.ptolemy_variety(3, obstruction_class = c[1])
        Traceback (most recent call last):
        ...
        AssertionError: PtolemyObstructionClass only makes sense for even N


        Hence, we get only one variety if we ask for all obstruction classes:
    
        >>> len(M.ptolemy_variety(3, obstruction_class = 'all'))
        1
        """

        return ptolemyManifoldMethods.get_ptolemy_obstruction_classes(self)

    def ptolemy_variety(self, N, obstruction_class = None, simplify = True):

        """		      
        Generates Ptolemy variety as described in
        (1) Garoufalidis, Thurston, Zickert
        The Complex Volume of SL(n,C)-Representations of 3-Manifolds
        http://arxiv.org/abs/1111.2828
        
        (2) Garoufalidis, Goerner, Zickert:
        Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
        http://arxiv.org/abs/1207.6711
        
        The variety can be exported to magma or sage and solved there. The
        solutions can be processed to compute invariants. See below.
        
        === Arguments ===
        
        N --- which SL(N,C) we want the variety.
        
        obstruction_class --- class from Definiton 1.7 of (1).
        None for trivial class or a value returned from get_ptolemy_obstruction_classes.
        Short cuts: obstruction_class = 'all' returns a list of Ptolemy varieties
        for each obstruction. For easier iteration, can set obstruction_class to 
        an integer.
        
        simplify --- boolean to indicate whether to simplify the equations which
        significantly reduces the number of variables.
        Simplifying means that several identified Ptolemy coordinates x = y = z = ...
        are eliminated instead of adding relations x - y = 0, y - z = 0, ...
        
        === Examples for 4_1 ===
        
        >>> M = Manifold("4_1")
        
        Get the varieties for all obstruction classes at once (use
        help(varieties[0]) for more information):
        
        >>> varieties = M.ptolemy_variety(2, obstruction_class = "all")
        
        Print the equations of the variety for the non-trivial class:
        
        >>> for eqn in varieties[1].equations:
        ...     print "    ", eqn
             1 - c_0101_0 + c_0101_0^2
             - 1 + c_0101_0 - c_0101_0^2
        
        Generate a magma file to compute Primary Decomposition for N = 3:
        
        >>> p = M.ptolemy_variety(3)
        >>> print p.to_magma()          #doctest: +ELLIPSIS
        P<t, c_0012_1, c_0102_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 7);
        I := ideal<P |
        c_0102_0 - c_0102_0 * c_1011_0 + c_1101_0,
        ...
        
        === If you have a magma installation ===
        
        Call p.compute_solutions() to automatically call magma on the above output
        and produce exact solutions!!!
        
        >>> try:
        ...     sols = p.compute_solutions()
        ... except:
        ...     sols = None     # magma failed, use precomputed output instead
        
        === If you do not have a magma installation ===
        
        Load a precomputed example from magma which is provided with the package:
        
        >>> from ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> print _magma_output_for_4_1__sl3      #doctest: +ELLIPSIS
        <BLANKLINE>
        ==TRIANGULATION=BEGINS==
        % Triangulation
        4_1
        ...
        
        Parse the file and produce solutions:
        
        >>> if sols is None:    # calling magma failed, so use precomputed example
        ...     sols = solutions_from_magma(_magma_output_for_4_1__sl3)
            
        === Continue here whether you have or do not have magma ===
        
        Pick the first solution of the three different solutions (up to Galois
        conjugates):
        
        >>> len(sols)
        3
        >>> solution = sols[0]
        
        Read the exact value for c_1020_0 (help(solution) for more information
        on how to compute cross ratios, volumes and other invariants):
        
        >>> solution['c_1020_0']
        Mod(1/2*x - 1, x^2 - x + 2)
        
        Example of simplified vs non-simplified variety for N = 4:
        
        >>> simplified = M.ptolemy_variety(4, obstruction_class = 1)
        >>> full = M.ptolemy_variety(4, obstruction_class = 1, simplify = False)
        >>> len(simplified.variables), len(full.variables)
        (17, 70)
        >>> len(simplified.equations), len(full.equations)
        (20, 79)
        """
        
        return ptolemyManifoldMethods.get_ptolemy_variety(
            self, N, obstruction_class, simplify)

    def gluing_equations(self,form='log'):
        """
        In the default mode, this function returns a matrix with rows
        of the form

                  a b c  d e f  ...

        which means

            a*log(z0) + b*log(1/(1-z0)) + c*log((z0-1)/z0) + d*log(z1) +... = 2 pi i

        for an edge equation, and (same) = 0 for a cusp equation.
        Here, the cusp equations come at the bottom of the matrix, and
        are listed in the form: meridian of cusp 0, longitude of cusp
        0, meridian of cusp 1, longitude of cusp 1,...

        In terms of the tetrahedra, a is the invariant of the edge
        (2,3), b the invariant of the edge (0,2) and c is the
        invariant of the edge (1,2).  See kernel_code/edge_classes.c
        for a detailed account of the convention used.  

        If the optional argument form='rect' is given, then this
        function returns a list of tuples of the form:

           ( [a0, a1,..,a_n], [b_0, b_1,...,b_n], c)

        where this corresponds to the equation

           z0^a0 (1 - z0)^b0 z1^a1(1 - z1)^b1 ...  = c

        where c = 1 or -1.

        >>> M = Triangulation('m004(2,3)')
        >>> M.gluing_equations()
        matrix([[ 2,  1,  0,  1,  0,  2],
                [ 0,  1,  2,  1,  2,  0],
                [ 2,  0,  0,  0, -8,  6]])
        >>> M.gluing_equations(form='rect')
        [([2, -1], [-1, 2], 1), ([-2, 1], [1, -2], 1), ([2, -6], [0, 14], 1)]
        """
        
        cdef int **c_eqns
        cdef int num_rows, num_cols
        cdef int* eqn

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        c_eqns = get_gluing_equations(self.c_triangulation,
                                      &num_rows, &num_cols)
        eqns = [ [ c_eqns[i][j] for j in range(num_cols) ]
                 for i in range(num_rows) ]
        free_gluing_equations(c_eqns, num_rows)

        for i in range(self.num_cusps()):
            cusp_info = self.cusp_info(i)
            if cusp_info.is_complete:
                to_do = [(1,0), (0,1)]
            else:
                to_do = [cusp_info.filling]
            for (m, l) in to_do:
                eqn = get_cusp_equation(self.c_triangulation,
                                        i, m, l, &num_rows)
                eqns.append([eqn[j] for j in range(num_rows)])
                free_cusp_equation(eqn)

        if form == 'log':
            return matrix(eqns)

        if form != 'rect':
            raise ValueError("Equations are available in 'log' and "
                             "'rect' forms only.")
        rect = []
        for row in eqns:
            n = self.num_tetrahedra()
            a, b = [0,]*n, [0,]*n
            c = 1
            for j in range(n):
                r = row[3*j + 2]
                a[j] = row[3*j] - r
                b[j] = -row[3*j + 1] + r
                c *= -1 if r % 2 else 1
            rect.append( (a, b, c) )
        return rect
                                                     
    def homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.
        
        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z

        """
        if 'homology' in self._cache.keys():
            return self._cache['homology']
        
        cdef c_AbelianGroup *H
        cdef RelationMatrix R
        cdef int m, n

        if self.c_triangulation is NULL:
            return AbelianGroup()
        H = homology(self.c_triangulation)
        if H != NULL:
            coefficient_list = []
            compress_abelian_group(H)
            for n from 0 <= n < H.num_torsion_coefficients:
                coefficient_list.append(H.torsion_coefficients[n])
            free_abelian_group(H)
            result = AbelianGroup(elementary_divisors=coefficient_list)
        else:
            homology_presentation(self.c_triangulation, &R)
            relations = []
            if R.relations != NULL:
                if R.num_rows == 0:
                    relations = [0,] * R.num_columns
                else:   
                    for m from 0 <= m < R.num_rows:
                        row = []
                        for n from 0 <= n < R.num_columns:
                            row.append(R.relations[m][n])
                        relations.append(row)
                    free_relations(&R)
            else:
                raise ValueError("The SnapPea kernel couldn't compute "
                                 "the homology presentation matrix")
            result = AbelianGroup(relations)
        self._cache['homology'] = result
        return result

    def fundamental_group(self,
                          simplify_presentation = True,
                          fillings_may_affect_generators = True,
                          minimize_number_of_generators = True):
        """
        Returns a FundamentalGroup object representing the fundamental
        group of the manifold.  If integer Dehn surgery parameters
        have been set, then the corresponding peripheral elements are
        killed.

        >>> M = Triangulation('m004')
        >>> G = M.fundamental_group()
        >>> G
        Generators:
           a,b
        Relators:
           aaabABBAb
        >>> G.peripheral_curves()
        [('ab', 'aBAbABab')]
        
        There are three optional arguments all of which default to True:

        - simplify_presentation
        - fillings_may_affect_generators
        - minimize_number_of_generators

        >>> M.fundamental_group(False, False, False)
        Generators:
           a,b,c
        Relators:
           CbAcB
           BacA
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        name_mangled = 'fundamental_group-%s-%s-%s' %\
                       (simplify_presentation,
                        fillings_may_affect_generators,
                        minimize_number_of_generators)
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = FundamentalGroup(self, simplify_presentation, fillings_may_affect_generators, minimize_number_of_generators)
        return self._cache[name_mangled]
    
    def cover(self, permutation_rep):
        """
        Returns a Triangulation representing the finite cover
        specified by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the simplified presentation of the
        fundamental group.  Each permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.

        >>> M = Triangulation('m004')
        >>> N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])
        >>> N0.homology()
        Z + Z + Z
        
        Within Sage the permutations can also be of type
        PermutationGroupElement, in which case they act on the set
        range(1, d + 1).  Or, you can specify a GAP or Magma subgroup
        of the fundamental group.  For examples, see the docstring for
        Manifold.cover
        """
        cdef RepresentationIntoSn* c_representation
        cdef c_Triangulation* c_triangulation
        cdef Triangulation cover

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        # For Sage, we need to check if we have been given some
        # alternate inputs
        
        if _within_sage:
            if is_GapElement(permutation_rep):
                if permutation_rep.IsSubgroupFpGroup():
                    GG = gap(self.fundamental_group())
                    coset_action = GG.FactorCosetAction(permutation_rep)
                    gap_image_gens = coset_action.Image().GeneratorsOfGroup()
                    Q = PermutationGroup(gap_image_gens)
                    return self.cover([Q(g) for g in gap_image_gens ])
                elif permutation_rep.IsToPermGroupHomomorphismByImages():
                    f = permutation_rep
                    return self.cover(f.PreImage(f.Image().Stabilizer(1)))
                    
            elif is_MagmaElement(permutation_rep):
                input_type = repr(permutation_rep.Type())
                if input_type == 'GrpFP':
                    GG = magma(self.fundamental_group())
                    f = GG.CosetAction(permutation_rep)
                elif input_type == 'HomGrp':
                    f = permutation_rep
                    if not repr(f.Image().Type()) == 'GrpPerm':
                        raise TypeError('The homomorphism image is not '
                                        'a permutation group.')
                else:
                    raise TypeError('That Magma type not recognized.')
                
                magma.eval("""\
                     FormatHomForSnapPea := function(f)
                         subone := function(L)   return [x - 1 : x in L]; end function;
                         return [subone(Eltseq(f(g))) : g in Generators(Domain(f))];
                       end function;""")
                permutation_rep = f.FormatHomForSnapPea().sage()

            # Not a useful GAP or MAGMA object, so let's try.  
            elif not False in [is_PermutationGroupElement(p)
                               for p in permutation_rep]:
                permutation_rep = [ [x - 1 for x in perm.list()]
                                   for perm in permutation_rep ]

        G = self.fundamental_group()
        c_representation = self.build_rep_into_Sn(permutation_rep)
        degree = len(permutation_rep[0])
        c_triangulation = construct_cover(self.c_triangulation,
                                          c_representation,
                                          degree)
        cover = Triangulation('empty')
        cover.set_c_triangulation(c_triangulation)
        cover.set_name(self.name() +'~')
        free_representation(c_representation,
                            G.num_original_generators(),
                            self.num_cusps())
        return cover

    def covers(self, degree, method=None, cover_type='all'):
        """
        Returns a list of Triangulations corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Triangulation('m003')
        >>> covers = M.covers(4)
        >>> [(N, N.homology()) for N in covers]
        [(m003~irr~0(0,0)(0,0), Z/5 + Z + Z), (m003~cyc~1(0,0), Z/3 + Z/15 + Z)]

        You can also look just at cyclic covers, which is much faster.

        >>> covers = M.covers(4, cover_type='cyclic')
        >>> [(N, N.homology()) for N in covers]
        [(m003~cyc~0(0,0), Z/3 + Z/15 + Z)]

        If you are using Sage, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional argument

        method = 'gap'

        If in addition you have Magma installed, you can use it to do
        the heavy-lifting by specifying method = 'magma'.
        """
        cdef RepresentationList* reps
        cdef RepresentationIntoSn* rep
        cdef c_Triangulation* cover
        cdef Triangulation T

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
                
        if cover_type == 'cyclic':
            method = None
            
        if method:
            if not _within_sage:
                raise RuntimeError('Only the default method of finding '
                                   'subgroups is available, as you are '
                                   'not using Sage.')
            if method == 'gap':
                G = gap(self.fundamental_group())
                return [self.cover(H)
                        for H in G.LowIndexSubgroupsFpGroup(degree)
                        if G.Index(H) == degree]
            if method == 'magma':
                G = magma(self.fundamental_group())
                return [self.cover(H)
                        for H in G.LowIndexSubgroups('<%d, %d>' %
                                                     (degree, degree))]

        if cover_type == 'all':
            reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Sn)
        elif cover_type == 'cyclic':
            reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Zn)
        else:
            raise ValueError("Supported cover_types are 'all' "
                             "and 'cyclic'.")

        covers = []
        rep = reps.list
        cover_count = 0
        while rep != NULL:
            cover = construct_cover(self.c_triangulation,
                                    rep,
                                    reps.num_sheets)
            T = Triangulation('empty')
            T.set_c_triangulation(cover)
            covers.append(T)
            cover_types = {1:"irr", 2:"reg", 3:"cyc"}
            T.set_name(self.name() + "~" + cover_types[rep.covering_type] + '~%d' % cover_count)
            cover_count += 1
            rep = rep.next
            
        free_representation_list(reps)
        return covers

    cdef RepresentationIntoSn *build_rep_into_Sn(self, perm_list) except ? NULL:
        """
        Build a SnapPea RepresentationIntoSn from a list of
        permutations, one for each generator of the simplified
        fundamental group.  A permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.  The representation constructed here is given in terms
        of the geometric generators, for use in construcing a covering
        space.
        """
        cdef c_Triangulation* cover
        cdef c_Triangulation* c_triangulation
        cdef c_GroupPresentation *c_group_presentation
        cdef RepresentationIntoSn* c_representation
        cdef RepresentationIntoSn* c_repn_in_original_gens = NULL
        cdef int i, j
        cdef num_generators, num_relators, num_orig_gens, num_cusps
        cdef int** c_original_generators
        cdef int** c_relators
        cdef int** c_meridians
        cdef int** c_longitudes

        degree = len(perm_list[0])

        # Sanity check
        S = set(range(degree))
        for permutation in perm_list:
            if set(permutation) != S:
                raise ValueError('The permutation list is invalid.')

        # Initialize
        num_cusps = self.num_cusps()
        c_triangulation = self.c_triangulation
        c_group_presentation = fundamental_group(c_triangulation,
                                             True, True, True)
        num_generators = fg_get_num_generators(c_group_presentation)
        num_relators = fg_get_num_relations(c_group_presentation)
        num_orig_gens = fg_get_num_orig_gens(c_group_presentation)

        # Allocate a whole bunch of memory, SnapPea and malloc.
        c_representation = initialize_new_representation(
            num_orig_gens,
            degree,
            num_cusps)
        for i from 0 <= i < num_generators:
            for j from 0 <= j < degree:
                c_representation.image[i][j] = perm_list[i][j]
        c_original_generators = <int**>malloc(num_orig_gens*sizeof(int*))
        for i from  0 <= i < num_orig_gens:
            c_original_generators[i] = fg_get_original_generator(
                c_group_presentation, i)
        c_relators = <int**>malloc(num_relators*sizeof(int*))
        for i from  0 <= i < num_relators:
            c_relators[i] = fg_get_relation(c_group_presentation, i)
        c_meridians = <int**>malloc(num_cusps*sizeof(int*))
        c_longitudes = <int**>malloc(num_cusps*sizeof(int*))
        for i from 0 <= i < num_cusps:
            c_meridians[i] = fg_get_meridian(c_group_presentation, i)
            c_longitudes[i] = fg_get_longitude(c_group_presentation, i)
        # Whew!

        failed = False
        if (candidateSn_is_valid(c_representation.image, 
                                 degree, c_relators, num_relators) and
            candidateSn_is_transitive(c_representation.image,
                                      num_generators, degree) ):
            c_repn_in_original_gens = convert_candidateSn_to_original_generators(
                c_representation.image,
                degree,
                num_orig_gens,
                c_original_generators,
                c_triangulation,
                c_meridians,
                c_longitudes)
        else:
            message = 'Invalid permutation data.'
            failed = True
        if c_repn_in_original_gens == NULL:
            message = 'Failed to construct permutation representation.'
            failed = True
    
        # Now free all that memory
        for i from 0 <= i < num_cusps:
            fg_free_relation(c_meridians[i])
            fg_free_relation(c_longitudes[i])
        free(c_meridians)
        free(c_longitudes)
        for i from 0 <= i < num_relators:
            fg_free_relation(c_relators[i])
        free(c_relators)
        for i from 0 <= i < num_orig_gens:
            fg_free_relation(c_original_generators[i])
        free(c_original_generators)
        free_representation(c_representation, num_generators, num_cusps)
        # Free at last!

        if failed:
            raise RuntimeError(message)
        return c_repn_in_original_gens

###  SnapPeaX not yet available in a usable form.  
##
##    def copy_to_SnapPeaX(manifold):
##        """
##        Copies the manifold over into the OS X SnapPea GUI.
##        """
##        file_name = tempfile.mktemp() + ".tri"
##        manifold.save(file_name)
##        script = """\
##        set f to POSIX file "%s"
##        tell application "%s"
##        open f
##        end tell""" % (file_name, SnapPeaX_name)
##        execute_applescript(script)
##        activate_SnapPeaX()

# Manifolds

cdef class Manifold(Triangulation):
    """
    A Manifold is a Triangulation together with a geometric structure.
    That is, a Manifold is an ideal triangulation of the interior of a
    compact 3-manifold with torus boundary, where each tetrahedron has
    has been assigned the geometry of an ideal tetrahedron in
    hyperbolic 3-space.  A Dehn-filling can be specified for each
    boundary component, allowing the description of closed 3-manifolds
    and some orbifolds.   Here's a quick example:

    >>> M = Manifold('9_42')
    >>> M.volume()
    4.05686022424
    >>> M.cusp_info('shape')
    [(-4.2789363159+1.9572867975j)]

    A Manifold can be specified in a number of ways, e.g.

    - Manifold('9_42') : The complement of the knot 9_42 in S^3.
    - Manifold('m125(1,2)(4,5)') : The SnapPea census manifold m125
       where the first cusp has Dehn filling (1,2) and the second cusp has
       filling (4,5).
    - Manifold() : Opens a link editor window where can you
       specify a link complement.
    
    In general, the specification can be from among the below, with
    information on Dehn fillings added.

    - SnapPea cusped census manifolds: e.g. 'm123', 's123', 'v123'.

    - Link complements:
       + Rolfsen's table: e.g. '4_1', '04_1', '5^2_6', '6_4^7', 'L20935', 'l104001'.
       + Hoste-Thistlethwaite Knotscape table:  e.g. '11a17' or '12n345'
       + Callahan-Dean-Weeks-Champanerkar-Kofman-Patterson knots: e.g. 'K6_21'.
       + Dowker-Thistlethwaite code: e.g. 'DT[6,8,2,4]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a_0, B_1])' or 'Splitting(S_{1,0}, [b_1, A_0], [a_0,B_1])'

      See the help for the 'twister.twister' function for more.  

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """
    def __init__(self, spec=None, DTcode=None):
        if self.c_triangulation != NULL:
            find_complete_hyperbolic_structure(self.c_triangulation)
            do_Dehn_filling(self.c_triangulation)

    def canonize(self):
        """
        Change the triangulation to an arbitrary retriangulation of
        the canonical cell decomposition.

        >>> M = Manifold('m007')
        >>> M.num_tetrahedra()
        3
        >>> M.canonize()
        >>> M.num_tetrahedra()
        4

        Note: do to rounding error, it is possible that this is not
        actually the canonical triangulation.          
        """
        cdef c_FuncResult result
        result = proto_canonize(self.c_triangulation)
        if FuncResult[result] != 'func_OK':
            raise RuntimeError('SnapPea failed to find the canonical '
                               'triangulation.')

    def _canonical_cells_are_tetrahedra(self):
        """
        Returns True if and only if the canonical
        cell decomposition is a triangulation.
        """
       
        cdef Manifold M
        M = self.copy()
        M.canonize()
        canonical_retriangulation(M.c_triangulation)
        return not B2B(mark_fake_cusps(M.c_triangulation))

    def copy(self):
        """
        Returns a copy of the manifold

        >>> M = Manifold('m015')
        >>> N = M.copy()
        """
        return Manifold_from_Triangulation(Triangulation.copy(self))

    def cusp_neighborhood(self):
        """
        Returns information about the cusp neighborhoods of the
        manifold, in the form of data about the corresponding horoball
        diagrams in hyperbolic 3-space.
        
        >>> M = Manifold('s000')
        >>> CN = M.cusp_neighborhood()
        >>> CN.volume()
        0.3247595264191645
        >>> len(CN.horoballs(0.01))
        178
        >>> CN.view()  # Opens 3-d picture of the horoballs 
        """
        return CuspNeighborhood(self)

    def dirichlet_domain(self,  vertex_epsilon=10.0**-8,
                      displacement = [0.0, 0.0, 0.0],
                      centroid_at_origin=True,
                      maximize_injectivity_radius=True):
        """
        Returns a DirichletDomain object representing a Dirichlet
        domain of the hyperbolic manifold, typically centered at a
        point which is a local maximum of injectivity radius.  It will
        have ideal vertices if the manifold is not closed.

        >>> M = Manifold('m015')
        >>> D = M.dirichlet_domain()
        >>> D
        32 finite vertices, 2 ideal vertices; 54 edges; 22 faces
        >>> D.view()   #Shows 3d-graphic of the DirichletDomain.  
        
        Other options can be provided to customize the computation;
        the default choices are shown below:

        >>> M.dirichlet_domain(vertex_epsilon=10.0**-8,  displacement = [0.0, 0.0, 0.0],
        ... centroid_at_origin=True, maximize_injectivity_radius=True)
        32 finite vertices, 2 ideal vertices; 54 edges; 22 faces
        """
        return DirichletDomain(self, vertex_epsilon, displacement, centroid_at_origin, maximize_injectivity_radius)

    def filled_triangulation(self, cusps_to_fill='all'):
        """
        Return a new Manifold where the specified cusps have been
        permanently filled in.  

        Filling all the cusps results in a Triangulation rather
        than a Manifold, since SnapPea can't deal with hyperbolic
        structures when there are no cusps. 

        Examples:
        
        >>> M = Manifold('m125(1,2)(3,4)')
        >>> N = M.filled_triangulation()
        >>> N.num_cusps()
        0
        >>> type(N) == Triangulation
        True

        Filling cusps 0 and 2 :

        >>> M = Manifold('v3227(1,2)(3,4)(5,6)')
        >>> M.filled_triangulation([0,2])
        v3227_filled(3,4)
        """
        filled = Triangulation.filled_triangulation(self, cusps_to_fill)
        if filled.num_cusps() == 0:
            return filled
        return Manifold_from_Triangulation(filled)
       
    def fundamental_group(self,
                   simplify_presentation = True,
                   fillings_may_affect_generators = True,
                   minimize_number_of_generators = True):
        """
        Return a HolonomyGroup representing the fundamental group of
        the manifold, together with its holonomy representation.  If
        integer Dehn surgery parameters have been set, then the
        corresponding peripheral elements are killed.

        >>> M = Manifold('m004')
        >>> G = M.fundamental_group()
        >>> G
        Generators:
           a,b
        Relators:
           aaabABBAb
        >>> G.peripheral_curves()
        [('ab', 'aBAbABab')]
        >>> G.SL2C('baaBA')
        matrix([[ (-2.5+2.59807621135j),   (6.06217782649+0.5j)],
                [(-0.866025403784+2.5j),     (4-1.73205080757j)]])

        There are three optional arguments all of which default to True:

        - simplify_presentation
        - fillings_may_affect_generators
        - minimize_number_of_generators

        >>> M.fundamental_group(False, False, False)
        Generators:
           a,b,c
        Relators:
           CbAcB
           BacA
        """
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        name_mangled = 'fundamental_group-%s-%s-%s' %\
                       (simplify_presentation,
                        fillings_may_affect_generators,
                        minimize_number_of_generators)
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = HolonomyGroup(
               self,
               simplify_presentation,
               fillings_may_affect_generators,
               minimize_number_of_generators)
        return self._cache[name_mangled]

    def symmetry_group(self, of_link=False):
        """
        Returns the symmetry group of the Manifold.
        If the flag "of_link" is set, then it only returns symmetries
        that preserves the meridians.
        """
        cdef c_SymmetryGroup* symmetries_of_manifold = NULL
        cdef c_SymmetryGroup* symmetries_of_link = NULL
        cdef c_Triangulation* c_symmetric_triangulation = NULL
        cdef Manifold symmetric_triangulation
        cdef Boolean is_full_group
        cdef c_FuncResult result
        cdef SymmetryGroup symmetry_group

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        name_mangled = 'symmetry_group-%s' % of_link
        if not name_mangled in self._cache.keys():
            result = compute_symmetry_group(self.c_triangulation,
                                            &symmetries_of_manifold,
                                            &symmetries_of_link,
                                            &c_symmetric_triangulation,
                                            &is_full_group)

            if result != func_OK:
                raise ValueError('SnapPea failed to compute any part '
                                 'of the symmetry group.')

            symmetric_triangulation = Manifold('empty')
            symmetric_triangulation.set_c_triangulation(c_symmetric_triangulation)
            self._cache['symmetric_triangulation'] = symmetric_triangulation

            symmetry_group = SymmetryGroup(B2B(is_full_group), True)
            if of_link:
                free_symmetry_group(symmetries_of_manifold)
                symmetry_group._set_c_symmetry_group(symmetries_of_link)
            else:
                free_symmetry_group(symmetries_of_link)
                symmetry_group._set_c_symmetry_group(symmetries_of_manifold)

            self._cache[name_mangled]= symmetry_group
        return self._cache[name_mangled] 
        

    def symmetric_triangulation(self):
        """
        Returns a Dehn filling description of the manifold realizing
        the symmetry group.

        >>> M = Manifold('m003(-3,1)')
        >>> M.symmetry_group()
        D6
        >>> N = M.symmetric_triangulation()
        >>> N
        m003(1,0)(1,0)(1,0)
        >>> N.dehn_fill( [(0,0), (0,0), (0,0)] )
        >>> N.symmetry_group(of_link=True)
        D6
        """
        if not 'symmetric_triangulation' in self._cache:
            self.symmetry_group()
        return self._cache['symmetric_triangulation']
            
    def cover(self, permutation_rep):
        """
        M.cover(permutation_rep)

        Returns a Manifold representing the finite cover
        specified by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the simplified presentation of the
        fundamental group.  Each permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.

        >>> M = Manifold('m004')
        >>> N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])
        >>> abs(N0.volume()/M.volume() - 5) < 0.0000000001
        True
        

        If within Sage, the permutations can also be of type
        PermutationGroupElement, in which case they act on the set
        range(1, d + 1).  Or, you can specify a GAP or Magma subgroup
        of the fundamental group.     Some examples::

          sage: M = snappy.Manifold('m004')

        The basic method::
        
          sage: N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])

        From a Gap subgroup::
        
          sage: G = gap(M.fundamental_group())
          sage: H = G.LowIndexSubgroupsFpGroup(5)[9]
          sage: N1 = M.cover(H)
          sage: N0 == N1
          True

        Or a homomorphism to a permutation group::

          sage: f = G.GQuotients(PSL(2,7))[1]
          sage: N2 = M.cover(f)
          sage: N2.volume()/M.volume()
          7.9999999999999947

        Or maybe we want larger cover coming from the kernel of this::

          sage: N3 = M.cover(f.Kernel())
          sage: N3.volume()/M.volume()
          167.99999999999858

        Check the homology against what Gap computes directly::
        
          sage: N3.homology().betti_number()
          32
          sage: len([ x for x in f.Kernel().AbelianInvariants().sage() if x == 0])
          32

        We can do the same for Magma::

          sage: G = magma(M.fundamental_group())
          sage: Q, f = G.pQuotient(5, 1, nvals = 2)
          sage: M.cover(f.Kernel()).volume()
          10.149416064096533
          sage: h = G.SimpleQuotients(1, 11, 2, 10^4)[1,1]
          sage: N4 = M.cover(h)
          sage: N2 == N4
          True
        """
        cover = Triangulation.cover(self, permutation_rep)
        return Manifold_from_Triangulation(cover, False)

    def covers(self, degree, method=None, cover_type='all'):
        """
        M.covers(degree, method=None)

        Returns a list of Manifolds corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Manifold('m003')
        >>> covers = M.covers(4)
        >>> [(N, N.homology()) for N in covers]
        [(m003~irr~0(0,0)(0,0), Z/5 + Z + Z), (m003~cyc~1(0,0), Z/3 + Z/15 + Z)]

        You can also look just at cyclic covers, which is much faster.

        >>> covers = M.covers(4, cover_type='cyclic')
        >>> [(N, N.homology()) for N in covers]
        [(m003~cyc~0(0,0), Z/3 + Z/15 + Z)]

        If you are using Sage, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional
        argument method = 'gap' If you have Magma installed, you can
        used it to do the heavy lifting by specifying method='magma'.
        """
        covers = Triangulation.covers(self, degree, method,cover_type)
        return [Manifold_from_Triangulation(cover, False) for cover in covers]
    
    def volume(self, accuracy=False, complex_volume=False):
        """
        Returns the volume of the current solution to the hyperbolic
        gluing equations; if the solution is sufficiently non-degenerate,
        this is the sum of the volumes of the hyperbolic pieces in
        the geometric decomposition of the manifold. 

        >>> M = Manifold('m004')
        >>> M.volume()
        2.02988321282
        >>> M.solution_type()
        'all tetrahedra positively oriented'

        The return value has an extra attribute, accuracy, which is the
        number of digits of accuracy as *estimated* by SnapPea.  When
        printing the volume, the result is rounded to 1 more than this
        number of digits.

        >>> M.volume().accuracy
        10
        """
        cdef int acc
        if complex_volume:
            return self.complex_volume()
        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('Solution type is: %s'%solution_type)
             
        vol = SnapPyFloat(volume(self.c_triangulation, &acc))
        vol.accuracy = acc
        if accuracy:
            return (vol, vol.accuracy)
        return vol
        
    def complex_volume(self):
        """
        Returns the complex volume, i.e.
            volume + i 2 pi^2 (chern simons)

        >>> M = Manifold('5_2')
        >>> M.complex_volume()
        (2.8281220883-3.0241283765j)
        """
        if True in self.cusp_info('is_complete'):
            return self.cusped_complex_volume()
        else:
            vol = self.volume()
            cs = self.chern_simons()
            result = SnapPyComplex(vol, 2*math.pi**2 * cs)
            result.accuracy = min(vol.accuracy, cs.accuracy)
            return result

    # cdef hides this method
    cdef cusped_complex_volume(self):
        """
        Returns the complex volume of the manifold, computed using
        Goerner's implementation of Zickert's algorithm.  This only
        works for manifolds with at least one cusp.  A ValueError
        is raised if all cusps are filled.

        >>> M = Manifold('5_2')
        >>> M.cusped_complex_volume()
        (2.8281220883-3.0241283765j)

        The return value has an extra attribute, accuracy, which is
        the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.cusped_complex_volume().accuracy
        9
        """
        cdef Complex vol
        cdef char* err_msg=NULL
        cdef c_Triangulation* copy_c_triangulation
        cdef int accuracy
        cdef c_Triangulation
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        volume = complex_volume(self.c_triangulation, &err_msg, &accuracy)
        # If at first you do not succeed, try again!
        if not err_msg is NULL:
            err_msg = NULL
            copy_triangulation(self.c_triangulation, &copy_c_triangulation)
            randomize_triangulation(copy_c_triangulation)
            volume = complex_volume(copy_c_triangulation, &err_msg, &accuracy)
            free_triangulation(copy_c_triangulation)
        if not err_msg is NULL:
            err_message = err_msg
            raise ValueError(err_message)

        result = SnapPyComplex(volume.real, volume.imag)
        result.accuracy = accuracy
        return result

    def without_hyperbolic_structure(self):
        """
        Returns self as a Triangulation, forgetting the hyperbolic
        structure in the process.

        >>> M = Manifold('9_42')
        >>> T = M.without_hyperbolic_structure()
        >>> type(T)
        <type 'snappy.SnapPy.Triangulation'>
        """
        return Triangulation_from_Manifold(self)

    def _two_to_three(self, tet_num, face_index):
        result = Triangulation._two_to_three(self, tet_num, face_index)
        polish_hyperbolic_structures(self.c_triangulation)
        return result

    def tetrahedra_shapes(self, part=None, fixed_alignment=True, bits_prec=None, dec_prec=None):
        """
        Gives the shapes of the tetrahedra in the current solution to
        the gluing equations.  Returns a list containing one info object
        for each tetrahedron.  The keys are:

        - rect : the shape of the tetrahedron, as a point in the
          complex plane.

        - log : the log of the shape

        - accuracies: a list of the approximate accuracies of the
          shapes, in order (rect re, rect im, log re, log im)

        If the optional variable 'part'is set to one of the above,
        then the function returns only that component of the data.
        
        If the flag 'fixed_alignment' is set to False, then the edges
        used to report the shape parameters are choosen so as to
        normalize the triangle.

        >>> M = Manifold('m015')
        >>> M.tetrahedra_shapes(part='rect')
        [(0.6623589786223731+0.5622795120623011j), (0.662358978622373+0.5622795120623011j), (0.6623589786223729+0.562279512062301j)]
        >>> M.tetrahedra_shapes()
        [{'accuracies': (11, 11, 12, 11), 'log': (-0.14059978716148094+0.7038577213014763j), 'rect': (0.6623589786223731+0.5622795120623011j)},
         {'accuracies': (11, 11, 11, 11), 'log': (-0.14059978716148103+0.7038577213014764j), 'rect': (0.662358978622373+0.5622795120623011j)},
         {'accuracies': (11, 11, 11, 11), 'log': (-0.14059978716148125+0.7038577213014764j), 'rect': (0.6623589786223729+0.562279512062301j)}]

        """        
        cdef double rect_re, rect_im, log_re, log_im
        cdef int acc_rec_re, acc_rec_im, acc_log_re, acc_log_im
        cdef Boolean is_geometric
        
        if self.c_triangulation is NULL: return []

        result = []
        if bits_prec or dec_prec:
            if fixed_alignment == False:
                raise ValueError('High precision shapes must be computed in the fixed alignment')
            shapes = snap.polished_tetrahedra_shapes(self, dec_prec=dec_prec, bits_prec=bits_prec, ignore_solution_type=True)
            for z in shapes:
                result.append(ShapeInfo(rect=z, log=z.log(), accuracies=(None, None, None, None)))
        else:
            for i in range(self.num_tetrahedra()):
                get_tet_shape(self.c_triangulation, i,  fixed_alignment,
                              &rect_re, &rect_im, &log_re, &log_im,
                              &acc_rec_re, &acc_rec_im, &acc_log_re, &acc_log_im,
                              &is_geometric)
                result.append(
                    ShapeInfo(
                        rect=(rect_re + rect_im*(1J)),
                        log=(log_re + log_im*(1J)),
                        accuracies=(acc_rec_re, acc_rec_im,
                                    acc_log_re, acc_log_im)))

        if part != None:
            try:
               return [a[part] for a in result]
            except KeyError:
                raise ValueError('A non-existent shape data type '
                                 'was specified.')
        else:
           return ListOnePerLine(result)

    def set_tetrahedra_shapes(self, shapes, fillings=[(1,0)]):
        """
        M.set_tetrahedra_shapes(shapes, fillings=[(1,0)]):

        Replaces the tetrahedron shapes with those in the list 'shapes'
        and sets the Dehn filling coefficients as specified by the
        fillings argument.
        """
        cdef int i, N
        cdef Complex *shape_array

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        N = get_num_tetrahedra(self.c_triangulation)
        shape_array = <Complex *>malloc(N*sizeof(Complex))
        set_cusps(self.c_triangulation, fillings)
        for i from 0 <= i < N:
            shape = complex(shapes[i]) 
            shape_array[i].real = shape.real
            shape_array[i].imag = shape.imag
        set_tet_shapes(self.c_triangulation, shape_array)
        free(shape_array)

    def solution_type(self, enum=False):
        """
        Returns the type of the current solution to the gluing
        equations, basically a summary of how degenerate the solution
        is.  If the flag enum=True is set, then an integer value is
        returned. The possible answers are:

        - 0: 'not attempted'
        
        - 1: 'all tetrahedra positively oriented' aka 'geometric_solution'
          Should correspond to a genuine hyperbolic structure

        - 2: 'contains negatively oriented tetrahedra' aka 'nongeometric_solution'
          Probably correponds to a hyperbolic structure but some
          simplices have reversed orientiations.  
             
        - 3: 'contains flat tetrahedra' Contains some tetrahedra with
          shapes in R - {0, 1}.

        - 4: 'contains degenerate tetrahedra' Some shapes are close to
          {0,1, or infinity}.  
        
        - 5: 'unrecognized solution type'
        
        - 6: 'no solution found'

        >>> M = Manifold('m007')
        >>> M.solution_type()
        'all tetrahedra positively oriented'
        >>> M.dehn_fill( (3,1) )
        >>> M.solution_type()
        'contains negatively oriented tetrahedra'
        >>> M.dehn_fill( (3,-1) )
        >>> M.solution_type()
        'contains degenerate tetrahedra'
        """
        cdef c_SolutionType solution_type

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        solution_type = get_filled_solution_type(self.c_triangulation)
        if enum:
            return solution_type
        else:
            return SolutionType[solution_type]

    def set_target_holonomy(self, target, which_cusp=0, recompute=True):
        """
        M.set_target_holonomy(target, which_cusp=0, recompute=True)

        Computes a geometric structure in which the Dehn filling curve
        on the specified cusp has holonomy equal to the target value.
        The holonomies of Dehn filling curves on other cusps are left
        unchanged.  If the 'recompute' flag is False, the Dehn filling
        equations are modified, but not solved.
        """
        cdef Complex c_target
        c_target.real = target.real
        c_target.imag = target.imag
        set_target_holonomy(self.c_triangulation, 
                            which_cusp, c_target, recompute)
        
    def cusp_info(self, data_spec=None):
        """
        Returns an info object containing information about the given
        cusp.   Usage:

        >>> M = Manifold('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0)

        To get more detailed information about the cusp, we do

        >>> c = M.cusp_info(0)
        >>> c.shape
        (0.110445017621+0.946770978498j)
        >>> c.modulus
        (-0.12155871955249957+1.042041282932261j)
        >>> sorted(c.keys())
        ['filling', 'holonomies', 'holonomy_accuracy', 'index', 'is_complete', 'modulus', 'shape', 'shape_accuracy', 'topology']

        Here 'shape' is the shape of the cusp, i.e.
        (longitude/meridian)
        and 'modulus' is its shape in the geometrically preferred
        basis, i.e.
        ( (second shortest translation)/(shortest translation)).
        For cusps that are filled, one instead cares about the
        holonomies:
        
        >>> M.cusp_info(-1)['holonomies']
        ((-0.598830888594+1.098125481710j), (0.89824633289+1.49440443102j))

        The complex numbers returned for the shape and for the two
        holonomies have an extra attribute, accuracy, which is
        SnapPea's *estimate* of their accuracy.
        
        You can also get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : complete torus cusp of shape (0.110445017621+0.946770978498j),
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('is_complete')
        [True, False, False]
        """
        return Triangulation.cusp_info(self, data_spec)
            
    def dehn_fill(self, filling_data, which_cusp=None):
        """
        Set the Dehn filling coefficients of the cusps.  This can be
        specified in the following ways, where the cusps are numbered
        by 0,1,...,(num_cusps - 1).  

        - Fill cusp 2:

          >>> M = Manifold('8^4_1')
          >>> M.dehn_fill((2,3), 2)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          8^4_1(0,0)(0,0)(2,3)(1,5)
        
        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          8^4_1(3,0)(1,-4)(2,3)(1,5)

        - When there is only one cusp, there's a shortcut

          >>> N = Manifold('m004')
          >>> N.dehn_fill( (-3,4) )
          >>> N
          m004(-3,4)
        
        Does not return a new Manifold.
        """
        Triangulation.dehn_fill(self, filling_data, which_cusp)
        do_Dehn_filling(self.c_triangulation)
        self._cache = {}

    def set_peripheral_curves(self, peripheral_data,
                              which_cusp=None, return_matrices=False):
        """
        Each cusp has a preferred marking. In the case of a torus
        cusp, this is pair of essential simple curves meeting in one
        point; equivalently, a basis of the first homology of the
        boundary torus. These curves are called the meridian and the
        longitude.

        This method changes these markings in various ways.  In many
        cases, if the flag return_matrices is True then it returns
        a list of change-of-basis matrices is returned, one per
        cusp, which will restore the original markings if passed
        as peripheral_data.

        - Make the shortest curves the meridians, and the second
          shortest curves the longitudes.  
          >>> M = Manifold('5_2')
          >>> M.cusp_info('shape')
          [(-2.49024466751+2.97944706648j)]
          >>> cob = M.set_peripheral_curves('shortest', return_matrices=True)
          >>> M.cusp_info('shape')
          [(-0.49024466751+2.97944706648j)]
          >>> cob
          [[[1, 0], [-2, 1]]]
          >>> M.set_peripheral_curves(cob)
          >>> M.cusp_info('shape')
          [(-2.49024466751+2.97944706648j)]

          You can also make just the meridians as short as 
          possible while fixing the longitudes via the option
          'shortest_meridians', and conversely with
          'shortest_longitudes'.  
          
        - If cusps are Dehn filled, make those curves meridians.  

          >>> M = Manifold('m125(0,0)(2,5)')
          >>> M.set_peripheral_curves('fillings')
          >>> M
          m125(0,0)(1,0)
        
        - Change the basis of a particular cusp, say the first one:

          >>> M.set_peripheral_curves( [ (1,2), (1,3) ] , 0)

          Here (1,2) is the new meridian written in the old basis, and
          (1,3) the new longitude.
          
        - Change the basis of all the cusps at once

          >>> new_curves = [ [(1,-1), (1,0)],  [(3,2), (-2,-1)] ]
          >>> M.set_peripheral_curves(new_curves)
          >>> M
          m125(0,0)(-1,-2)
        """
        cdef int a,b,c,d
        cdef MatrixInt22 *matrices
        cdef c_FuncResult result 

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty')

        if which_cusp != None:
           try:
              which_cusp = range(self.num_cusps())[which_cusp]
           except IndexError:
              raise IndexError('The specified cusp (%s) does not '
                               'exist.'%which_cusp)

        self._cache = {}

        if peripheral_data == 'shortest_meridians':
            # For each cusp, replace its current meridian with the
            # shortest one possible.
            cusps = [which_cusp] if which_cusp else range(self.num_cusps())
            cobs = []
            for i in cusps:
                d = round((1/self.cusp_info(i)['shape']).real, 0)
                self.set_peripheral_curves([(1, -d), (0,1)], i)
                cobs.append([[1,d],[0,1]])
            if return_matrices:
                return cobs
            
        elif peripheral_data == 'shortest_longitudes':
            # For each cusp, replace its current longitude with the
            # shortest one possible.
            cusps = [which_cusp] if which_cusp else range(self.num_cusps())
            cobs = []
            for i in cusps:
                d = round(self.cusp_info(i)['shape'].real, 0)
                self.set_peripheral_curves([(1, 0), (-d,1)], i)
                cobs.append([[1,0],[d,1]])
            if return_matrices:
                return cobs
                    
        elif peripheral_data == 'fillings':
            if which_cusp != None:
                raise ValueError("You must apply 'fillings' to all "
                                 "of the cusps.")
            install_current_curve_bases(self.c_triangulation)
            return

        elif peripheral_data == 'shortest':
            if which_cusp != None:
                raise ValueError("You must apply 'shortest' to all "
                                 "of the cusps.")
            if return_matrices:
                matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                                 sizeof(MatrixInt22))
                install_shortest_with_matrices(self.c_triangulation, matrices)
                cobs = []
                for n in range(self.num_cusps()):
                    cobs.append([ [ matrices[n][1][1], -matrices[n][0][1] ],
                                  [ -matrices[n][1][0], matrices[n][0][0] ] ])
                free(matrices)
                return cobs
            else:
                install_shortest_bases(self.c_triangulation)

        elif peripheral_data == 'combinatorial':
            if return_matrices:
                matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                                 sizeof(MatrixInt22))
                install_combinatorial_bases(self.c_triangulation, matrices)
                cobs = []
                for n in range(self.num_cusps()):
                    cobs.append([ [ matrices[n][0][0], matrices[n][0][1] ],
                                  [ matrices[n][1][0], matrices[n][1][1] ] ])
                free(matrices)
                return cobs
            else:
                peripheral_curves(self.c_triangulation)
                
        elif which_cusp != None:
            meridian, longitude = peripheral_data
            a, b = meridian
            c, d = longitude
            if a*d - b*c != 1:
                raise ValueError('The data provided does not give a '
                                 '(positively oriented) basis.')

            matrices = <MatrixInt22 *>malloc(self.num_cusps() *
                                             sizeof(MatrixInt22))

            for n in range(self.num_cusps()):
                for i,j in [(0,0),(0,1),(1,0),(1,1)]:
                    matrices[n][i][j] = 1 if i == j else 0

            matrices[which_cusp][0][0] = a
            matrices[which_cusp][0][1] = b
            matrices[which_cusp][1][0] = c
            matrices[which_cusp][1][1] = d
            result = change_peripheral_curves(self.c_triangulation, matrices)
            if result == func_bad_input:
                raise ValueError('The peripheral data ((%d, %d), (%d,%d)) '
                                 'is not acceptable.' % (a,b,c,d))
            free(matrices)
            
        else:
            if self.num_cusps() == 1 and len(peripheral_data) == 2:
                self.set_peripheral_curves(peripheral_data, 0)
                return 
            if len(peripheral_data) > self.num_cusps():
                raise IndexError('You provided more peripheral data '
                                 'than there are cusps.')
            for i, basis in enumerate(peripheral_data):
                self.set_peripheral_curves(basis, i)

    def orientation_cover(self):
        """
        For a non-orientable manifold, returns the 2-fold cover which
        is orientable.

        >>> X = Manifold('x123')
        >>> Y = X.orientation_cover()
        >>> (X.is_orientable(), Y.is_orientable())
        (False, True)
        >>> Y
        x123~(0,0)(0,0)
        """ 
        return self.without_hyperbolic_structure().orientation_cover().with_hyperbolic_structure()


    def dual_curves(self, max_segments=6):
        """
        Constructs a *reasonable* selection of simple closed curves in
        a manifold's dual 1-skeleton.  In particular, it returns those
        that appear to represent geodesics. The resulting curves can
        be drilled out.

        >>> M = Manifold('m015')
        >>> curves = M.dual_curves()
        >>> curves
        [  0: orientation-preserving curve of length (0.5623991486459233-2.815430885205906j),
           1: orientation-preserving curve of length (1.124798297291847+0.6523235367677742j),
           2: orientation-preserving curve of length (1.260804017474151+1.978046890227184j),
           3: orientation-preserving curve of length (1.5882693259837328+1.6734716736926436j),
           4: orientation-preserving curve of length (1.6871974459377679+2.8154308852059073j)]

        Each curve is returned as an info object with these keys
        
        >>> curves[0].keys()
        ['index', 'filled_length', 'complete_length', 'max_segments', 'parity']
        
        We can drill out any of these curves to get a new manifold
        with one more cusp.

        >>> N = M.drill(curves[0])
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        
        By default, this function only finds curves of length 6; this
        can be changed by specifying the optional argument
        max_segments

        >>> M.dual_curves(max_segments=2)
        [  0: orientation-preserving curve of length (0.5623991486459233-2.815430885205906j)]
        """
        cdef int i, num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_MatrixParity parity
        cdef Complex complete_length, filled_length

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)
        result = []
        for i from 0 <= i < num_curves:
            get_dual_curve_info(
                curve_list[i], 
                &complete_length,
                &filled_length,
                &parity
                )
            result.append(
                DualCurveInfo(
                    index=SnapPyInt(i),
                    parity=SnapPyInt(parity),
                    filled_length=SnapPyComplex(C2C(filled_length)),
                    complete_length=SnapPyComplex(C2C(complete_length)),
                    max_segments=SnapPyInt(max_segments)
                  )
               )
        free_dual_curves(num_curves, curve_list)
        return ListOnePerLine(result)
    
    def length_spectrum(self, cutoff=1.0):
        """
        M.length_spectrum(cutoff=1.0)

        Returns a list of geodesics (with multiplicities) of length
        up to the specified cutoff value. (The default cutoff is 1.0.)
        """
        try:
            D = DirichletDomain(self)
        except:
            raise RuntimeError('The length spectrum not available: '
                                'no Dirichlet Domain.')
        return D.length_spectrum_dicts(cutoff_length=cutoff)

    # cdef will hide this method.
    cdef old_chern_simons(self):
        """
        Compute the Chern-Simons invariant using SnapPea's original
        algorithm, which is based on Meyerhoff-Hodgson-Neumann.
        """
        cdef Boolean is_known, requires_initialization
        cdef double CS
        cdef int accuracy

        if self.c_triangulation is NULL: return 0
        get_CS_value(self.c_triangulation,
                     &is_known,
                     &CS,
                     &accuracy,
                     &requires_initialization)
        if not is_known:
           raise ValueError("The Chern-Simons invariant isn't "
                            "currently known.") 
        cs = SnapPyFloat(CS)
        cs.accuracy = accuracy
        return cs

    def chern_simons(self):
        """
        Returns the Chern-Simons invariant of the manifold, if it is known.

        >>> M = Manifold('m015')
        >>> M.chern_simons()
        -0.1532041333

        The return value has an extra attribute, accuracy, which
        is the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.chern_simons().accuracy
        9

        By default, when the manifold has at least one cusp, Zickert's
        algorithm is used; when the manifold is closed we use SnapPea's
        original algorithm, which is based on Meyerhoff-Hodgson-Neumann.
        
        Note: When computing the Chern-Simons invariant of a closed
        manifold, one must sometimes compute it first for the unfilled
        manifold so as to initialize SnapPea's internals.  For instance,

        >>> M = Manifold('5_2')
        >>> M.chern_simons()
        -0.1532041333
        >>> M.dehn_fill( (1,2) )
        >>> M.chern_simons()
        0.0773178713861

        works, but will fail with 'Chern-Simons invariant not
        currently known' if the first call to chern_simons is not
        made.
        """
        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('The solution type is: %s'%solution_type)

        if not True in self.cusp_info('is_complete'):
           cs = self.old_chern_simons()
        else:
            volume = self.cusped_complex_volume()
            cs = SnapPyFloat( volume.imag/(2.0*math.pi**2) )
            set_CS_value(self.c_triangulation, cs)
            cs.accuracy = volume.accuracy
        return cs
        
    def drill(self, which_curve, max_segments=6):
        """
        Drills out the specified dual curve from among all dual curves
        with at most max_segments, which defaults to 6. The method
        dual_curve allows one to see the properties of curves before
        chosing which one to drill out.

        >>> M = Manifold('v3000')
        >>> N = M.drill(0, max_segments=3)
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        """
        
        cdef int num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_Triangulation *c_triangulation
        cdef Triangulation result
        cdef char* c_new_name

        if isinstance(which_curve, DualCurveInfo):
            max_segments = which_curve.max_segments
            which_curve = which_curve.index

        new_name = to_byte_str(self.name() + '-%d'%which_curve)
        c_new_name = new_name

        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)

        if which_curve not in range(num_curves):
            raise IndexError('The drilling curve requested is not '
                             'in range(%d).' % num_curves)
        
        c_triangulation = drill_cusp(self.c_triangulation,
                                     curve_list[which_curve],
                                     c_new_name)
        free_dual_curves(num_curves, curve_list)

        if c_triangulation == NULL:
            raise RuntimeError('The curve is not isotopic to a geodesic.')
        else:
            result = Manifold('empty')
            result.set_c_triangulation(c_triangulation)
            return result

    def is_isometric_to(self, Manifold other, return_isometries=False):
        """
        Returns True if M and N are isometric, False otherwise.

        >>> M = Manifold('m004')
        >>> N = Manifold('4_1')
        >>> K = Manifold('5_2')
        >>> M.is_isometric_to(N)
        True
        >>> N.is_isometric_to(K)
        False

        We can also get a complete list of isometries between the two
        manifolds:

        >>> M = Manifold('5^2_1')  # The Whitehead link
        >>> N = Manifold('m129')
        >>> isoms = M.is_isometric_to(N, return_isometries = True)
        >>> isoms[6]  # Includes action on cusps
        0 -> 1  1 -> 0 
        [1  2]  [-1 -2]
        [0 -1]  [ 0  1]
        Extends to link

        Each transformation between cusps is given by a matrix which
        acts on the left.  That is, the two *columns* of the matrix
        give the image of the meridian and longitude respectively.  In
        the above example, the meridian of cusp 0 is sent to the
        meridian of cusp 1.
        
        Note: The answer True is rigorous, but the answer False may
        not be as there could be numerical errors resulting in finding
        an incorrect canonical triangulation.
        """
        cdef Boolean are_isometric
        cdef c_FuncResult result
        cdef IsometryList *isometries = NULL

        if self.c_triangulation is NULL or other.c_triangulation is NULL:
            raise ValueError('Manifolds must be non-empty.')

        if return_isometries:
            result = compute_isometries(self.c_triangulation,
                                        other.c_triangulation, 
                                        &are_isometric,
                                        &isometries, NULL)
        else:
            result = compute_isometries(self.c_triangulation,
                                        other.c_triangulation, 
                                        &are_isometric, NULL, NULL)
            
        if FuncResult[result] == 'func_bad_input':
            raise ValueError('The Dehn filling coefficients must be '
                             'relatively prime integers.')

        if FuncResult[result] == 'func_failed':
            raise RuntimeError('SnapPea failed to determine whether '
                               'the manifolds are isometric.')

        ans = bool(are_isometric)
        if return_isometries:
            if not ans:
                return []
            else:
                ans = IsometryListToIsometries(isometries)
                free_isometry_list(isometries)

        return ans 

    def is_two_bridge(self):
        """
        If the manifold is the complement of a two-bridge knot or link
        in S^3, then this method returns (p,q) where p/q is the
        fraction describing the link.  Otherwise, returns False.

        >>> M = Manifold('m004')
        >>> M.is_two_bridge()
        (2, 5)
        >>> M = Manifold('m016')
        >>> M.is_two_bridge()
        False
        
        Note: An answer of 'True' is rigorous, but not the answer
        'False', as there could be numerical errors resulting in
        finding an incorrect canonical triangulation.
        """
        cdef Boolean is_two_bridge
        cdef long int p, q
        cdef c_Triangulation *c_canonized_triangulation
        
        if self.c_triangulation is NULL: return False

        copy_triangulation(self.c_triangulation, &c_canonized_triangulation)
        proto_canonize(c_canonized_triangulation)
        two_bridge(c_canonized_triangulation, &is_two_bridge, &p, &q)        
        return (p,q) if  is_two_bridge else False

    def _choose_generators(self, compute_corners, centroid_at_origin):
        choose_generators(self.c_triangulation,
                          compute_corners,
                          centroid_at_origin)

    def _choose_generators_info(self):
        """
        Extracts, from the bowels of SnapPea, the information about the
        underlying generators of the fundamental group.  Returns a
        list with one entry for each tetrahedra.
        """
        cdef int generator_path, face0_gen, face1_gen, face2_gen, face3_gen
        cdef Complex c0, c1, c2, c3
        ans = []
        for i in range(self.num_tetrahedra()):
            choose_gen_tetrahedron_info(self.c_triangulation,
                                        i, &generator_path,
                                        &face0_gen, &face1_gen,
                                        &face2_gen, &face3_gen,
                                        &c0, &c1, &c2, &c3)
            ans.append(
                {'index':i,
                 'generators':(face0_gen, face1_gen, face2_gen, face3_gen),
                 'corners': (C2C(c0), C2C(c1), C2C(c2), C2C(c3)),
                 'generator_path':generator_path}
                )
        return ans

# Conversion functions Manifold <-> Triangulation

def Manifold_from_Triangulation(Triangulation T, recompute=True):
    cdef c_Triangulation *c_triangulation
    cdef Manifold M

    copy_triangulation(T.c_triangulation, &c_triangulation)
    M = Manifold('empty')
    M.set_c_triangulation(c_triangulation)
    if recompute:
        find_complete_hyperbolic_structure(c_triangulation)
        do_Dehn_filling(c_triangulation)
    M.set_name(T.name())
    return M

def Triangulation_from_Manifold(Manifold M):
    cdef c_Triangulation *c_triangulation
    cdef Triangulation T

    copy_triangulation(M.c_triangulation, &c_triangulation)
    remove_hyperbolic_structures(c_triangulation)
    T = Triangulation('empty')
    T.set_c_triangulation(c_triangulation)
    T.set_name(M.name())
    return T


# Fundamental Groups

Alphabet = '$abcdefghijklmnopqrstuvwxyzZYXWVUTSRQPONMLKJIHGFEDCBA'

# Helper functions for manipulating fg. words

def inverse_word(word):
    parts = list(word.swapcase())
    parts.reverse()
    return ''.join(parts)

reduce_word_regexp = re.compile('|'.join([x + x.swapcase()
                                          for x in string.ascii_letters]))

def reduce_word(word):
    """
    Cancels inverse generators.
    """
    ans, progress = word, 1

    while progress:
        ans, progress =  reduce_word_regexp.subn('', ans)
    return ans

def format_word(word, verbose_form):
    return word if not verbose_form else '*'.join(
        [a if a.islower() else a.lower() + '^-1' for a in list(word)]
        )

cdef class CFundamentalGroup:
    cdef c_GroupPresentation *c_group_presentation
    cdef c_Triangulation *c_triangulation
    cdef readonly num_cusps

    cdef c_word_as_int_list(self, int *word):
        cdef int n = 0
        word_list = []
        while word[n] != 0:
            word_list.append(word[n])
            n += 1
        return word_list

    cdef c_word_as_string(self, int *word):
        cdef int n = 0
        word_list = []
        while word[n] != 0:
            word_list.append(Alphabet[word[n]])
            n += 1
        return ''.join(word_list)

    cdef int *c_word_from_list(self, word_list):
        cdef int *c_word, length, size, n
        length = len(word_list)
        size = sizeof(int)*(1+length)
        c_word = <int *>malloc(size)
        for n from 0 <= n < length:
            c_word[n] = word_list[n]
        c_word[length] = 0
        return c_word

    def __cinit__(self, Triangulation triangulation,
                  simplify_presentation = True,
                  fillings_may_affect_generators = True,
                  minimize_number_of_generators = True):
        if triangulation.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        copy_triangulation(triangulation.c_triangulation,
                           &self.c_triangulation)
        self.c_group_presentation = fundamental_group(
            self.c_triangulation,
            simplify_presentation,
            fillings_may_affect_generators,
            minimize_number_of_generators)
        self.num_cusps = triangulation.num_cusps()

    def __dealloc__(self):
        free_triangulation(self.c_triangulation)
        free_group_presentation(self.c_group_presentation)

    def __repr__(self):
        return 'Generators:\n   %s\nRelators:\n   %s'%(
            ','.join(self.generators()),
            '\n   '.join(self.relators()))

    def _word_as_list(self, word):
        if not isinstance(word, basestring):
            raise TypeError('Words must be represented '
                            'as Python strings.')
        word_list = []
        generators = self.generators()
        for letter in word:
            try:
                if letter.islower():
                    word_list.append(1 + generators.index(letter))
                else:
                    word_list.append(-1 - generators.index(letter.lower()))
            except ValueError:
                raise RuntimeError('The word contains a non-generator.')
        return word_list

    def num_generators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_generators(self.c_group_presentation)

    def num_relators(self):
        """
        Return the number of generators for the presentation.
        """
        return fg_get_num_relations(self.c_group_presentation)
                            
    def num_original_generators(self):
        """
        Return the number of geometric generators (before simplification).
        """
        return fg_get_num_orig_gens(self.c_group_presentation)

    def original_generators(self, verbose_form=False):
        """
        Return the original geometric generators (before
        simplification) in terms of the current generators.
        """
        cdef int n
        cdef int *gen
        orig_gens = []
        num_orig_gens = fg_get_num_orig_gens(self.c_group_presentation)
        for n from 0 <= n < num_orig_gens:
            gen = fg_get_original_generator(self.c_group_presentation, n)
            word = format_word(self.c_word_as_string(gen), verbose_form)
            orig_gens.append(word)
            fg_free_relation(gen)
        return orig_gens

    def generators_in_originals(self, verbose_form=False, raw_form =False):
        """
        Return the current generators in terms of the original
        geometric generators (before simplification).

        If the flag "raw_form" is set to True, it returns a sequence of
        instructions for expressing the current generators in terms of
        the orignal ones.  This is sometimes much more concise, though
        the format is somewhat obscure.  See the source code of this
        function in SnapPy.pyx for details. 
        """
        moves = self._word_moves()
        if raw_form:
            return moves

        n = self.num_original_generators()
        if n > 26:
            raise ValueError('Too many generators.')

        words = [None] + list(string.ascii_letters[:n])

        while len(moves) > 0:
            a = moves.pop(0)
            if a >= len(words): # new generator added
                n = moves.index(a)  # end symbol location
                # word is the expression of the new generator in terms
                # of the old ones
                word, moves = moves[:n], moves[n+1:]
                words.append( reduce_word(''.join(
                    [words[g] if g > 0 else inverse_word(words[-g])
                     for g in word]
                    )))
            else:
                b = moves.pop(0)
                if a == b:  # generator removed
                    words[a] = words[-1]
                    words = words[:-1]
                elif a == -b: # invert generator
                    words[a] = inverse_word(words[a])
                else: #handle slide
                    A, B = words[abs(a)], words[abs(b)]
                    if a*b < 0:
                        B = inverse_word(B)
                    words[abs(a)] = reduce_word(  A+B if a > 0 else B+A ) 

        return [format_word( w, verbose_form) for w in words[1:]]

    def _word_moves(self):
        cdef int *c_moves
        c_moves = fg_get_word_moves(self.c_group_presentation)
        moves = self.c_word_as_int_list(c_moves)
        fg_free_relation(c_moves)
        return moves
        
    def generators(self):
        """
        Return the letters representing the generators in the presentation.
        """
        return [ Alphabet[i] for i in range(1, 1 + self.num_generators()) ]

    def relators(self, verbose_form = False, as_int_list = False):
        """
        Return a list of words representing the relators in the presentation.

        If the optional argument verbose_form is True, then the
        relator is returned in the form "a*b*a^-1*b^-1" instead of "abAB".  
        """
        cdef int n
        cdef int *relation
        relation_list = []
        num_relations = fg_get_num_relations(self.c_group_presentation)
        for n from 0 <= n < num_relations:
            relation = fg_get_relation(self.c_group_presentation, n)
            if as_int_list:
                word = self.c_word_as_int_list(relation)
            else:
                word = format_word(self.c_word_as_string(relation), verbose_form)
            relation_list.append(word)
            fg_free_relation(relation)
        return relation_list

    def meridian(self, int which_cusp=0, as_int_list = False):
        """
        Returns a word representing a conjugate of the current
        meridian for the given cusp.  Guaranteed to commute with the
        longitude for the same cusp.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.meridian(0)
        'aaba'
        >>> G.meridian(-1)  # The last cusp
        'baaba'
        """
        try:
            which_cusp = range(self.num_cusps)[which_cusp]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)
        if as_int_list:
            return self.c_word_as_int_list(
               fg_get_meridian(self.c_group_presentation, which_cusp))
        else:
            return self.c_word_as_string(
               fg_get_meridian(self.c_group_presentation, which_cusp))

    def longitude(self, int which_cusp=0, as_int_list = False):
        """
        Returns a word representing a conjugate of the current
        longitude for the given cusp.  Guaranteed to commute with the
        meridian for the same cusp.  Note: for Klein bottle cusps,
        the longitude must be defined carefully.

        >>> G = Manifold('m004').fundamental_group()
        >>> G.longitude(0)
        'aBAbABab'
        >>> G.longitude()   # shortcut for the above.  
        'aBAbABab'
        """
        try:
            which_cusp = range(self.num_cusps)[which_cusp]
        except IndexError:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)

        if as_int_list:
            return self.c_word_as_int_list(
               fg_get_longitude(self.c_group_presentation, which_cusp))
        else:
            return self.c_word_as_string(
               fg_get_longitude(self.c_group_presentation, which_cusp))

    def peripheral_curves(self, as_int_list = False):
        """
        Returns a list of meridian-longitude pairs for all cusps.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.peripheral_curves()
        [('aaba', 'abb'), ('baaba', 'Ba')]
        """
        return [ (self.meridian(n, as_int_list),
                  self.longitude(n, as_int_list))
                 for n in range(self.num_cusps) ]

    def magma_string(self):
        """
        Returns a string which will define this group within MAGMA.
        """
        return ('Group<' + ','.join(self.generators()) + '|' +
                ', '.join(self.relators(verbose_form = True)) + '>')

    def gap_string(self):
        """
        Returns a string which will define this group within GAP.
        """
        gens = ', '.join(self.generators())
        gen_names = ', '.join(['"' + x + '"' for x in self.generators()])
        relators = ', '.join(self.relators(verbose_form = True))
        assignments = ''.join(
            ['%s := F.%d; ' % (x, i+1)
             for (i, x) in enumerate(self.generators())]
            )
        return ('CallFuncList(function() local F, %s; '
                'F := FreeGroup(%s); %s  return F/[%s]; end,[])'
                % (gens, gen_names, assignments, relators)
                )

    def _gap_init_(self):
        return self.gap_string()

    def _magma_init_(self, magma):
        return self.magma_string()

class FundamentalGroup(CFundamentalGroup):
    """
    A FundamentalGroup represents a presentation of the fundamental
    group of a SnapPea Triangulation.  Group elements are described as
    words in the generators a,b,..., where the inverse of a is denoted
    A.  Words are represented by Python strings (and the concatenation
    operator is named "+", according to Python conventions).

    Instantiate as T.fundamental_group(), where T is a Triangulation.

    Methods:
        num_generators() --> number of generators
        num_relators()   --> number of relators
        generators()     --> list of generators
        relators()       --> list of relators
        meridian(n)      --> word representing the meridian on cusp #n
        longitude(n)     --> word representing the longitude on cusp #n
    """

if _within_sage:
    FundamentalGroup.__bases__ += (sage.structure.sage_object.SageObject,)

# Holonomy Groups

cdef C2C(Complex C):
    return complex(C.real, C.imag)

cdef B2B(Boolean B):
    return B != 0

cdef class CHolonomyGroup(CFundamentalGroup):
    def _matrices(self, word):
        """
        Returns (M,O,L) where M = SL2C(word), O = O31(word), and L is
        the complex length.
        """
        cdef MoebiusTransformation M 
        cdef O31Matrix O
        cdef int *c_word
        cdef c_FuncResult result
        word_list = self._word_as_list(word)
        c_word = self.c_word_from_list(word_list)
        result = fg_word_to_matrix(self.c_group_presentation, c_word, O, &M)
        if result == 0:
            sl2 = matrix([[C2C(M.matrix[0][0]), C2C(M.matrix[0][1])],
                           [C2C(M.matrix[1][0]), C2C(M.matrix[1][1])]]) 
            o31 = matrix([
                [O[0][0], O[0][1], O[0][2], O[0][3]],
                [O[1][0], O[1][1], O[1][2], O[1][3]],
                [O[2][0], O[2][1], O[2][2], O[2][3]],
                [O[3][0], O[3][1], O[3][2], O[3][3]]
                ])
            L = C2C(complex_length_mt(&M))
            return sl2, o31, L
        else:
            return None

    def SL2C(self, word):
        """
        Return the image of the element represented by the input word
        under some SL(2,C) representation that lifts the holonomy
        representation.  Note: the choice of lift is not guaranteed to
        vary continuously when filling coefficients are changed.
        """
        return self._matrices(word)[0]

    def O31(self, word):
        """
        Return the image of the element represented by the input word
        under the holonomy representation, where Isom(H^3) is
        identified with SO(3,1).
        """
        return self._matrices(word)[1]

    def complex_length(self, word):
        """
        Return the complex length of the isometry represented by the
        input word.
        """
        return self._matrices(word)[2]

    def _choose_generators_info(self):
        """
        Extracts from the bowels of SnapPea the information about the
        underlying generators of the fundamental group.  Returns a
        list with one entry for each tetrahedra.
        """
        
        cdef int face0_gen, face1_gen, face2_gen, face3_gen
        cdef int generator_path, N, i
        cdef Complex c0, c1, c2, c3
        ans = []
        N = get_num_tetrahedra(self.c_triangulation)
        for i from 0 <= i < N:
            choose_gen_tetrahedron_info(
                self.c_triangulation, i, &generator_path,
                &face0_gen, &face1_gen, &face2_gen, &face3_gen,
                &c0, &c1, &c2, &c3)
            ans.append(
                {'index':i,
                 'generators':(face0_gen, face1_gen, face2_gen, face3_gen),
                 'corners': (C2C(c0), C2C(c1), C2C(c2), C2C(c3)),
                 'generator_path':generator_path
                 }
                )
        return ans


class HolonomyGroup(CHolonomyGroup):
    """
    A HolonomyGroup is a FundamentalGroup with added structure
    consisting of a holonomy representation into O(3,1), and an
    arbitrarily chosen lift of the holonomy representation to SL(2,C).
    The holonomy is determined by the shapes of the tetrahedra, so a
    HolonomyGroup is associated to a Manifold, while a Triangulation
    only has a FundamentalGroup.  Methods are provided to evaluate the
    representations on a group element.

    A FundamentalGroup represents a presentation of the fundamental
    group of a SnapPea Triangulation.  Group elements are described as
    words in the generators a,b,..., where the inverse of a is denoted
    A.  Words are represented by python strings (and the concatenation
    operator is named '+', according to Python conventions).

    Instantiate via M.fundamental_group(), where M is a Manifold.
    """

if _within_sage:
    HolonomyGroup.__bases__ += (sage.structure.sage_object.SageObject,)


# Dirichlet Domains

cdef WEPolyhedron *read_generators_from_file(file_name,
                                             vertex_epsilon=10.0**-8):
    DET_ERROR_EPSILON = 10.0**-3
    
    data = open(file_name).readlines()
    if data[0].strip() != '% Generators':
        raise ValueError('The generator file does not start with '
                         '"% Generators"')
    nums = []
    for line in data[1:]:
        nums +=  line.split()
    num_gens = int(nums[0])
    nums = [float(f) for f in nums[1:]]

    cdef O31Matrix *generators
    cdef MoebiusTransformation *temp_gens
    if len(nums) == 16 * num_gens:
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            for j in range(4):
                for k in range(4):
                    generators[i][j][k] =  nums.pop(0)
    elif len(nums) == 8*num_gens:
        temp_gens = <MoebiusTransformation *>malloc(
            num_gens*sizeof(MoebiusTransformation))
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            temp_gens[i].parity = orientation_preserving
            for j in range(2):
                for k in range(2):
                    temp_gens[i].matrix[j][k].real = nums.pop(0)
                    temp_gens[i].matrix[j][k].imag = nums.pop(0)
            #a = C2C(temp_gens[i].matrix[0][0])
            #b = C2C(temp_gens[i].matrix[0][1])
            #c = C2C(temp_gens[i].matrix[1][0])
            #d = C2C(temp_gens[i].matrix[1][1])
            #print a, b
            #print c, d 
            #print a*d - b*c
        Moebius_array_to_O31_array(temp_gens, generators, num_gens)
        free(temp_gens)
    else:
        raise ValueError('The amount of data given is not consistent '
                         'with %d O31 or SL2C matrices.' % num_gens)

    if not O31_determinants_OK(generators, num_gens, DET_ERROR_EPSILON):
        raise ValueError('The data given do not have the '
                         'right determinants.')
        
    cdef WEPolyhedron *dirichlet_domain
    dirichlet_domain = Dirichlet_from_generators(generators,
                                                 num_gens,
                                                 vertex_epsilon,
                                                 Dirichlet_keep_going,
                                                 True);
    return dirichlet_domain

cdef class CDirichletDomain:
    cdef WEPolyhedron *c_dirichlet_domain
    cdef c_Triangulation *c_triangulation

    def __cinit__(self, Manifold manifold=None,
                      vertex_epsilon=10.0**-8,
                      displacement = [0.0, 0.0, 0.0],
                      centroid_at_origin=True,
                      maximize_injectivity_radius=True, generator_file = None):
        cdef double c_displacement[3]

        if generator_file != None:
            self.c_dirichlet_domain = read_generators_from_file(
                generator_file)
            self.manifold_name = generator_file
        else:
            if manifold.c_triangulation is NULL:
                raise ValueError('The Triangulation is empty.')
            for n from 0 <= n < 3:
                c_displacement[n] = displacement[n] 
            copy_triangulation(manifold.c_triangulation,
                               &self.c_triangulation)
            self.c_dirichlet_domain = Dirichlet_with_displacement(
                self.c_triangulation,
                c_displacement,
                vertex_epsilon,
                centroid_at_origin,
                Dirichlet_keep_going,
                maximize_injectivity_radius )
            if self.c_dirichlet_domain == NULL:
                raise RuntimeError('The Dirichlet construction failed.')
            self.manifold_name = manifold.name()

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_dirichlet_domain != NULL:
            free_Dirichlet_domain(self.c_dirichlet_domain)

    def __repr__(self):
        return '%d finite vertices, %d ideal vertices; %d edges; %d faces'%(
            self.num_finite_vertices(),
            self.num_ideal_vertices(),
            self.num_edges(),
            self.num_faces(),
            )

    def num_vertices(self):
        """
        Return the number of vertices.
        """
        return self.c_dirichlet_domain.num_vertices

    def num_finite_vertices(self):
        """
        Return the number of finite (non-ideal) vertices.
        """
        return self.c_dirichlet_domain.num_finite_vertices

    def num_ideal_vertices(self):
        """
        Return the number of ideal vertices.
        """
        return self.c_dirichlet_domain.num_ideal_vertices

    def num_edges(self):
        """
        Return the number of edges.
        """
        return self.c_dirichlet_domain.num_edges

    def num_faces(self):
        """
        Return the number of faces.
        """
        return self.c_dirichlet_domain.num_faces

    def in_radius(self):
        """
        Return the radius of the largest inscribed sphere.
        """
        return self.c_dirichlet_domain.inradius

    def out_radius(self):
        """
        Return the radius of the smallest circubscribed sphere.
        """
        return self.c_dirichlet_domain.outradius

    def length_spectrum_dicts(self, cutoff_length=1.0,
                        full_rigor=True,
                        multiplicities=True,
                        user_radius=0.0):
        """
        Return a list of info objects describing the short
        geodesics up to the specified cutoff length.  The keys are
        'length', 'parity', 'topology', and 'multiplicity'.  The
        length is the complex length; the parity specifies whether
        orientation is preserved; and topology distinguishes between
        circles and mirrored intervals.
        """
        cdef int num_lengths
        cdef MultiLength* geodesics
        length_spectrum(self.c_dirichlet_domain,
                        cutoff_length,
                        full_rigor,
                        multiplicities,
                        user_radius,
                        &geodesics,
                        &num_lengths)
        spectrum = []
        for n from 0 <= n < num_lengths:
            spectrum.append(
               LengthSpectrumInfo(
                  length=SnapPyComplex(C2C(geodesics[n].length)),
                  parity=SnapPyStr(MatrixParity[geodesics[n].parity]),
                  topology=SnapPyStr(Orbifold1[geodesics[n].topology]),
                  multiplicity=SnapPyInt(geodesics[n].multiplicity)
                  )
               )
        free_length_spectrum(geodesics)
        return LengthSpectrum(spectrum)

    def vertex_list(self):
        """
        Return a list of the coordinates of the vertices.  These are
        the three space coordinates of a point in the time=1 slice of
        Minkowski space.  That is to say, these are the coordinates of
        the image of the point under projection into the Klein model.
        """
        cdef WEVertex *vertex = &self.c_dirichlet_domain.vertex_list_begin
        vertices = []
        vertex = vertex.next
        while vertex != &self.c_dirichlet_domain.vertex_list_end:
          vertices.append( (vertex.x[1], vertex.x[2], vertex.x[3]) )
          vertex = vertex.next
        return vertices

    def face_list(self):
        """
        Return a list of faces, each represented as a dictionary with
        keys 'vertices', 'distance', 'closest', 'hue'.  The distance
        from the origin is the value for 'distance', and the value for
        'closest' is the orthogonal projection of the origin to the
        plane containing the face.  The vertices of each face are
        listed in clockwise order, as viewed from outside the
        polyhedron.
        """
        cdef WEFace *face = &self.c_dirichlet_domain.face_list_begin
        cdef WEEdge *edge
        cdef WEVertex *vertex
        # WE enums -- see winged_edge.h
        left, right, tail, tip = 0, 1, 0, 1
        faces = []
        face = face.next
        while face != &self.c_dirichlet_domain.face_list_end:
            vertices = []
            edge = face.some_edge
            while True:
                # find the vertex at the counter-clockwise end
                if edge.f[left] == face:
                    vertex = edge.v[tip]
                else:
                    vertex = edge.v[tail]
                vertices.append( (vertex.x[1], vertex.x[2], vertex.x[3]) )
                # get the next edge
                if edge.f[left] == face:
                    edge = edge.e[tip][left];
                else:
                    edge = edge.e[tail][right];
                if edge == face.some_edge:
                    break
            faces.append(
                {'vertices' : vertices,
                 'distance' : face.dist,
                 'closest'  : [face.closest_point[i] for i in range(1,4)],
                 'hue'      : face.f_class.hue })
            face = face.next
        return faces

    def view(self):
        if PolyhedronViewer:
            self.viewer = PolyhedronViewer(
                self.face_list(),
                title='Dirichlet Domain of %s'%self.manifold_name)
        else:
            raise RuntimeError('The PolyhedronViewer class '
                               'was not imported.')

    def triangulation(self):
        """
        Returns the corresponding manifold as computed directly from
        the Dirichlet domain, regarded as polyhedron with faces
        identified in pairs.  Only works if this gives a manifold not
        an orbifold.
        """
        cdef c_Triangulation *c_manifold
        cdef Manifold M
        c_manifold = Dirichlet_to_triangulation(self.c_dirichlet_domain)
        if c_manifold is NULL:
            raise ValueError('The Dirichlet domain could not be '
                             'triangulated; perhaps this is an '
                             'orbifold group?')
        M = Manifold('empty')
        M.set_c_triangulation(c_manifold)
        M.set_name(self.manifold_name)
        return M

    def volume(self):
        """
        Returns the approximate volume of the DirichletDomain.
        Because matrices in O(3,1) tend to accumulate roundoff error,
        it's hard to get a good bound on the accuracy of the computed
        volume.  Nevertheless, the kernel computes the best value it
        can, with the hope that it will aid the user in recognizing
        manifolds defined by a set of generators.
        """
        return self.c_dirichlet_domain.approximate_volume
    


class DirichletDomain(CDirichletDomain):
    """
    A DirichletDomain object represents a Dirichlet Domain of 
    a hyperbolic manifold, typically centered at a point which
    is a local maximum of injectivity radius.  It will have ideal
    vertices if the manifold is not closed.
    
    Instantiate as M.dirichlet_domain() where M is a Manifold to
    obtain a Dirichlet Domain centered at a point which maximizes
    injectivity radius.
    
    Other options can be provided to customize the computation, with
    the default values shown here
    
    >>> M = Manifold('m003(3,-4)')
    >>> M.dirichlet_domain(vertex_epsilon=10.0**-8, displacement = [0.0, 0.0, 0.0], centroid_at_origin=True, maximize_injectivity_radius=True)
    40 finite vertices, 0 ideal vertices; 60 edges; 22 faces

    You can also create a Dirichlet Domain from a file listing matrix
    generators for the group, in SnapPea's "% Generator" format, via

       D = DirichletDomain(generator_file='test.gens')
    """
    pass

# Cusp Neighborhoods

cdef class CCuspNeighborhood:
    cdef c_CuspNeighborhoods *c_cusp_neighborhood
    cdef c_Triangulation *c_triangulation

    def __cinit__(self, Manifold manifold):
        if manifold.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        copy_triangulation(manifold.c_triangulation,
                           &self.c_triangulation)
        self.c_cusp_neighborhood = initialize_cusp_neighborhoods(
            self.c_triangulation)
        if self.c_cusp_neighborhood == NULL:
            raise RuntimeError('The cusp neighborhood construction failed.')
        self.manifold_name = manifold.name()

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_cusp_neighborhood != NULL:
            free_cusp_neighborhoods(self.c_cusp_neighborhood)

    def __repr__(self):
        N = self.num_cusps()
        return 'Cusp Neighborhood with %d cusp%s'%(
            N, N != 1 and 's' or '')

    def manifold(self):
        """
        Return a Manifold built from the current canonical triangulation.
        """
        cdef c_Triangulation *c_triangulation
        cdef Manifold M

        copy_triangulation(self.c_cusp_neighborhood.its_triangulation,
                           &c_triangulation)
        M = Manifold('empty')
        M.set_c_triangulation(c_triangulation)
        M.set_name(self.manifold_name + '_canonical')
        return M

    def check_index(self, which_cusp):
        N = int(which_cusp)
        if 0 <= N < self.num_cusps():
            return N
        else:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)
        
    def num_cusps(self):
        """
        Return the number of cusps.
        """
        return get_num_cusp_neighborhoods(self.c_cusp_neighborhood)

    def topology(self, which_cusp = 0):
        """
        Return the topological type of the specified cusp.
        """
        N = self.check_index(which_cusp)
        topology = get_cusp_neighborhood_topology(self.c_cusp_neighborhood,N)
        return CuspTopology[topology]

    def get_displacement(self, which_cusp = 0):
        """
        Return the displacement of the horospherical boundary of the
        specified cusp. The displacement is the hyperbolic distance
        that the horospherical boundary has been displaced from its
        "home" position, at which the area of the boundary is
        3sqrt(3)/8.  (The translates of all of the horospheres are
        guaranteed to be pairwise disjoint when each cusp has
        displacement 0.)
        """
        N = self.check_index(which_cusp)
        return get_cusp_neighborhood_displacement(self.c_cusp_neighborhood, N)

    def set_displacement(self, new_displacement, which_cusp=0):
        """
        Set the displacement of the specified cusp.
        """
        N = self.check_index(which_cusp)
        set_cusp_neighborhood_displacement(self.c_cusp_neighborhood,
                                           N, new_displacement)

    def stopping_displacement(self, which_cusp=0):
        """
        Return the displacement at which the specified cusp
        neighborhood bumps into itself or another cusp neighborhood.
        (Assumes the other displacements are fixed.)
        """
        return get_cusp_neighborhood_stopping_displacement(
            self.c_cusp_neighborhood, which_cusp)

    def stopper(self, which_cusp):
        """
        Return the index of the cusp which will be the first one that
        the specified cusp neighborhood bumps into.
        (Assumes the other displacements are fixed.)
        """
        return get_cusp_neighborhood_stopper_cusp_index(
            self.c_cusp_neighborhood, which_cusp)

    def reach(self, which_cusp=0):
        """
        Return the displacement at which the specified cusp
        neighborhood bumps into itself.  (This is twice the
        distance between nearest horoball lifts.)
        """
        return get_cusp_neighborhood_reach(
            self.c_cusp_neighborhood, which_cusp)

    def max_reach(self):
        """
        Return the maximum reach over all cusps.
        """
        return get_cusp_neighborhood_max_reach(
            self.c_cusp_neighborhood)

    def get_tie(self, which_cusp):
        """
        Return True if the specified cusp is a member of the tied group. 
        The displacements of the tied cusps are all the same.        
        """
        N = self.check_index(which_cusp)
        return get_cusp_neighborhood_tie(self.c_cusp_neighborhood, N)

    def set_tie(self, which_cusp, new_tie):
        """
        Mark the specified cusp as a member of the tied group. 
        """
        N = self.check_index(which_cusp)
        set_cusp_neighborhood_tie(self.c_cusp_neighborhood, N, new_tie)

    def volume(self, which_cusp=0):
        """
        Return the volume of the horoball neighborhood of the specified
        cusp.
        """
        N = self.check_index(which_cusp)
        return get_cusp_neighborhood_cusp_volume(self.c_cusp_neighborhood, N)
    
    def translations(self, which_cusp=0):
        """
        Return the (complex) Euclidean translations of the meridian
        and longitude of the specified cusp.
        """
        cdef Complex meridian
        cdef Complex longitude
        N = self.check_index(which_cusp)
        get_cusp_neighborhood_translations(self.c_cusp_neighborhood,
                                           N,
                                           &meridian,
                                           &longitude)
        return C2C(meridian), C2C(longitude)

    def horoballs(self, cutoff=0.1, which_cusp=0, full_list=True):
        """
        Return a list of dictionaries describing the horoballs with
        height at least cutoff.  The keys are 'center', 'radius', 'index'.
        """
        cdef CuspNbhdHoroballList* list
        cdef CuspNbhdHoroball ball
        list = get_cusp_neighborhood_horoballs(self.c_cusp_neighborhood,
                                                which_cusp,
                                                full_list,
                                                cutoff)
        if list == NULL:
            raise RuntimeError('The horoball construction failed.')
        result = []
        for n from 0 <= n < list.num_horoballs:
            ball = list.horoball[n]
            dict = {'center' : C2C(ball.center),
                    'radius' : ball.radius,
                    'index'  : ball.cusp_index}
            result.append(dict)
        free_cusp_neighborhood_horoball_list(list)
        return result

    def Ford_domain(self, which_cusp=0):
        """
        Return a list of pairs of complex numbers describing the
        endpoins of the segments obtained by projecting the edges of
        the Ford domain to the xy-plane in the upper half space model.
        """
        cdef CuspNbhdSegmentList* list
        cdef CuspNbhdSegment segment
        list = get_cusp_neighborhood_Ford_domain(self.c_cusp_neighborhood,
                                                 which_cusp)
        if list == NULL:
            raise RuntimeError('The Ford domain construction failed.')
        result = []
        for n from 0 <= n < list.num_segments:
            segment = list.segment[n]
            pair = ( C2C(segment.endpoint[0]), C2C(segment.endpoint[1]) )
            result.append(pair)
        free_cusp_neighborhood_segment_list(list)
        return result

    def triangulation(self, which_cusp=0):
        """
        Return a list of dictionaries describing the endpoints of the
        segments obtained by projecting the edges of the triangulation
        dual to the Ford domain into the xy-plane in the upper half
        space model.  The keys are 'endpoints' and 'indices'.
        """
        cdef CuspNbhdSegmentList* list
        cdef CuspNbhdSegment segment
        list = get_cusp_neighborhood_triangulation(self.c_cusp_neighborhood,
                                                   which_cusp)
        if list == NULL:
            raise RuntimeError('The triangulation construction failed.')
        result = []
        for n from 0 <= n < list.num_segments:
            segment = list.segment[n]
            endpoints = ( C2C(segment.endpoint[0]), C2C(segment.endpoint[1]) )
            indices = (segment.start_index,
                       segment.middle_index,
                       segment.end_index)
            result.append({'endpoints' : endpoints, 'indices' : indices})
        free_cusp_neighborhood_segment_list(list)
        return result

    def view(self, which_cusp=0, cutoff=None):
        """
        Create a 3D picture of the horoball packing.  One can specify
        which cusp to put at infinity and how large of horoballs to
        look at, e.g.

        >>> M = Manifold('m125')
        >>> C = M.cusp_neighborhood()
        >>> C.view(which_cusp = 1, cutoff=0.2)
        """
        if HoroballViewer:
            self.viewer = HoroballViewer(
                self, which_cusp=which_cusp, cutoff=cutoff,
                title='Cusp neighborhood #%s of %s'%(
                    which_cusp,
                    self.manifold_name
                    ))
        else:
            raise RuntimeError('The HoroballViewer class was not imported.')
        
class CuspNeighborhood(CCuspNeighborhood):
    """
    A CuspNeighborhood object represents an equivariant collection of
    disjoint horoballs that project to cusp neighborhoods.

    Instantiate as M.cusp_neighborhood()
    """
    pass

#  Symmetry_group

cdef class SymmetryGroup:
    """
    A SymmetryGroup is a group of self-isometries of hyperbolic
    3-manifold.  Instantiate as follows:

    >>> M = Manifold('m004')
    >>> M.symmetry_group()
    D4
    """
    cdef c_SymmetryGroup *c_symmetry_group
    cdef readonly _is_full_group
    cdef readonly _owns_c_symmetry_group
    
    def __cinit__(self, is_full_group, owns_c_symmetry_group):
        self.c_symmetry_group = NULL 
        self._is_full_group = is_full_group
        self._owns_c_symmetry_group = owns_c_symmetry_group

    def __dealloc__(self):
        if self._owns_c_symmetry_group:
            free_symmetry_group(self.c_symmetry_group)

    cdef _set_c_symmetry_group(self, c_SymmetryGroup * c_symmetry_group):
        if c_symmetry_group is NULL:
            raise ValueError('You tried to create an *empty* SymmetryGroup.')
        self.c_symmetry_group = c_symmetry_group

    def is_full_group(self):
        """
        Return whether the full symmetry group has been found.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_full_group()
        True
        """
        return self._is_full_group

    def __repr__(self):
        if self.is_full_group():
            thePretext = ''
        else:
            thePretext = 'at least '

        if self.is_abelian():
            theText = repr(self.abelian_description())
        elif self.is_dihedral():
            theText = 'D%d'%(self.order()//2)
        elif self.is_polyhedral():
            theText = self.polyhedral_description()
        elif self.is_S5():
            theText = 'S5'
        elif self.is_direct_product():
            theText =     '%s x %s' % self.direct_product_description()
        else:
            theText = 'nonabelian group of order %d'%self.order()
        
        return thePretext + theText
    
    def order(self):
        """
        Return the order of the symmetry group

        >>> S = Manifold('s000').symmetry_group()
        >>> S.order()
        4
        """
        return symmetry_group_order(self.c_symmetry_group)
    
    def is_abelian(self):
        """
        Return whether the symmetry group is abelian.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_abelian()
        False
        """
        cdef c_AbelianGroup* abelian_description = NULL
        ans = B2B(symmetry_group_is_abelian(
                self.c_symmetry_group, &abelian_description))
        return ans

    def abelian_description(self):
        """
        If the symmetry group is abelian, return it as an AbelianGroup

        >>> S = Manifold('v3379').symmetry_group()
        >>> S.abelian_description()
        Z/2 + Z/2 + Z/2
        """
        cdef c_AbelianGroup* A
        cdef int n
        is_abelian = B2B(symmetry_group_is_abelian(self.c_symmetry_group, &A))
        if not is_abelian:
            raise ValueError('The symmetry group is not abelian.')

        coeffs = []
        for n from 0 <= n < A.num_torsion_coefficients:
                coeffs.append(A.torsion_coefficients[n])

        # Don't need to free A as it is attached to the symmetry group object
        return AbelianGroup(elementary_divisors=coeffs)
            
    
    def is_dihedral(self):
        """
        Return whether the symmetry group is dihedral.
        
        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_dihedral()
        True
        """
        return B2B(symmetry_group_is_dihedral(self.c_symmetry_group))
    
    def is_polyhedral(self):
        """
        Returns whether the symmetry group is a (possibly binary)
        polyhedral group.
        """
        return B2B(symmetry_group_is_polyhedral(self.c_symmetry_group,
                                                NULL, NULL, NULL, NULL))
    
    def polyhedral_description(self):
        """
        If the symmetry group is a (possibly binary)
        polyhedral group, return a description of it.  
        """
        cdef Boolean is_binary_group
        cdef int p,q,r
        
        if not self.is_polyhedral():
            raise ValueError('The symmetry group is not polyhedral.')

        symmetry_group_is_polyhedral(self.c_symmetry_group,
                                     &is_binary_group, &p, &q, &r)

        assert p == 2

        if q == 2:
            assert(is_binary_group)
            name = 'binary dihedral group <2,2,%d>' % r
        elif q== 3:
            name = 'binary ' if is_binary_group else ''
            name += {3: 'tetrahedral group',
                     4: 'octahedral group',
                     5: 'icosohedral group'}[r]

        return name 
    
    def is_S5(self):
        """
        Returns whether the group is the symmetric group on five things.  
        """
        return B2B(symmetry_group_is_S5(self.c_symmetry_group))
    
    def is_direct_product(self):
        """
        Return whether the SymmetryGroup is a nontrivial direct
        product with at least one nonabelian factor.  
        
        >>> S = Manifold('s960').symmetry_group()
        >>> S.is_direct_product()
        True
        >>> S
        Z/4 x D3
        """
        return B2B(symmetry_group_is_direct_product(self.c_symmetry_group))
    
    def direct_product_description(self):
        """
        If the SymmetryGroup is a nontrivial direct product with at
        least one nonabelian factor, return a pair of SymmetryGroups
        consisting of the (two) factors.

        >>> S = Manifold('s960').symmetry_group()
        >>> S.direct_product_description()
        (Z/4, D3)
        """
        if not self.is_direct_product():
            raise ValueError('The symmetry group is not a nontrivial, '
                             'nonabelian direct product.')

        cdef c_SymmetryGroup* c_factor_0
        cdef c_SymmetryGroup* c_factor_1
        cdef SymmetryGroup factor_0
        cdef SymmetryGroup factor_1
        
        c_factor_0 = get_symmetry_group_factor(self.c_symmetry_group, 0)
        c_factor_1 = get_symmetry_group_factor(self.c_symmetry_group, 1)
        
        factor_0 = SymmetryGroup(True, False)
        factor_1 = SymmetryGroup(True, False)
        factor_0._set_c_symmetry_group(c_factor_0)
        factor_1._set_c_symmetry_group(c_factor_1)
        return (factor_0, factor_1)
    
    def is_amphicheiral(self):
        """
        Return whether the manifold has an orientation reversing symmetry.

        >>> S = Manifold('m004').symmetry_group()
        >>> S.is_amphicheiral()
        True
        """        
        return B2B(symmetry_group_is_amphicheiral(self.c_symmetry_group))
    
    def is_invertible_knot(self):
        """
        Return whether a one-cusped has a symmetry that acts on the
        cusp via the matrix -I.

        >>> S = Manifold('m015').symmetry_group()
        >>> S.is_invertible_knot()
        True
        """
        return B2B(symmetry_group_invertible_knot(self.c_symmetry_group))
        
    def commutator_subgroup(self):
        """
        Return the commutator subgroup of the SymmetryGroup

        >>> S = Manifold('m004').symmetry_group()
        >>> S
        D4
        >>> S.commutator_subgroup()
        Z/2
        """
        cdef c_SymmetryGroup* c_comm_subgroup
        cdef SymmetryGroup comm_subgroup

        c_comm_subgroup = get_commutator_subgroup(self.c_symmetry_group)
        comm_subgroup = SymmetryGroup(self.is_full_group(), True)
        comm_subgroup._set_c_symmetry_group(c_comm_subgroup)
        return comm_subgroup

    def abelianization(self):
        """
        Return the abelianization of the symmetry group

        >>> S = Manifold('m004').symmetry_group()
        >>> S.abelianization()
        Z/2 + Z/2
        """
        
        if not self.is_full_group():
            raise ValueError('The full symmetry group is not known.')

        cdef c_SymmetryGroup* c_abelianization
        cdef SymmetryGroup abelianization

        c_abelianization = get_abelianization(self.c_symmetry_group)
        abelianization = SymmetryGroup(self.is_full_group(), True)
        abelianization._set_c_symmetry_group(c_abelianization)
        return abelianization.abelian_description()

    def center(self):
        """
        Return the center of the symmetry group

        >>> S = Manifold('m004').symmetry_group()
        >>> S.center()
        Z/2
        """
        
        if not self.is_full_group():
            raise ValueError('The full symmetry group not known.')

        cdef c_SymmetryGroup* c_center
        cdef SymmetryGroup center

        c_center = get_center(self.c_symmetry_group)
        center = SymmetryGroup(self.is_full_group(), True)
        center._set_c_symmetry_group(c_center)
        return center.abelian_description()

    def multiply_elements(self, i, j):
        """
        Returns the product of group elements i and j.  The convention
        is that products of symmetries read right to left.  That is,
        the composition (symmetry[i] o symmetry[j]) acts by first
        doing symmetry[j], then symmetry[i].

        >>> S = Manifold('m004').symmetry_group()
        >>> S.multiply_elements(2, 3)
        1
        """
        cdef int prod
        order = self.order()
        for x in [i,j]:
            if not (0 <= x < order):
                raise ValueError('The symmetry group has only %d '
                                 'elements.' % order)

        return symmetry_group_product(self.c_symmetry_group, i, j)

    def isometries(self):
        """
        Return a detailed list of all the isometries in the symmetry group.

        >>> S = Manifold('s959').symmetry_group()
        >>> isoms = S.isometries()
        >>> isoms[8]
        0 -> 1   1 -> 0 
        [-1 -1]  [ 0  1]
        [ 1  0]  [-1 -1]
        Does not extend to link
        """
        cdef IsometryList *isometries

        isometries = get_symmetry_list(self.c_symmetry_group)
        ans = IsometryListToIsometries(isometries)
        return ans

# get_triangulation

split_filling_info = re.compile('(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))+)')
is_census_manifold = re.compile('([msvtxy])([0-9]+)$')
is_torus_bundle = re.compile('b([+-no])([+-])([lLrR]+)$')
is_knot_complement = re.compile('([0-9]+_[0-9]+)$')
is_link_complement1 = re.compile('(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$')
is_link_complement2 = re.compile('(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$')
is_link_complement3 = re.compile('[lL](?P<components>[0-9]{1})(?P<crossings>[0-9]{2})(?P<index>[0-9]+)$')
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)$')
is_HT_link = re.compile('[KL][0-9]+[an]([0-9]+)$')
is_braid_complement = re.compile('braid(\[[0-9, -]+\])$')
is_int_DT_exterior = re.compile('DT(\[[0-9, -]+\])$')
is_alpha_DT_exterior = re.compile('DT\[([a-zA-Z]+)\]$')
is_census_knot = re.compile('[kK][2-7]_([0-9]+)$')
is_twister_bundle = re.compile('Bundle\(S_\{(\d+),(\d+)\},\[*([abcABC_\d!,*]*)\]*\)')
is_twister_splitting = re.compile('Splitting\(S_\{(\d+),(\d+)\},\[*([abcABC_\d!,*]*)\]*,*\[*([abcABC_\d!,*]*)\]*\)')

def bundle_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_bundle.match(desc)
    if m:
        g, n, monodromy = m.groups()
        g, n = int(g), int(n)
        monodromy = monodromy.replace(',', '*')
        return twister.twister(surface=(g,n), monodromy=monodromy,with_hyperbolic_structure=True)
    
def splitting_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_splitting.match(desc)
    if m:
        g, n, gluing, handles = m.groups()
        g, n = int(g), int(n)
        gluing, handles = gluing.replace(',', '*'), handles.replace(',', '*')
        return twister.twister(surface=(g,n),gluing=gluing, handles=handles,with_hyperbolic_structure=True)

#Orientability.orientable = 0
rev_spec_dict = {(5, 0) : 'm',
                 (5, 1) : 'm',
                 (6, 0) : 's',
                 (7, 0) : 'v',
                 (6, 1) : 'x',
                 (7, 1) : 'y'}

triangulation_help =  """ 
    A %s is specified by a string, according to the
    conventions detailed in its docstring.  
    """

cdef c_Triangulation* triangulation_from_bytes(bytestring) except ? NULL:
    cdef c_Triangulation* c_triangulation = NULL
    cdef TerseTriangulation c_terse
    cdef int N=0, n=0, m=1
    byteseq = bytearray(bytestring)
    c_terse.num_tetrahedra = N = byteseq[0]
    c_terse.glues_to_old_tet = <Boolean*>malloc(2*N*sizeof(Boolean))
    c_terse.which_old_tet = <int*>malloc((N+1)*sizeof(int))
    c_terse.which_gluing = <Permutation*>malloc((N+1)*sizeof(Permutation))
    bit = 0
    for n from 0 <= n < 2*N:
        if byteseq[m] & (1 << bit):
            c_terse.glues_to_old_tet[n] = 1
        else:
            c_terse.glues_to_old_tet[n] = 0
        bit += 1
        if bit%8 == 0:
            m += 1
            bit = 0
    if bit:
        m += 1
    for n from 0 <= n < 1 + N:
        c_terse.which_old_tet[n] = int(byteseq[m])
        m += 1
    for n from 0 <= n < 1 + N:
        c_terse.which_gluing[n] = int(byteseq[m])
        m += 1
    c_triangulation = terse_to_tri(&c_terse)
    free(c_terse.glues_to_old_tet)
    free(c_terse.which_old_tet)
    free(c_terse.which_gluing)
    return c_triangulation

cdef c_Triangulation* triangulation_from_database(data_object, name) except ? NULL:
    cdef c_Triangulation* c_triangulation=NULL
    cdef char* c_name
    cdef int n
    use_string, cobs, perm, bytestring = data_object(name)
    if use_string:
        c_triangulation = read_triangulation_from_string(bytestring)
    else:
        c_triangulation = triangulation_from_bytes(bytestring)
        if cobs:
            n = len(cobs)
            matrices = <MatrixInt22 *>malloc(n*sizeof(MatrixInt22))
            install_combinatorial_bases(c_triangulation, matrices)
            for i in range(n):
                matrices[i][0][0]=cobs[i][0][0]
                matrices[i][0][1]=cobs[i][0][1]
                matrices[i][1][0]=cobs[i][1][0]
                matrices[i][1][1]=cobs[i][1][1]
            change_peripheral_curves(c_triangulation, matrices)
            if perm:
                indices = <int *>malloc(n*sizeof(int))
                for i in range(n):
                    indices[i] = perm[i]
                reindex_cusps(c_triangulation, indices)
                free(indices)
            free(matrices)
    c_name = b_name = to_byte_str(name)
    set_triangulation_name(c_triangulation, c_name)
    return c_triangulation

cdef c_Triangulation* get_triangulation(spec) except ? NULL:
    cdef c_Triangulation* c_triangulation = NULL
    cdef LRFactorization* gluing
    cdef char* LRstring
    cdef char* c_name
    cdef int LRlength, i, n
    cdef Triangulation T

    # Step -1 Check for an entire-triangulation-file-in-a-string
    if spec.startswith('% Triangulation'):
        b_spec = to_byte_str(spec)
        return read_triangulation_from_string(b_spec)

    # get filling info, if any
    m = split_filling_info.match(spec)
    if m:
        real_name = m.group(1)
        fillings = re.subn('\)\(', '),(', m.group(2))[0]
        fillings = eval( '[' + fillings + ']', {})
    else:
        real_name = spec
        fillings = ()

    # Step 1. Check for a census manifold
    m = is_census_manifold.match(real_name)
    if m:
        c_triangulation = triangulation_from_database(
            database.CuspedManifoldData, real_name)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

     # Step 2. Check for a punctured torus bundle 
    m = is_torus_bundle.match(real_name)
    if m:
        LRpart = m.group(3).upper()
        LRlength = len(LRpart)
        LRstring = LRpart
        negative_determinant = negative_trace = 0

        if m.group(1) == '-' or m.group(1) == 'n':
            negative_determinant = 1
            
        if m.group(2) == '+':
            negative_trace = 0
        else:
            negative_trace = 1
        gluing = alloc_LR_factorization(LRlength)
        gluing.is_available = True
        gluing.negative_determinant = negative_determinant
        gluing.negative_trace = negative_trace
        for i from 0 <= i < LRlength:
           gluing.LR_factors[i] = LRstring[i] 
        c_triangulation =  triangulate_punctured_torus_bundle(gluing)
        free_LR_factorization(gluing)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 3. Check for a Rolfsen link complement
    name = None
    m = is_knot_complement.match(real_name)
    if m:
        name = real_name
    for regex in [is_link_complement1,
                  is_link_complement2,
                  is_link_complement3]:
        m = regex.match(real_name)
        if m:
            if int(m.group('components')) > 1:
                name = '%d^%d_%d' % (int(m.group('crossings')),
                                     int(m.group('components')),
                                     int(m.group('index')))
            else:
                name = '%d_%d' % (int(m.group('crossings')),
                                  int(m.group('index')))
            
            break
    if name:
        try:
            c_triangulation = triangulation_from_database(
                database.LinkExteriorData, name)
        except: 
            raise KeyError('The link complement %s was not found.'%
                           real_name)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 4. Check for an HT link.
    m = is_HT_link.match(real_name)
    if m:
        try:
            c_triangulation = triangulation_from_database(
                database.HTLinkExteriorData, real_name)
        except: 
            raise IOError('The HT link %s was not found.'%real_name)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 5. Check for a Hoste-Thistlethwaite knot.
    m = is_HT_knot.match(real_name)
    if m:
        c_triangulation = get_HT_knot(int(m.group('crossings')),
                             m.group('alternation'),
                             int(m.group('index')))
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 6. Check for a census knot.
    m = is_census_knot.match(real_name)
    if m:
        try:
            c_triangulation = triangulation_from_database(
                database.CensusKnotData, real_name)
        except: 
            raise IOError('The census knot %s was not found.'%real_name)
        set_cusps(c_triangulation, fillings)
        return c_triangulation
        
    # Step 7. See if a (fibered) braid complement is requested

    m = is_braid_complement.match(real_name)
    if m:
        word = eval(m.group(1), {})
        num_strands = max([abs(x) for x in word]) + 1
        c_triangulation = get_fibered_manifold_associated_to_braid(
            num_strands, word)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 8.  See if a knot exterior is requested via its
    # Dowker-Thistlethwaite code:

    m = is_int_DT_exterior.match(real_name)
    if m:
        word = eval(m.group(1), {})
        c_triangulation = get_link_exterior_from_DT(word)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    m = is_alpha_DT_exterior.match(real_name)
    if m:
        c_triangulation = get_link_exterior_from_alpha_DT(m.group(1))
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 9.  See if a bundle or splitting is given in Twister's notation

    shortened_name = real_name.replace(' ', '')
    mb = is_twister_bundle.match(shortened_name)
    ms = is_twister_splitting.match(shortened_name)
    if mb or ms:
        func = bundle_from_string if mb else splitting_from_string
        T = func(shortened_name)
        copy_triangulation(T.c_triangulation, &c_triangulation)
        return c_triangulation

    # Step 10. If all else fails, try to load a manifold from a file.
    try:
        locations = [os.curdir, os.environ['SNAPPEA_MANIFOLD_DIRECTORY']]
    except KeyError:
        locations = [os.curdir]
    found = 0
    for location in locations:
        pathname = os.path.join(location, real_name)
        if os.path.isfile(pathname):
            file = open(pathname, 'r')
            first_line = file.readline()[:-1]
            file.close()
            if first_line.find('% Link Projection') > -1:
                c_triangulation = triangulate_link_complement_from_file(
                    pathname, '')
            else:
                c_triangulation = read_triangulation(pathname)
            set_cusps(c_triangulation, fillings)
            return c_triangulation

    # Step 9. Give up.
    raise IOError('The manifold file %s was not found.\n%s'%
                  (real_name, triangulation_help%
                   'Triangulation or Manifold'))
        
cdef int set_cusps(c_Triangulation* c_triangulation, fillings) except -1:
    if c_triangulation == NULL:
        return 0
    if len(fillings) > 0:
        num_cusps = get_num_cusps(c_triangulation) 
        if len(fillings) > num_cusps:
            raise ValueError('The number of fillings specified exceeds '
                             'the number of cusps.')
        for i in range(len(fillings)):
            meridian, longitude = fillings[i]
            is_complete = (meridian == 0 and longitude == 0)
            set_cusp_info(c_triangulation, i,
                          is_complete, meridian, longitude)
    return 0

# Support for Hoste-Thistethwaite tables

# These dictionaries are used in accessing the tables.  The key is the
# number of crossings; the value is the number of knots with that many
# crossings.

Alternating_numbers = { 3:1, 4:1, 5:2, 6:3, 7:7, 8:18, 9:41, 10:123, 11:367,
                        12:1288, 13:4878, 14:19536, 15:85263, 16:379799 }

Nonalternating_numbers = { 8:3, 9:8, 10:42, 11:185, 12:888, 13:5110,
                           14:27436, 15:168030, 16:1008906 }

Alternating_offsets = {}
offset = 0
for i in range(3,17):
    Alternating_offsets[i] = offset
    offset +=  Alternating_numbers[i]
Num_Alternating = offset

Nonalternating_offsets = {}
offset = 0
for i in range(8,17):
    Nonalternating_offsets[i] = offset
    offset += Nonalternating_numbers[i]
Num_Nonalternating = offset

def extract_HT_knot(record, crossings, alternation):
    DT=[]
    size = (1+crossings)//2
    for byte in record[:size]:
        first_nybble = (byte & 0xf0) >> 4
        if first_nybble == 0: first_nybble = 16
        DT.append(2*first_nybble)
        second_nybble = byte & 0x0f
        if second_nybble == 0: second_nybble = 16
        DT.append(2*second_nybble)
    if alternation == 'n':
        signs = record[-2]<<8 | record[-1]
        mask = 0x8000
        for i in range(crossings):
            if (signs & (mask >> i)) == 0:
                DT[i] = -DT[i]
    return DT[:crossings]

def get_HT_knot_DT(crossings, alternation, index):
    size = (1 + crossings)//2
    index -= 1
    if ( alternation == 'a'
         and crossings in Alternating_numbers.keys()
         and 0 <= index < Alternating_numbers[crossings] ):
        offset = 8*(Alternating_offsets[crossings] +  index)
        Alternating_table.seek(offset)
        data = Alternating_table.read(size)
        record = struct.unpack('%dB'%size, data)
    elif ( alternation == 'n'
         and crossings in Nonalternating_numbers.keys()
         and 0 <= index < Nonalternating_numbers[crossings] ):
        offset = 10*(Nonalternating_offsets[crossings] +  index)
        Nonalternating_table.seek(offset)
        data = Nonalternating_table.read(size+2)
        record = struct.unpack('%dB'%(size+2), data)
    else:
        raise ValueError('You have specified a Hoste-Thistlethwaite '
                         'knot with an \n'
                         'inappropriate index or number of crossings.')

    DT = extract_HT_knot(record, crossings, alternation)
    return DT

cdef c_Triangulation* get_link_exterior_from_DT(DT) except ? NULL:
    cdef int* DT_array
    cdef int i
    cdef c_Triangulation* c_triangulation
    DT_array = <int*>malloc(len(DT)*sizeof(int))
    for i from 0 <= i < len(DT):
        DT_array[i] = DT[i]
    c_triangulation = DT_int_to_triangulation(len(DT), DT_array)
    name = to_byte_str('DT'+repr(DT))
    set_triangulation_name(c_triangulation, name)
    free(DT_array)
    return c_triangulation


def DT_alpha_to_int(x):
        return string.ascii_lowercase.index(x) + 1

cdef c_Triangulation* get_link_exterior_from_alpha_DT(DT) except ? NULL:
    """
    Load the link exterior specified by the alpha DT code in the
    extended Snap DT style.
    The format is:

    Creates a link complement from a Dowker-Thistlethwaite code.
    For knots this is just the Dowker code preceded by
    <crossings>a<crossings> where <crossings> is a single letter
    code for a number "a=1,b=2...". E.g. figure 8 knot is "dadbcda".
    More generally we have:

    <num-crossings><num-cpts>
    <num-cross-cpt-1><num-cross-cpt-2>...<num-cross-cpt-n>
    <cpt-1-code><cpt-2-code>...<cpt-n-code>

    I'm guessing that this code was only tested with (or only written
    to work with) a certain subset of DT codes.  Looking at the DT
    codes supplied with Morwen's table of links, I think that it
    assumes at least the following, in the context of 2-crossing
    links:

    There is a crossing (1, n + 2) where n = (# of crossings on the
    first comp)

    More broadly, the code may assume that the input is
    lexiographically first among some class of such DT codes.
    """    
    cdef c_Triangulation* c_triangulation
    cdef char* c_DT
    DTbytes = bytes(DT.encode('ascii'))
    c_DT = DTbytes
    
    # Let's do a rudimentary check that the DT code is well-formed.
    crossings = DT_alpha_to_int(DT[0])
    components = DT_alpha_to_int(DT[1])
    if (len(DT) != 2 + components + crossings
        or sum(map(DT_alpha_to_int, DT[2:2+components])) != crossings):
        raise ValueError('The DT string %s is not well-formed.'  % DT)
    rest = list(DT[2 + components:].lower())
    rest.sort()
    if ''.join(rest) != string.ascii_lowercase[:crossings]:
        raise ValueError('The DT string %s is not well-formed.'  % DT)

    c_triangulation = DT2Triangulation(c_DT)
    if c_triangulation == NULL:
        raise ValueError("The DT string %s doesn't seem to be "
                         "realizable." % DT)
    name = to_byte_str('DT['+ DT + ']')
    set_triangulation_name(c_triangulation, name)
    return c_triangulation


cdef c_Triangulation* get_HT_knot(crossings, alternation, index) except ? NULL:
    cdef c_Triangulation* c_triangulation
    DT = get_HT_knot_DT(crossings, alternation, index)
    c_triangulation = get_link_exterior_from_DT(DT)
    name = to_byte_str('%d'%crossings + alternation + '%d'%index)
    set_triangulation_name(c_triangulation, name)
    return c_triangulation

def get_HT_knot_by_index(alternation, index):
    DT=[]
    crossings = 16
    if alternation == 'a':
        for i in range(3,17):
            if Alternating_offsets[i] > index:
                crossings = i-1
                break
        index_within_crossings = index - Alternating_offsets[crossings]
    elif alternation == 'n':
        for i in range(8,17):
            if Nonalternating_offsets[i] > index:
                crossings = i-1
                break
        index_within_crossings = index - Nonalternating_offsets[crossings]
    name = '%d'%crossings + alternation + '%d'%(index_within_crossings + 1)
    return Manifold(name)

#   Iterators

class Census:
    """
    Base class for manifold Iterators/Sequences.
    """
    # subclasses redefine this
    length = 0

    def __init__(self, indices=(0,0,0)):
        myslice = slice(*indices)
        self.start, self.stop, self.step = myslice.indices(self.length)
        self.index = self.start

    def __iter__(self):
        return self

    def __contains__(self, item):
        raise NotImplementedError("This census does not support manifold lookup")

    def next(self):
        if self.index >= self.stop:
            raise StopIteration
        self.index = self.index + self.step
        return self[self.index-self.step]

    __next__ = next

    def __len__(self):
        return self.length
    
    # subclasses override this
    def __getitem__(self, n):
        pass

#  Cusped Census

Orientable_lengths = (301, 962, 3552, 12846, 301+962+3552+12846)
Nonorientable_lengths = (114, 259, 887, 114+259+887) 

five_tet_orientable = (
   3, 4, 6, 7, 9, 10, 11, 15, 16, 17, 19, 22, 23, 26, 27, 29, 30, 32, 33,
   34, 35, 36, 37, 38, 39, 40, 43, 44, 45, 46, 47, 49, 52, 53, 54, 55, 58,
   59, 60, 61, 62, 64, 66, 67, 69, 70, 71, 72, 73, 74, 76, 77, 78, 79, 80,
   81, 82, 83, 84, 85, 87, 89, 90, 93, 94, 95, 96, 98, 99, 100, 102, 103,
   104, 105, 106, 108, 110, 111, 112, 113, 115, 116, 117, 118, 119, 120,
   121, 122, 123, 125, 129, 130, 135, 136, 137, 139, 140, 141, 142, 143,
   144, 145, 146, 147, 148, 149, 150, 151, 154, 155, 157, 159, 160, 161,
   162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 175, 178,
   179, 180, 181, 182, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193,
   194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 206, 207, 208, 209,
   210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
   224, 225, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 238, 239,
   240, 241, 242, 243, 244, 246, 247, 248, 249, 250, 251, 252, 253, 255,
   256, 257, 258, 259, 260, 261, 262, 263, 266, 267, 269, 270, 271, 272,
   273, 275, 276, 277, 278, 279, 280, 281, 282, 284, 285, 286, 287, 288,
   289, 290, 291, 292, 293, 294, 295, 296, 297, 299, 300, 302, 303, 304,
   305, 306, 307, 310, 312, 316, 318, 319, 320, 322, 326, 328, 329, 335,
   336, 337, 338, 339, 340, 342, 343, 345, 346, 349, 350, 351, 352, 356,
   357, 358, 359, 360, 361, 363, 364, 366, 367, 368, 369, 370, 371, 372,
   373, 374, 376, 378, 385, 388, 389, 390, 391, 392, 393, 395, 397, 400,
   401, 402, 403, 410, 412)

five_tet_nonorientable = (
   0, 1, 2, 5, 8, 12, 13, 14, 18, 20, 21, 24, 25, 28, 31, 41, 42, 48, 50,
   51, 56, 57, 63, 65, 68, 75, 86, 88, 91, 92, 97, 101, 107, 109, 114, 124,
   126, 127, 128, 131, 132, 133, 134, 138, 152, 153, 156, 158, 174, 176,
   177, 187, 204, 205, 226, 237, 245, 254, 264, 265, 268, 274, 283, 298,
   301, 308, 309, 311, 313, 314, 315, 317, 321, 323, 324, 325, 327, 330,
   331, 332, 333, 334, 341, 344, 347, 348, 353, 354, 355, 362, 365, 375,
   377, 379, 380, 381, 382, 383, 384, 386, 387, 394, 396, 398, 399, 404,
   405, 406, 407, 408, 409, 411, 413, 414)

class CuspedCensus(Census):
    """
    Base class for Iterators/Sequences for manifolds in the SnapPea
    Cusped Census.
    """
    five_length, six_length, seven_length, eight_length, length = Orientable_lengths
    orientability = Orientability.index('orientable')

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)
        self.Census_Morwen8 = tarfile.open(os.path.join(manifold_path, 'morwen8.tgz'), 'r:*')

    # Override
    def lookup(self, n):
        return five_tet_orientable[n]

    def __getitem__(self, n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        if n < 0:
            n = self.length - n
        if n < self.five_length:
            num_tet = 5
            census_index = self.lookup(n)
        elif n - self.five_length < self.six_length:
            num_tet = 6
            census_index = n - self.five_length
        elif n - self.five_length - self.six_length < self.seven_length:
            num_tet = 7
            census_index = n - self.five_length - self.six_length
        elif (self.orientability == Orientability.index('orientable')
              and (n - self.five_length - self.six_length -
                   self.seven_length < self.eight_length)):
              census_index = (n - self.five_length - self.six_length
                              - self.seven_length)
              # Make this work without passing the spec to Manifold()
              num = repr(census_index)
              spec =  "t" + "0"*(5 - len(num)) + num
              tarpath = "morwen8/" + spec
              try:
                  filedata = self.Census_Morwen8.extractfile(tarpath).read()
                  c_triangulation = read_triangulation_from_string(filedata)
              except: 
                  raise IOError('The Morwen 8 tetrahedra manifold %s '
                                'was not found.'% spec)
              result = Manifold(spec='empty')
              result.set_c_triangulation(c_triangulation)
              return result              
              ###return Manifold('t%d' % census_index)
        else:
            raise IndexError('Index is out of range.')
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, num_tet, self.orientability, census_index)
        if c_triangulation == NULL:
            print(num_tet, census_index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        return result

class ObsOrientableCuspedCensus(CuspedCensus):
    """
    Obsolete.
    """

class ObsNonorientableCuspedCensus(CuspedCensus):
    """
    Obsolete.
    """
    five_length, six_length, seven_length, length = Nonorientable_lengths
    orientability = Orientability.index('nonorientable')

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)

    def lookup(self, n):
        return five_tet_nonorientable[n]

# Closed Census (Obsolete)

class ObsOrientableClosedCensus(Census):
    """
    Obsolete.
    """
    data = None
    orientability = Orientability.index('orientable')

    def __init__(self, indices=(0,11031,1)):
        if ObsOrientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedOrientableDistinct.txt')
            closed_orientable = open(datafile)
            ObsOrientableClosedCensus.data = closed_orientable.readlines()
            closed_orientable.close()
        self.length = len(ObsOrientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = ObsOrientableClosedCensus.data[n].split()
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill(( int(m),int(l)) )
        return result

class ObsNonorientableClosedCensus(Census):
    """
    Obsolete.
    """
    data = None
    orientability = Orientability.index('nonorientable')
    
    def __init__(self, indices=(0,17,1)):
        if ObsNonorientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedNonorientableDistinct.txt')
            closed_nonorientable = open(datafile)
            ObsNonorientableClosedCensus.data = closed_nonorientable.readlines()
            closed_nonorientable.close()
        self.length = len(ObsNonorientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        cdef c_Triangulation* c_triangulation
        cdef Manifold result
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = ObsNonorientableClosedCensus.data[n].split()
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill( (int(m),int(l)) )
        return result

# Knot tables

class KnotExteriors(Census):
    """
    Base class for Iterators/Sequences for knots from the
    Hoste-Thistlethwaite tables.
    """
    length = sum(Alternating_numbers.values())
    alternation = 'a'

    def __init__(self, indices=(0, sum(Alternating_numbers.values()), 1)):
        Census.__init__(self, indices)

    def __getitem__(self, n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        else:
            return get_HT_knot_by_index(self.alternation, n)

class AlternatingKnotExteriors(KnotExteriors):
    """
    Iterator/Sequence for Alternating knot exteriors from the
    Hoste-Thistlethwaite tables.   Goes through 16 crossings. 
    """

class NonalternatingKnotExteriors(KnotExteriors):
    """
    Iterator/Sequence for nonAlternating knot exteriors from the
    Hoste-Thistlethwaite tables.  Goes through 16 crossings. 
    """
    length = sum(Nonalternating_numbers.values())
    alternation = 'n'

    def __init__(self, indices=(0, sum(Nonalternating_numbers.values()), 1)):
        Census.__init__(self, indices)

census_knot_numbers = [0, 0, 1, 2, 4, 22, 43, 129]

class ObsCensusKnots(Census):
    """
    Obsolete
    """
    length = sum(census_knot_numbers)

    def __init__(self, indices=(0, sum(census_knot_numbers), 1)):
        Census.__init__(self, indices)
        self.Census_Knots = tarfile.open(census_knot_archive, 'r:*')

    def __repr__(self):
        return 'Knots in S^3 which appear in the SnapPea Census'
    
    def __getitem__(self, n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        else:
            total = 0
            for m in range(2,8):
                if total + census_knot_numbers[m] <= n :
                    total += census_knot_numbers[m]
                    continue
                else:
                    name = 'K%s_%s'%(m, n - total + 1)
                    break
            if name:
                tarpath =  'CensusKnots/%s'%name
                try:
                    filedata = self.Census_Knots.extractfile(tarpath).read()
                    c_triangulation = read_triangulation_from_string(filedata)
                except: 
                    raise IOError, "The census knot %s was not found."%name
                result =  Manifold('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result              
            else:
                raise IndexError('There are only 201 census knots.')

class ObsLinkExteriors(Census):
    """
    Obsolete.
    """
    # num_links[component][crossings] is the number of links with
    # specified number of components and crossings.
    num_links= [ [], [0, 0, 0, 1, 1, 2, 3, 7, 21, 49, 166, 552],
                 [0, 0, 0, 0, 1, 1, 3, 8, 16, 61, 183],
                 [0, 0, 0, 0, 0, 0, 3, 1, 10, 21, 74],
                 [0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 24],
                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]

    max_crossings = 11

    def __init__(self, components, indices=(0,10000,1)):
         self.Christy_links = tarfile.open(link_archive, 'r:*')

         if not (1 <= components < len(self.num_links) ):
            raise IndexError('SnapPy has no data on links with '
                             '%s components.' % components)

         self.components = components

         self.length = sum(self.num_links[components])

         Census.__init__(self, (indices[0],
                                min(self.length, indices[1]),
                                indices[2]))

    def __repr__(self):
        return ('Christy census of %s-component link complements '
                'in S^3' % self.components)

    def __getitem__(self,j):
        if isinstance(j, slice):
            return self.__class__(j.indices(self.length))
        so_far = 0
        for k in range(self.max_crossings + 1):
            n =  self.num_links[self.components][k]
            so_far = so_far + n
            if so_far > j:
                l = j - so_far + n + 1
                filename = 'L%d%.2d%.3d' % (self.components, k, l)
                if self.components > 1:
                    name = "%d^%d_%d" % (k, self.components, l)
                else:        
                    name = "%d_%d" % (k,  l)
                tarpath =  'ChristyLinks/%s'%filename
                try:
                    filedata = self.Christy_links.extractfile(tarpath).read()
                    c_triangulation = read_triangulation_from_string(filedata)
                except: 
                    raise IOError('The link complement %s was not '
                                  'found.'%filename)
                result =  Manifold('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result              

#----------------------------------------------------------------
#
#  Morwen's link table (data taken from Snap)
#
#----------------------------------------------------------------

def left_pad_string(str,  length, pad=' '):
    return pad*(length - len(str)) + str

class MorwenLinks(Census):
    """
    Morwen Thistlethwaite's table of links with at most 14 crossings
    (about 180k links).  For instance, to look at first few
    2-component links do:

    >>> C = MorwenLinks(2)
    >>> for M in C[:3]:
    ...     print( M, M.volume() )
    ... 
    DT[ebbccdaeb](0,0)(0,0) 3.66386237671
    DT[fbbdceafbd](0,0)(0,0) 5.33348956690
    DT[fbccdefacb](0,0)(0,0) 4.05976642564

    To look at those with 3 components and 11 crossings do:

    >>> C = MorwenLinks(3, 11)
    >>> len(C)   # How many such links are there?
    329
    """
    def __init__(self, num_components, num_crossings = None):
        if num_components < 0:
            return
        
        crossings = range(4, 15) if num_crossings == None else [num_crossings]
        files = []
        for c in crossings:
            n = left_pad_string('%d' % c, 2, '0')
            files.append('hyperbolic_data_%sa.gz' % n)
            if c > 5:
                files.append('hyperbolic_data_%sn.gz' % n)

        self.files = files

        self.DT_codes = []
        second_char = string.ascii_lowercase[num_components-1]
        for file in files:
            data = gzip.open(morwen_link_directory + os.sep + file).readlines()
            for line in data:
                strline = line.decode()
                if strline[1] == second_char:
                    self.DT_codes.append(strline.strip())

        self.length =  len(self.DT_codes)
        Census.__init__(self, indices=(0, len(self.DT_codes), 1))

    def __getitem__(self, n):
        if isinstance(n, slice):
            SC = self.__class__(-1)
            SC.DT_codes = self.DT_codes[n]
            SC.length = len(SC.DT_codes)
            Census.__init__(SC, indices=(0, len(SC.DT_codes), 1))
            return SC
        else:
            return Manifold( 'DT[%s]' % self.DT_codes[n])


# Creating fibered manifolds from braids

cdef c_Triangulation* get_fibered_manifold_associated_to_braid(num_strands,
                                                               braid_word):
    if num_strands < 2:
        raise ValueError('Braids must have at least 2 strands.')
    allowed_letters = ([x for x in range(1,num_strands)] +
                       [x for x in range(-num_strands+1, 0)])
    if False in [b in allowed_letters for b in braid_word]:
        raise ValueError('The braid word is invalid.')

    cdef int* word
    cdef c_Triangulation* c_triangulation

    n = len(braid_word)
    word = <int*>malloc(n*sizeof(int))
    for  i, w in enumerate(braid_word):
        word[i] = w
    c_triangulation = fibered_manifold_associated_to_braid(num_strands, n, word)
    free(word)
    name = to_byte_str('braid' + repr(braid_word))
    set_triangulation_name(c_triangulation,name)
    return c_triangulation

# Link Editor support

strandtype = {'X': KLPStrandX,     'Y': KLPStrandY}
signtype =   {'R': KLPHalfTwistCL, 'L': KLPHalfTwistCCL}

cdef c_Triangulation* get_triangulation_from_PythonKLP(pythonklp) except *:
    cdef KLPProjection P
    cdef c_Triangulation  *c_triangulation

    P.num_crossings, P.num_free_loops, P.num_components, Pycrossings = pythonklp
    P.crossings = <KLPCrossing *>malloc(P.num_crossings * sizeof(KLPCrossing))
    if P.crossings == NULL:
        raise RuntimeError('Could not allocate a crossing table.')

    for i from 0 <= i < P.num_crossings:
        cr_dict = Pycrossings[i]
        neighbor = cr_dict['Xbackward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandX][<int>KLPBackward] = &P.crossings[neighbor] 
        neighbor =  cr_dict['Xforward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandX][<int>KLPForward] = &P.crossings[neighbor]
        neighbor = cr_dict['Ybackward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandY][<int>KLPBackward] = &P.crossings[neighbor]
        neighbor = cr_dict['Yforward_neighbor']
        P.crossings[i].neighbor[<int>KLPStrandY][<int>KLPForward] = &P.crossings[neighbor]

        strand = cr_dict['Xbackward_strand']
        P.crossings[i].strand[<int>KLPStrandX][<int>KLPBackward] = strandtype[strand]
        strand = cr_dict['Xforward_strand']
        P.crossings[i].strand[<int>KLPStrandX][<int>KLPForward] = strandtype[strand]
        strand = cr_dict['Ybackward_strand']
        P.crossings[i].strand[<int>KLPStrandY][<int>KLPBackward] = strandtype[strand]
        strand = cr_dict['Yforward_strand']
        P.crossings[i].strand[<int>KLPStrandY][<int>KLPForward] = strandtype[strand]

        sign = cr_dict['sign']
        P.crossings[i].handedness = signtype[sign[0]]
        P.crossings[i].component[<int>KLPStrandX] = cr_dict['Xcomponent']
        P.crossings[i].component[<int>KLPStrandY] = cr_dict['Ycomponent']

    c_triangulation = triangulate_link_complement(&P);
    free(P.crossings)

    if c_triangulation == NULL:
        raise RuntimeError('Could not create the triangulation.')

    # The triangulation must have a name or SnapPea will segfault when
    #   trying to copy the triangulation.
    tri_name = to_byte_str('unnamed link')
    set_triangulation_name(c_triangulation, tri_name)
    return c_triangulation

def triangulate_link_complement_from_data(data):
    cdef c_Triangulation* c_triangulation
    M = Manifold('empty')
    c_triangulation = get_triangulation_from_PythonKLP(data)
    M.set_c_triangulation(c_triangulation)
    return M

