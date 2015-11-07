# Python modules
import os, sys, operator, types, re, gzip, struct, tempfile
import tarfile, atexit, math, string, time
python_major_version = sys.version_info[0]

# Exceptions from the SnapPea kernel
class SnapPeaFatalError(Exception):
    """
    This exception is raised by SnapPy when the SnapPea kernel
    encounters a fatal error.
    """
# Sage interaction
try:
    import sage.all
    import sage.structure.sage_object
    from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.groups.free_group import FreeGroup
    from sage.interfaces.gap import gap
    from sage.interfaces.gap import is_GapElement
    from sage.interfaces.magma import magma
    from sage.interfaces.magma import is_MagmaElement
    from sage.matrix.constructor import matrix
    # for testing:
    from sage.matrix.constructor import matrix as sage_matrix
    try:
        from sage.libs.pari.gen import pari as pari
    except ImportError:
        from sage.libs.pari.pari_instance import pari as pari
    from sage.libs.pari.gen import gen as gen
    _within_sage = True
except ImportError:
    from cypari.gen import pari as pari
    from cypari.gen import gen as gen
    _within_sage = False

## SnapPy components
import spherogram
from .manifolds import __path__ as manifold_paths
from . import database
from . import twister
from . import snap
from . import verify
from . import decorated_isosig
from .ptolemy import manifoldMethods as ptolemyManifoldMethods
try:
    from plink import LinkEditor, LinkManager
except:
    LinkEditor, LinkManager = None, None
try:
    from .polyviewer import PolyhedronViewer
except ImportError:
    PolyhedronViewer = None
try:
    from .horoviewer import HoroballViewer
except ImportError:
    HoroballViewer = None
try:
    from .browser import Browser
except ImportError:
    Browser = None
try:
    from .filedialog import asksaveasfile
except ImportError:
    asksaveasfile = None

# This is part of the UCS2 hack.
cdef public UCS2_hack (char *string, Py_ssize_t length, char *errors) :   
    return string

# A stream for asynchronous messages
class MsgIO:
    def __init__(self):
        self.write = sys.stdout.write
        self.flush = sys.stdout.flush()

msg_stream = MsgIO()

# A very basic matrix class
class SimpleMatrix:
    """
    A very simple matrix class that wraps a list of lists.  It has
    two indices and can print itself.  Nothing more.
    """
    def __init__(self, list_of_lists):
        self.data = list_of_lists
        try:
            self.type = type(self.data[0][0])
            self.shape = (len(list_of_lists), len(list_of_lists[0]))
        except IndexError:
            self.type = type(0)
            self.shape = (0,0)

    def __iter__(self):
        return self.data.__iter__()

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
            return self.data[key]
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

    def list(self):
        return self.entries()

    __add__ = __sub__ = __mul__ = __div__ = __inv = _noalgebra

if not _within_sage:
    matrix = SimpleMatrix

# Uniform string testing for Python 2 and 3.
if python_major_version == 2:
    def to_str(s):
        return s
    def bytearray_to_bytes(x):
        return str(x)
if python_major_version == 3:
    basestring = unicode = str
    def to_str(s):
        return s.decode()
    def bytearray_to_bytes(x):
        return bytes(x)

def to_byte_str(s):
    return s.encode('utf-8') if type(s) != bytes else s

# Types of covering spaces
cover_types = {1:"irregular", 2:"regular", 3:"cyclic"}

# Paths
manifold_path = manifold_paths[0] + os.sep

# Obsolete data
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
cdef public void uFatalError(const_char_ptr function,
                             const_char_ptr file) except *:
    raise SnapPeaFatalError('SnapPea crashed in function %s(), '
                            'defined in %s.c.'%(function, file))

# Global variables used for interrupt processing 
cdef public Boolean gLongComputationInProgress
cdef public Boolean gLongComputationCancelled
cdef public gLongComputationTicker

# If not None, this will be called in gLongComputationContinues.
# This enables a GUI to do updates during long computations.
UI_callback = None

def SnapPea_interrupt():
    """
    The UI can call this to stop SnapPea.  Returns True if SnapPea s busy.
    If SnapPea is busy, the side effect is to set the gLongComputationCancelled
    flag.
    """
    global gLongComputationCancelled
    global gLongComputationInProgress
    if gLongComputationInProgress:
        gLongComputationCancelled = True
    return gLongComputationInProgress

cdef public void uLongComputationBegins(const_char_ptr message,
                                        Boolean is_abortable):
    global gLongComputationCancelled
    global gLongComputationInProgress
    global gLongComputationTicker
    gLongComputationCancelled = False
    gLongComputationInProgress = True
    gLongComputationTicker = time.time()

cdef public c_FuncResult uLongComputationContinues() except *:
    global gLongComputationCancelled
    global gLongComputationInProgress
    global gLongComputationTicker
    cdef now
    if gLongComputationCancelled:
        return func_cancelled
    elif UI_callback is not None:
        now = time.time()
        if now - gLongComputationTicker > 0.2:
            UI_callback()
            gLongComputationTicker = now
    return func_OK

cdef public void uLongComputationEnds() except*:
    global gLongComputationCancelled
    global gLongComputationInProgress
    gLongComputationInProgress = False
    if gLongComputationCancelled:
        gLongComputationCancelled = False
        if UI_callback is not None:
            UI_callback(interrupted=True)

cdef public void uAcknowledge(const_char_ptr message):
#    sys.stderr.write(<char *> message)
#    sys.stderr.write('\n')
    return

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
            M = sage_matrix(M) 
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
                'contains negatively oriented tetrahedra',
                'contains flat tetrahedra', 'contains degenerate tetrahedra',
                'unrecognized solution type', 'no solution found',
                'tetrahedra shapes were inserted']

# SnapPea memory usage
def check_SnapPea_memory():
    verify_my_malloc_usage()

# Ptolemy utility functions
# convert and free an identification of variables structure
cdef convert_and_free_identification_of_variables(
    Identification_of_variables c_vars):
    var_list = []

    if c_vars.variables:
        for i in range(c_vars.num_identifications):
            var_list.append(
                  (c_vars.signs[i],
                   c_vars.powers[i],
                   to_str(c_vars.variables[i][0]),
                   to_str(c_vars.variables[i][1])))
    return var_list

# convert and free an integer matrix from C
cdef convert_and_free_integer_matrix(
    Integer_matrix_with_explanations c_matrix):
    if not c_matrix.entries:
        return []

    python_matrix = [
        [ c_matrix.entries[i][j] for j in range(c_matrix.num_cols)]
                                 for i in range(c_matrix.num_rows)]

    explain_row = []

    for i in range(c_matrix.num_rows):
        if c_matrix.explain_row[i]:
            explain_row.append(to_str(c_matrix.explain_row[i]))
        else:
            explain_row.append(None)

    explain_column = []

    for i in range(c_matrix.num_cols):
        if c_matrix.explain_column[i]:
            explain_column.append(to_str(c_matrix.explain_column[i]))
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

# Conversions between various numerical types.
IF HIGH_PRECISION:
    from snappy.number import Number as NumberHP
    class Number(NumberHP):
        _default_precision=212
ELSE:
    from snappy.number import Number as NumberLP
    class Number(NumberLP):
        _default_precision=53

cdef Real2gen_direct(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.  This constructs the gen
    directly, but requires the non-sage cypari method
    pari._real_coerced_to_bits_prec.

    A high precision real with 212 bits of precision is converted to
    a gen with 256 bits of precision since pari numbers have precision
    divisible by 32.

    """
    IF HIGH_PRECISION: # Real = qd_real
        cdef double* qd = <double*>&R
        cdef int i
        # The value of a qd_real is the sum of the values of its four doubles.
        cdef result = pari._real_coerced_to_bits_prec(qd[0], 256)
        for i in range(1,4):
            result += pari._real_coerced_to_bits_prec(qd[i], 256)
        return result
    ELSE:              # Real = double
        return pari(R)

cdef Real2gen_string(Real R):
    """
    Convert a Real to a pari gen of type t_REAL.
    This constructs the gen from the string representation of the real.
    """
    return pari(real_to_string(R))

# C type for a function of Real returning an object
ctypedef object (*func_real_to_obj)(Real)

# Convert Real to gen in an appropriate manner for this environment
cdef func_real_to_obj Real2gen

if hasattr(pari, '_real_coerced_to_bits_prec'): # Cypari
    Real2gen = Real2gen_direct
else:
    Real2gen = Real2gen_string

cdef Complex2gen(Complex C):
    """
    Convert a Complex to a pari gen.
    """
    cdef real_part = Real2gen(C.real)
    cdef imag_part = Real2gen(C.imag)
    return pari.complex(real_part, imag_part)

cdef RealImag2gen(Real R, Real I):
        return pari.complex(Real2gen(R), Real2gen(I))

cdef Complex2complex(Complex C):
    """
    Convert a Complex to a python complex.
    """
    return complex( float(<double>C.real), float(<double>C.imag) )

cdef Real2float(Real R):
    """
    Convert a Real to a python float.
    """
    return float(<double>R)

cdef Complex complex2Complex(complex z):
    """
    Convert a python complex to a Complex.
    """
    cdef Complex result
    result.real = <Real>z.real
    result.imag = <Real>z.imag
    return result

cdef Real Object2Real(obj):
    cdef char* c_string
    try:
        string = obj.as_string() if isinstance(obj, Number) else str(obj)
        float(string)
    except:
        raise ValueError('Cannot convert %s to a Real.'%type(obj))
    string = to_byte_str(string)
    c_string = string
    return Real_from_string(c_string)

cdef Complex Object2Complex(obj):
    cdef Real real, imag
    cdef Complex result
    if hasattr(obj, 'real') and hasattr(obj, 'imag'):
        try:
            float(obj.real)
            real = Object2Real(obj.real)
        except TypeError:  # Probably Sage type
            real = Object2Real(obj.real())
        try:
            float(obj.imag)
            imag = Object2Real(obj.imag)
        except TypeError:  # Probably Sage type
            imag = Object2Real(obj.imag())
    else:
        real = Object2Real(obj)
        imag = <Real>0.0
    result.real = real
    result.imag = imag
    return result
             
    

cdef double Real2double(Real R):
    cdef double* quad = <double *>&R
    return quad[0]

cdef Complex gen2Complex(g):
    cdef Complex result
    IF HIGH_PRECISION: # Real = qd_real; 212 bits of precision
        cdef py_string
        cdef char* c_string
        cdef Real real_part, imag_part
        old_precision = pari.set_real_precision(64)

        py_string = to_byte_str(str(g.real()).replace(' E','E')) # save a reference
        c_string = py_string
        real_part = <Real>c_string 
        py_string = to_byte_str(str(g.imag()).replace(' E','E')) # save a reference
        c_string = py_string
        imag_part = <Real>c_string
        result.real, result.imag = real_part, imag_part

        pari.set_real_precision(old_precision)
    ELSE: # Real = double
        result.real, result.imag = g.real(), g.imag()
    return result

#IF HIGH_PRECISION:
cdef Real2Number(Real R):
    return Number(Real2gen(R))
cdef Complex2Number(Complex C):
    return Number(Complex2gen(C))
#ELSE:
#    cdef Real2Number(Real R):
#        return Number(R)
#    cdef Complex2Number(Complex C):
#        return Number(C)

cdef B2B(Boolean B):
    return B != 0

# Callback function which uses Pari to compute the dilogarithm.
# Used by Manifold.complex_volume

IF HIGH_PRECISION:
    cdef Complex dilog_callback(Complex z):
        g = Complex2gen(z)
        # Sometimes we will get a value (e.g. 0.5) that pari decides is
        # "exact" and should be converted to standard precision.  By
        # supplying the precision argument we make sure it gets converted
        # to high precision, i.e. 256 bits in pari.
        li2 = g.dilog(precision=256)
        return gen2Complex(li2)

# Infos are immutable containers which hold information about SnapPy
# objects.  The base class for these is Info. Subclasses should
# override __repr__ to appropriately display the information they
# contain.

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

class NormalSurfaceInfo(Info):
    def __repr__(self):
        orient = 'Orientable' if self.orientable else 'Non-orientable'
        sided = 'two-sided' if self.two_sided else 'one-sided'
        return '%s %s with euler = %s' % (orient, sided, self.euler)
    
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
        return self.divisors

    def rank(self):
        """
        The rank of the group.
        """
        return len(self.divisors)

    def betti_number(self):
        """
        The rank of the maximal free abelian subgroup.
        """
        return len([n for n in self.divisors if n == 0])
    def order(self):
        """
        The order of the group.  Returns the string 'infinite' if the
        group is infinite.        
        """
        det = 1
        for c in self.divisors:
            det = det * c
        return 'infinite' if det == 0 else det

cdef class PresentationMatrix:
    """
    A sparse representation of the presentation matrix of an abelian group.
    """
    cdef rows, cols, _row_support, _col_support, _entries, _units, dead_columns
    def __init__(self, rows, cols):
        self.rows = rows
        self.cols = cols
        self._row_support = {}
        self._col_support = {}
        self._entries = {}
        self._units = set()
        self.dead_columns = set()

    def __setitem__(self, ij, value):
        i, j = ij
        # check bounds
        if not 0 <= i < self.rows or not 0 <= j < self.cols:
            raise IndexError
        self._set(ij, value)

    cdef _set(self, ij, value):
        # cdef'ed to allow for faster calling
        i, j = ij
        # keep track of units
        if value == 1 or value == -1:
            self._units.add(ij)
        else:
            try:
                self._units.remove(ij)
            except KeyError:
                pass
        # don't store any zeros
        if value == 0:
            try:
                self._entries.pop(ij)
                self._row_support[i].remove(j)
                self._col_support[j].remove(i)
            except KeyError:
                pass
            return
        # keep track of row and column supports
        try:
            self._row_support[i].add(j)
        except KeyError:
            self._row_support[i] = set([j])
        try:
            self._col_support[j].add(i)
        except KeyError:
            self._col_support[j] = set([i])
        # set the value
        self._entries[ij] = value
    
    def __getitem__(self, ij):
        return self._entries.get(ij, 0)

    def __repr__(self):
        return repr(self.explode())

    def add_rows(self, n):
        self.rows += n

    def row_op(self, i, j, m):
        """
        Subtract m * row_i from row_j
        """
        for k in list(self._row_support[i]):
            self[j,k] -= m*self[i,k]

    def explode(self):
        """
        Return the full matrix, including dead columns, as a list of lists.
        """
        return [
            [self._entries.get((i,j), 0) for j in xrange(self.cols)]
            for i in xrange(self.rows)]

    
    def simplify(self):
        """
        If any entry is a unit, eliminate the corresponding generator.
        Continue until no units remain.  When a generator is removed,
        remember its column index.
        """
        cdef temp, m, i, j, k, l
        while len(self._units) > 0:
            for i,j in self._units: break 
            col_support = [k for k in self._col_support[j] if k != i] + [i]
            row_entries = [(l, self._entries.get((i,l), 0))
                           for l in self._row_support[i]]
            unit = self[i,j]
            for k in col_support:
                m = unit*self[k,j]
                # this is an "inlined" call to self.row_op(k,i,m)
                # (avoids calling python functions in the loop)
                for l, a_il in row_entries:
                    kl = (k,l)
                    temp = self._entries.get(kl, 0)
                    self._set(kl, temp - m*a_il )
            self.dead_columns.add(j)
    
    def simplified_matrix(self):
        """
        Return the simplified presentation as a matrix.
        """
        self.simplify()
        columns = [j for j in xrange(self.cols) if j not in self.dead_columns]
        rows = [i for i in xrange(self.rows) if self._row_support.get(i,None)]
        if len(rows) == 0:
            presentation = [ [0 for j in columns] ]
        else:
            presentation = [ [self._entries.get((i,j), 0) for j in columns]
                             for i in rows ]
        return matrix(presentation)

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
        + Dowker-Thistlethwaite code: e.g. 'DT:[(6,8,2,4)]', 'DT:dadbcda'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'Braid:[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a0, B1])' or 'Splitting(S_{1,0}, [b1, A0], [a0,B1])'

      See the help for the 'twister' module for more.  

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A Regina-style isomorphism signature, such as 'dLQbcccdxwb'.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """
    cdef c_Triangulation* c_triangulation
    cdef _DTcode
    cdef _cover_info
    cdef readonly _cache
    cdef readonly LE
    cdef _link_file_full_path
    cdef hyperbolic_structure_initialized

    def __cinit__(self, spec=None, remove_finite_vertices=True):
        if UI_callback is not None:
            uLongComputationBegins('Constructing a manifold', 1)
            UI_callback()
            uLongComputationEnds()
        # Answers to potentially hard computations are cached
        self._cache = {}
        self._DTcode = None
        self._cover_info = None
        self.LE = None
        self.hyperbolic_structure_initialized = False
        self._link_file_full_path = None
        if spec is not None and spec != 'empty':
            if not isinstance(spec, basestring):
                raise TypeError(triangulation_help%
                                self.__class__.__name__)
            self.get_triangulation(spec, remove_finite_vertices)
            if self.c_triangulation == NULL:
                raise RuntimeError, 'An empty triangulation was generated.'
        elif spec is None:
            self.get_from_new_plink()

        if self.c_triangulation != NULL and not self.hyperbolic_structure_initialized:    
            remove_hyperbolic_structures(self.c_triangulation)

    cdef get_from_new_plink(self, file_name=None):
        if LinkEditor is None:
            raise RuntimeError, 'PLink was not imported.'
        self.LE = LinkEditor(no_arcs=True,
                             callback=_plink_callback,
                             cb_menu='Send to SnapPy',
                             file_name=file_name,
                             manifold=self)
        print('Starting the link editor.\n'
              'Select Tools->Send to SnapPy to load the '
              'link complement.')

    cdef get_triangulation(self, spec, remove_finite_vertices=True):
        cdef Triangulation T
        
        # Step -1 Check for an entire-triangulation-file-in-a-string
        if spec.startswith('% Triangulation'):
            return self._from_string(spec)

        # Get fillings, if any
        m = split_filling_info.match(spec)
        name = m.group(1)
        fillings = eval( '[' + m.group(2).replace(')(', '),(')+ ']', {})

        # Step 1. Cusped census manifolds
        if is_census_manifold.match(name):
            try:
                database.OrientableCuspedCensus._one_manifold(name, self)
            except KeyError:
                database.NonorientableCuspedCensus._one_manifold(name, self)
            
        # Step 2. The easy databases
        databases = [ (is_HT_link, database.HTLinkExteriors),
                      (is_census_knot, database.CensusKnots),
                      (is_knot_complement, database.LinkExteriors)]

        for regex, db in databases:
            if regex.match(name):
                db._one_manifold(name, self)
                break

        # Step 3. Rolfsen links
        for regex in rolfsen_link_regexs:
            m = regex.match(name)
            if m:
                if int(m.group('components')) > 1:
                    rolfsen_name = '%d^%d_%d' % (int(m.group('crossings')),
                        int(m.group('components')), int(m.group('index')))
                else:
                    rolfsen_name = '%d_%d' % (int(m.group('crossings')),
                                  int(m.group('index')))
                database.LinkExteriors._one_manifold(rolfsen_name, self)

        # Step 4. Hoste-Thistlethwaite knots
        m = is_HT_knot.match(name)
        if m:
            self.get_HT_knot(int(m.group('crossings')), m.group('alternation'),
                        int(m.group('index')))
            
        # Step 5. Once-punctured torus bundles
        m = is_torus_bundle.match(name)
        if m:
            self.get_punctured_torus_bundle(m)

        # Step 6. (fibered) braid complements
        m = is_braid_complement.match(name)
        if m:
            word = eval(m.group(1), {})
            num_strands = max([abs(x) for x in word]) + 1
            self.set_c_triangulation(get_fibered_manifold_associated_to_braid(num_strands, word))

        # Step 7. Dowker-Thistlethwaite codes
        m = is_int_DT_exterior.match(name)
        if m:
            code = eval(m.group(1), {})
            if isinstance(code, tuple):
                knot = spherogram.DTcodec(*code)
            elif isinstance(code, list) and isinstance(code[0], int):
                knot = spherogram.DTcodec([tuple(code)])
            else:
                knot = spherogram.DTcodec(code)
            klp = knot.KLPProjection()
            self.set_c_triangulation(get_triangulation_from_PythonKLP(klp))
            self.set_name(name)
            self._set_DTcode(knot)
            
            
        m = is_alpha_DT_exterior.match(name)
        if m:
            knot = spherogram.DTcodec(m.group(1))
            klp=knot.KLPProjection()
            self.set_c_triangulation(get_triangulation_from_PythonKLP(klp))
            self.set_name(name)
            self._set_DTcode(knot)
            
        # Step 8.  Bundle or splitting is given in Twister's notation

        shortened_name = name.replace(' ', '')
        mb = is_twister_bundle.match(shortened_name)
        ms = is_twister_splitting.match(shortened_name)
        if mb or ms:
            func = bundle_from_string if mb else splitting_from_string
            T = func(shortened_name)
            copy_triangulation(T.c_triangulation, &self.c_triangulation)
            if remove_finite_vertices:
                self._remove_finite_vertices()

        # Step 9. Regina/Burton isomorphism signatures.
        if self.c_triangulation == NULL:
            if is_isosig.match(name):
                self.set_c_triangulation(
                    triangulation_from_isomorphism_signature(name))
                if remove_finite_vertices:
                    self._remove_finite_vertices()
            if is_decorated_isosig.match(name):
                isosig, decoration = name.split('_')
                self.set_c_triangulation(
                    triangulation_from_isomorphism_signature(isosig))
                if remove_finite_vertices:
                    self._remove_finite_vertices()
                decorated_isosig.set_peripheral_from_decoration(self, decoration)
            
        # Step 10. If all else fails, try to load a manifold from a file.
        if self.c_triangulation == NULL:
            self.get_from_file(name, remove_finite_vertices)
        
        # Set the dehn fillings
        Triangulation.dehn_fill(self, fillings)

    cdef get_HT_knot(self, crossings, alternation, index):
        cdef c_Triangulation* c_triangulation
        DT = [get_HT_knot_DT(crossings, alternation, index)]
        knot = spherogram.DTcodec(DT)
        c_triangulation = get_triangulation_from_PythonKLP(knot.KLPProjection())
        name = to_byte_str('%d'%crossings + alternation + '%d'%index)
        set_triangulation_name(c_triangulation, name)
        self._set_DTcode(knot)
        self.set_c_triangulation(c_triangulation)

    cdef get_punctured_torus_bundle(self, match):
        cdef LRFactorization* gluing
        cdef int LRlength, i
        cdef char* LRstring

        LRpart = to_byte_str(match.group(3).upper())
        LRlength = len(LRpart)
        LRstring = LRpart
        gluing = alloc_LR_factorization(LRlength)
        gluing.is_available = True
        gluing.negative_determinant = 1 if match.group(1) in ['-', 'n'] else 0
        gluing.negative_trace = 0 if match.group(2) == '+' else 1
        for i from 0 <= i < LRlength:
            gluing.LR_factors[i] = LRstring[i]
        self.set_c_triangulation(triangulate_punctured_torus_bundle(gluing))
        free_LR_factorization(gluing)

    def _get_from_link_data(self, data):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        self.c_triangulation = get_triangulation_from_PythonKLP(data)

    cdef get_from_file(self, name, remove_finite_vertices=True):
        try:
            locations = [os.curdir, os.environ['SNAPPEA_MANIFOLD_DIRECTORY']]
        except KeyError:
            locations = [os.curdir]
        found = 0
        for location in locations:
            pathname = os.path.join(location, name)
            if os.path.isfile(pathname):
                file = open(pathname, 'r')
                first_line = file.readline()[:-1]
                file.close()
                if first_line.find('% Link Projection') > -1:
                    LM = LinkManager()
                    LM._from_string(open(pathname, 'r').read())
                    klp = LM.SnapPea_KLPProjection()
                    self._link_file_full_path = os.path.abspath(pathname)
                    self._set_DTcode(spherogram.DTcodec(*LM.DT_code()))
                    self.set_c_triangulation(get_triangulation_from_PythonKLP(klp))
                else:
                    self.set_c_triangulation(read_triangulation(pathname))

        if self.c_triangulation == NULL:
            raise IOError('The manifold file %s was not found.\n%s'%
                          (name, triangulation_help % 'Triangulation or Manifold'))
        else:
            self._remove_finite_vertices()

    def _remove_finite_vertices(self):
        """
        Removing any finite vertices by simplification and drilling.

        >>> isosig = 'kLLLLMQkccfigghjijjlnabnwnpsii'
        >>> T = Triangulation(isosig, remove_finite_vertices=False)
        >>> T.triangulation_isosig() == isosig
        True
        >>> T.num_cusps(),  T._num_fake_cusps()
        (0, 1)
        >>> T._remove_finite_vertices()
        >>> T.num_cusps(),  T._num_fake_cusps()
        (1, 0)
        """
        if self.c_triangulation != NULL:
            remove_finite_vertices(self.c_triangulation)

    def cover_info(self):
        """
        If this is a manifold or triangulation which was constructed as
        a covering space, return a dictionary describing the cover.  Otherwise
        return 0.  The dictionary keys are 'base', 'type' and 'degree'.
        """
        if self._cover_info:
            return dict(self._cover_info)
        
    def _clear_cache(self, key=None, message=''):
        # (Cache debugging) print '_clear_cache: %s'%message
        if not key: 
            self._cache.clear()
        else:
            self._cache.pop(key)

    def plink(self):
        """
        Brings up a link editor window if there is a link known to be associated
        with the manifold.
        """
        if self.LE is not None:
            self.LE.reopen()
        elif self._link_file_full_path is not None:
            self.get_from_new_plink(self._link_file_full_path)
        elif self.DT_code() is not None:
            self.get_from_new_plink()
            L = spherogram.DTcodec(self.DT_code()).link()
            L.view(self.LE)
        else:
            raise ValueError('No associated link known.')

    def link(self):
        if self.LE is not None:
            return spherogram.Link(self.LE.PD_code())
        elif self.DT_code() is not None:
            return spherogram.DTcodec(self.DT_code()).link()
        else:
            raise ValueError('No associated link known.')

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

    def _num_fake_cusps(self):
        """
        Returns the number of "fake" cusps, which is typically the number
        of finite vertices.

        >>> M = Triangulation('m004(1,2)')
        >>> F = M.filled_triangulation()
        >>> F.num_cusps(),  F._num_fake_cusps()
        (0, 1)
        >>> S = Triangulation('bkaagb', remove_finite_vertices=False)
        >>> S.num_cusps(),  S._num_fake_cusps()
        (0, 2)
        """
        return get_num_fake_cusps(self.c_triangulation)

    def orientation_cover(self):
        """
        For a non-orientable Triangulation, returns the 2-fold cover which
        is orientable.

        >>> X = Triangulation('x123')
        >>> Y = X.orientation_cover()
        >>> (X.is_orientable(), Y.is_orientable())
        (False, True)
        >>> Y
        x123~(0,0)(0,0)
        >>> Y.cover_info()['type']
        'cyclic'
        """ 
        if self.is_orientable():
            raise ValueError, 'The Triangulation is already orientable.'

        cdef c_Triangulation* cover_c_triangulation = NULL
        cdef Triangulation new_tri

        cover_c_triangulation = double_cover(self.c_triangulation)
        new_tri = self.__class__('empty')
        new_tri.set_c_triangulation(cover_c_triangulation)
        new_tri.set_name(self.name() + '~')
        new_tri._cover_info = {'base'   : self.name(),
                               'type'   : 'cyclic',
                               'degree' : 2}
        return new_tri

    def is_orientable(self):
        """
        Return whether the underlying 3-manifold is orientable.

        >>> M = Triangulation('x124')
        >>> M.is_orientable()
        False
        """
        orientability = Orientability[get_orientability(self.c_triangulation)]
        if orientability == 'orientable': return True
        elif orientability == 'nonorientable': return False
        else: return None

    def copy(self):
        """
        Returns a copy of the triangulation.

        >>> M = Triangulation('m125')
        >>> N = M.copy()
        """
        cdef c_Triangulation* copy_c_triangulation = NULL
        cdef Triangulation new_tri
        
        if self.c_triangulation is NULL:
            return self.__class__('empty')
        copy_triangulation(self.c_triangulation, &copy_c_triangulation)
        new_tri = self.__class__('empty')
        new_tri.set_c_triangulation(copy_c_triangulation)
        return new_tri

    def randomize(self):
        """
        Perform random Pachner moves on the underlying triangulation.

        >>> M = Triangulation('Braid:[1,2,-3,-3,1,2]')
        >>> M.randomize()
        """
        if self.c_triangulation is NULL: return
        randomize_triangulation(self.c_triangulation)
        self._clear_cache(message='randomize')

    def simplify(self):
        """
        Try to simplify the triangulation by doing Pachner moves.

        >>> M = Triangulation('12n123')
        >>> M.simplify()
        """
        if self.c_triangulation is NULL: return
        basic_simplification(self.c_triangulation)
        self._clear_cache(message='simplify')

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
        2.02988321
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
        >>> M.save('fig-eight.tri')     #doctest: +SKIP
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
        WARNING: Users should not use this function directly.  To
        create a Triangulation or Manifold or ManifoldHP from a string
        containing the contents of a triangulation file, simply do:

        >>> M = Manifold('7_4')
        >>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N == M
        True
        
        Fill an empty Triangulation from a string generated by
        _to_string.
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

    def _from_bytes(self, bytestring, *args):
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

    def __reduce__(self):
        """
        Used to pickle Manifolds.
        """
        return (self.__class__, (self._to_string(),))

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

    def isomorphisms_to(self, Triangulation other not None):
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
        return to_str(get_triangulation_name(self.c_triangulation))

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
    
    def DT_code(self, alpha=False, flips=False):
        """
        Return the Dowker-Thistlethwaite code of this link complement,
        if it is a link complement. The DT code is intended to be an
        immutable attribute, for use with knot and link exteriors
        only, which is set only when the manifold was created.
        
        Here is the Whitehead link:
        
        >>> M = Manifold('L5a1')
        >>> M.DT_code()
        [(6, 8), (2, 10, 4)]
        >>> M.DT_code(alpha=True)
        'ebbccdaeb'
        >>> M.DT_code(alpha=True, flips=True)
        'ebbccdaeb.01110'
        >>> M.DT_code(flips=True)
        ([(6, 8), (2, 10, 4)], [0, 1, 1, 1, 0])
        """
        codec = self._DTcode
        if self._DTcode is None:
            return None
        if alpha:
            return codec.encode(header=False, flips=flips)
        else:
            if flips:
                return codec.code, [int(x) for x in codec.flips]
            else:
                return codec.code

    def _set_DTcode(self, code):
        assert isinstance(code, spherogram.DTcodec)
        self._DTcode = code

    def num_tetrahedra(self):
        """
        Return the number of tetrahedra in the triangulation.

        >>> M = Triangulation('m004')
        >>> M.num_tetrahedra()
        2
        """
        if self.c_triangulation is NULL: return 0
        return get_num_tetrahedra(self.c_triangulation)
    
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

        num_cusps = self.num_cusps()

        if which_cusp != None:
            try:
                which_cusp = range(num_cusps)[which_cusp]
            except IndexError:
                raise IndexError('The specified cusp (%s) does not '
                                 'exist.'%which_cusp)

            meridian, longitude = filling_data
            complete = ( meridian == 0 and longitude == 0)
            set_cusp_info(self.c_triangulation,
                          which_cusp, complete,
                          Object2Real(meridian),
                          Object2Real(longitude))
            self._clear_cache(message='dehn_fill')
        else:
            if num_cusps > 1 and len(filling_data) == 2:
                if ( not hasattr(filling_data, '__getitem__')
                     or not hasattr(filling_data[0], '__getitem__') ):
                    raise IndexError('If there is more than one cusp '
                                     'you must specify which one you\n'
                                     'are filling, e.g. M.dehn_fill((2,3),1)')
            if num_cusps == 1 and len(filling_data) == 2:
                Triangulation.dehn_fill(self, filling_data, 0)
                return 
            if len(filling_data) > num_cusps:
                raise IndexError('You provided filling data for too '
                                 'many cusps.  There are only %s.'%
                                 num_cusps)
            for i, fill in enumerate(filling_data):
                Triangulation.dehn_fill(self, fill, i)

    # When doctesting, the M,L coefficients acquire an accuracy of 8.
    # So we have to include the zeros in the doctest string.
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
        >>> sorted(c.keys())
        ['filling', 'index', 'is_complete', 'topology']

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
        cdef Real m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_accuracy, current_shape_accuracy,
        cdef Complex initial_modulus, current_modulus

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
           'filling' : (Real2float(m), Real2float(l))
           }
                
        return CuspInfo(**info)

    def reverse_orientation(self):
        """
        Reverses the orientation of the Triangulation, presuming that
        it is orientable.

        >>> M = Manifold('m015')
        >>> cs = M.chern_simons()
        >>> M.reverse_orientation()
        >>> round(abs(cs + M.chern_simons()), 15)
        0.0
        """
        if not self.is_orientable():
            raise ValueError("The Manifold is not orientable, so its "
                             "orientation can't be reversed.")
        reorient(self.c_triangulation)
        self._clear_cache(message='reverse_orientation')
            
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
            fill_cusp_spec[i] = True if i in cusps_to_fill else False
        fill_all = True if not False in [i in cusps_to_fill
                                      for i in range(n)] else False
        c_filled_tri = fill_cusps(self.c_triangulation,
                                  fill_cusp_spec, '', fill_all)
        free(fill_cusp_spec)
        filled_tri = self.__class__('empty')
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
        M.gluing_equations_pgl(N = 2, equation_type='all')

        Returns a NeumannZagierTypeEquations object that contains a matrix
        encoding the gluing equations for boundary-parabolic PGL(N,C)
        representations together with explanations of the meaning
        of the rows and the columns of the matrix.

        This method generalizes gluing_equations() to PGL(N,C)-representations
        as described in
        Stavros Garoufalidis, Matthias Goerner, Christian K. Zickert: 
        "Gluing Equations for PGL(n,C)-Representations of 3-Manifolds"
        (http://arxiv.org/abs/1207.6711).

        The result of the traditional gluing_equations() can be obtained from
        the general method by:

        >>> M = Triangulation('m004')
        >>> M.gluing_equations_pgl().matrix
        matrix([[ 2,  1,  0,  1,  0,  2],
                [ 0,  1,  2,  1,  2,  0],
                [ 1,  0,  0,  0, -1,  0],
                [ 0,  0,  0,  0, -2,  2]])

        But besides the matrix, the method also returns explanations of
        the columns and rows:

        >>> M = Triangulation("m004")
        >>> M.gluing_equations_pgl()
        NeumannZagierTypeEquations(
          matrix([[ 2,  1,  0,  1,  0,  2],
                  [ 0,  1,  2,  1,  2,  0],
                  [ 1,  0,  0,  0, -1,  0],
                  [ 0,  0,  0,  0, -2,  2]]),
          explain_columns = ['z_0000_0', 'zp_0000_0', 'zpp_0000_0', 'z_0000_1', 'zp_0000_1', 'zpp_0000_1'],
          explain_rows = ['edge_0_0', 'edge_0_1', 'meridian_0_0', 'longitude_0_0'])

        The first row of the matrix means that the edge equation for
        edge 0 is
        
           z_0000_0 ^ 2 * zp_0000_0 * z_0000_1 * zpp_0000_1 ^ 2 = 1.

        Similarly, the next row encodes the edge equation for the other edge
        and the next two rows encode peripheral equations.

        Following the SnapPy convention, a z denotes the cross ratio z at the
        edge (0,1), a zp the cross ratio z' at the edge (0,2) and a zpp the cross
        ratio z" at the edge (1,2). The entire symbol z_xxxx_y then
        denotes the cross ratio belonging to the subsimplex at integral
        point xxxx (always 0000 for N = 2) of the simplex y. Note: the
        SnapPy convention is different from the paper
        mentioned above, e.g., compare 
        kernel_code/edge_classes.c with Figure 3. We follow the SnapPy
        convention here so that all computations done in SnapPy are
        consistent.

        The explanations of the rows and columns can be obtained explicitly by:

        >>> M.gluing_equations_pgl(N = 3, equation_type = 'peripheral').explain_rows
        ['meridian_0_0', 'meridian_1_0', 'longitude_0_0', 'longitude_1_0']
        >>> M.gluing_equations_pgl(N = 2).explain_columns
        ['z_0000_0', 'zp_0000_0', 'zpp_0000_0', 'z_0000_1', 'zp_0000_1', 'zpp_0000_1']

        A subset of all gluing equations can be obtained by setting the
        equation_type:

        * all gluing equations: 'all'
        * non-peripheral equations: 'non_peripheral'

          * edge gluing equations: 'edge'
          * face gluing equations: 'face'
          * internal gluing equations: 'internal'

        * cusp gluing equations: 'peripheral'

          * cusp gluing equations for meridians: 'meridian'
          * cusp gluing equations for longitudes: 'longitude'
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
                
    def _ptolemy_equations_identified_coordinates(self, N,
                                                  obstruction_class = None):

        """
        Ptolemy coordinates that need to be identified for the given
        triangulation when computing pSL(N,C) representations. 
        """
 
        cdef Identification_of_variables c_vars
        cdef int *c_obstruction_class = NULL
        cdef int T

        if N < 2 or N > 15:
            raise ValueError('N has to be 2...15')

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        if obstruction_class:
            T = get_num_tetrahedra(self.c_triangulation)
            if 2 * T != len(obstruction_class):
                raise ValueError('Obstruction class has wrong length')
            c_obstruction_class = <int *>malloc(2*T*sizeof(int))
            for i, c in enumerate(obstruction_class):
                c_obstruction_class[i] = c

        get_ptolemy_equations_identified_coordinates(
            self.c_triangulation, &c_vars, N, c_obstruction_class)

        free(c_obstruction_class)

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

    def _ptolemy_equations_boundary_map_1(self):
        """
        Boundary map C_1 -> C_0 in cellular homology represented as matrix.
        This will compute the homology of the cell complex obtained when
	gluing together the tetrahedra and not of the cusped manifold.

        Also see _ptolemy_equations_boundary_map_3.
        """
        
        cdef Integer_matrix_with_explanations c_matrix

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        get_ptolemy_equations_boundary_map_1(self.c_triangulation, &c_matrix)

        m, explain_rows, explain_cols = convert_and_free_integer_matrix(
            c_matrix)

        return m, explain_rows, explain_cols

    def ptolemy_obstruction_classes(self):

        """
        Returns the obstruction classes needed to compute
        pSL(N,C) = SL(N,C)/{+1,-1} representations for even N, i.e., it
        returns a list with a representative cocycle for each class in
        H^2(M, boundary M; Z/2). The first element in the list is always
        representing the trivial obstruction class.

        For example, 4_1 has two obstruction classes:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_obstruction_classes()
        >>> len(c)
        2
        
        The primary use of these obstruction classes is to construct
        the Ptolemy variety as described in Definition 1.7 of 
        Stavros Garoufalidis, Dylan Thurston, Christian K. Zickert:
        "The Complex Volume of SL(n,C)-Representations of 3-Manifolds"
        (http://arxiv.org/abs/1111.2828).

        For example, to construct the Ptolemy variety for 
        PSL(2,C)-representations of 4_1 that do not lift to boundary-parabolic
        SL(2,C)-representations, use:
    
        >>> p = M.ptolemy_variety(N = 2, obstruction_class = c[1])

        Or the following short-cut:
    
        >>> p = M.ptolemy_variety(2, obstruction_class = 1)

        Note that this obstruction class only makes sense for even N:

        >>> p = M.ptolemy_variety(3, obstruction_class = c[1])
        Traceback (most recent call last):
        ...
        AssertionError: PtolemyObstructionClass only makes sense for even N, try PtolemyGeneralizedObstructionClass
                
        To obtain PGL(N,C)-representations for N > 2, use the generalized
        obstruction class:

        >>> c = M.ptolemy_generalized_obstruction_classes(3)
        >>> p = M.ptolemy_variety(3, obstruction_class = c[1])
    
        The orginal obstruction class encodes a representing cocycle in Z/2 as follows:
        
        >>> c = M.ptolemy_obstruction_classes()
        >>> c[1]
        PtolemyObstructionClass(s_0_0 + 1, s_1_0 - 1, s_2_0 - 1, s_3_0 + 1, s_0_0 - s_0_1, s_1_0 - s_3_1, s_2_0 - s_2_1, s_3_0 - s_1_1)

        This means that the cocycle to represent this obstruction class in Z/2
        takes value 1 in Z/2 on face 0 of tetrahedra 0 (because s_0_0 = -1)
        and value 0 in Z/2 on face 1 of tetrahedra 0 (because s_1_0 = +1).

        Face 3 of tetrahedra 0 and face 1 of tetrahedra 1 are identified,
        hence the cocycle takes the same value on those two faces (s_3_0 = s_1_1).

        """

        return ptolemyManifoldMethods.get_ptolemy_obstruction_classes(self)

    def ptolemy_generalized_obstruction_classes(self, N):

        """
        M.ptolemy_generalized_obstruction_classes(N)

        Returns the obstruction classes needed to compute
        PGL(N,C)-representations for any N, i.e., it returns a list with
        a representative cocycle for each element in
        H^2(M, boundary M; Z/N) / (Z/N)^* where (Z/N)^* are the units in Z/N.
        The first element in the list always corresponds to the trivial
        obstruction class.
        The generalized ptolemy obstruction classes are thus a generalization
        of the ptolemy obstruction classes that allow to find all
        boundary-unipotent
        PGL(N,C)-representations including those that do not lift to
        boundary-unipotent SL(N,C)-representations for N odd or
        SL(N,C)/{+1,-1}-representations for N even.
        
        For example, 4_1 has three obstruction classes up to equivalence:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_generalized_obstruction_classes(4)
        >>> len(c)
        3

        For 4_1, we only get three obstruction classes even though we have
        H^2(M, boundary M; Z/4) = Z/4 because the two obstruction classes 
        1 in Z/4 and -1 in Z/4 are related by a unit and thus give
        isomorphic Ptolemy varieties.

        The primary use of an obstruction class sigma is to construct the
        Ptolemy variety of sigma. This variety computes boundary-unipotent
        PGL(N,C)-representations whose obstruction class to a
        boundary-unipotent lift to SL(N,C) is sigma.

        For example for 4_1, there are 2 obstruction classes for N = 3:

        >>> M = Manifold("4_1")
        >>> c = M.ptolemy_generalized_obstruction_classes(3)
        >>> len(c)
        2

        The Ptolemy variety parametrizing boundary-unipotent
        SL(3,C)-representations of 4_1 is obtained by

        >>> p = M.ptolemy_variety(N = 3, obstruction_class = c[0])

        and the Ptolemy variety parametrizing boundary-unipotent
        PSL(3,C)-representations of 4_1 that do not lift to
        boundary-unipotent SL(3,C)-representations is obtained by

        >>> p = M.ptolemy_variety(N = 3, obstruction_class = c[1])

        The cocycle representing the non-trivial obstruction class looks as
        follows:

        >>> c[1]
        PtolemyGeneralizedObstructionClass([2, 0, 0, 1])

        This means that the cocycle takes the value -1 in Z/3 on the first face
        class and 1 on the fourth face class but zero on every other of the
        four face classes.
        """


        return (
            ptolemyManifoldMethods.get_generalized_ptolemy_obstruction_classes(
                self, N))

    def ptolemy_variety(self, N, obstruction_class = None,
                        simplify = True, eliminate_fixed_ptolemys = False):

        """
        M.ptolemy_variety(N, obstruction_class = None, simplify = True, eliminate_fixed_ptolemys = False)

        Returns a Ptolemy variety as described in

        * Stavros Garoufalidis, Dyland Thurston, Christian K. Zickert: 
          "The Complex Volume of SL(n,C)-Representations of 3-Manifolds"
          (http://arxiv.org/abs/1111.2828)
        * Stavros Garoufalidis, Matthias Goerner, Christian K. Zickert:
          "Gluing Equations for PGL(n,C)-Representations of 3-Manifolds "
          (http://arxiv.org/abs/1207.6711)
        
        The variety can be exported to magma or sage and solved there. The
        solutions can be processed to compute invariants. The method can also
        be used to automatically look up precomputed solutions from the
        database at http://ptolemy.unhyperbolic.org/data .

        Example for m011 and PSL(2,C)-representations:

        >>> M = Manifold("m011")

        Obtain all Ptolemy varieties for PSL(2,C)-representations:

        >>> p = M.ptolemy_variety(2, obstruction_class = 'all')
        
        There are two Ptolemy varieties for the two obstruction classes:

        >>> len(p)
        2

        Retrieve the solutions from the database

        >>> sols = p.retrieve_solutions() #doctest: +SKIP

        Compute the solutions using magma (default in SnapPy)

        >>> sols = p.compute_solutions(engine = 'magma') #doctest: +SKIP
        
        Compute the solutions using singular (default in sage)

        >>> sols = p.compute_solutions(engine = 'sage') #doctest: +SKIP
        
        Note that magma is significantly faster.

        Compute all resulting complex volumes

        >>> cvols = sols.complex_volume_numerical() #doctest: +SKIP
        >>> cvols  #doctest: +SKIP
        [[[-4.29405713186238 E-16 + 0.725471193740844*I,
           -0.942707362776931 + 0.459731436553693*I,
           0.942707362776931 + 0.459731436553693*I]],
         [[3.94159248086745 E-15 + 0.312682687518267*I,
           4.64549527022581 E-15 + 0.680993020093457*I,
           -2.78183391239608 - 0.496837853805869*I,
           2.78183391239608 - 0.496837853805869*I]]]

        Show complex volumes as a non-nested list:

        >>> cvols.flatten(depth=2) #doctest: +SKIP
        [-4.29405713186238 E-16 + 0.725471193740844*I,
         -0.942707362776931 + 0.459731436553693*I,
         0.942707362776931 + 0.459731436553693*I,
         3.94159248086745 E-15 + 0.312682687518267*I,
         4.64549527022581 E-15 + 0.680993020093457*I,
         -2.78183391239608 - 0.496837853805869*I,
         2.78183391239608 - 0.496837853805869*I]
        
        For more examples, go to http://ptolemy.unhyperbolic.org/

        === Optional Arguments ===
        
        obstruction_class --- class from Definiton 1.7 of (1).
        None for trivial class or a value returned from ptolemy_obstruction_classes.
        Short cuts: obstruction_class = 'all' returns a list of Ptolemy varieties
        for each obstruction. For easier iteration, can set obstruction_class to 
        an integer.
        
        simplify --- boolean to indicate whether to simplify the equations which
        significantly reduces the number of variables.
        Simplifying means that several identified Ptolemy coordinates x = y = z = ...
        are eliminated instead of adding relations x - y = 0, y - z = 0, ...
        
        eliminate_fixed_ptolemys --- boolean to indicate whether to eliminate
        the Ptolemy coordinates that are set to 1 for fixing the decoration.
        Even though this simplifies the resulting representation, setting it to
        True can cause magma to run longer when finding a Groebner basis.

        === Examples for 4_1 ===
        
        >>> M = Manifold("4_1")
        
        Get the varieties for all obstruction classes at once (use
        help(varieties[0]) for more information):
        
        >>> varieties = M.ptolemy_variety(2, obstruction_class = "all")
        
        Print the variety as an ideal (sage object) for the non-trivial class:

        >>> varieties[1].ideal    #doctest: +SKIP
        Ideal (-c_0011_0^2 + c_0011_0*c_0101_0 + c_0101_0^2, -c_0011_0^2 - c_0011_0*c_0101_0 + c_0101_0^2, c_0011_0 - 1) of Multivariate Polynomial Ring in c_0011_0, c_0101_0 over Rational Field                                                       

        Print the equations of the variety for the non-trivial class:
        
        >>> for eqn in varieties[1].equations:
        ...     print(eqn)          #doctest: +NORMALIZE_WHITESPACE
             - c_0011_0 * c_0101_0 + c_0011_0^2 + c_0101_0^2
             c_0011_0 * c_0101_0 - c_0011_0^2 - c_0101_0^2
             - 1 + c_0011_0
 
        Generate a magma file to compute Primary Decomposition for N = 3:
        
        >>> p = M.ptolemy_variety(3)
        >>> s = p.to_magma()
        >>> print(s.split("ring and ideal")[1].strip())     #doctest: +ELLIPSIS
        R<c_0012_0, c_0012_1, c_0102_0, c_0111_0, c_0201_0, c_1011_0, c_1011_1, c_1101_0> := PolynomialRing(RationalField(), 8, "grevlex");
        MyIdeal := ideal<R |
                  c_0012_0 * c_1101_0 + c_0102_0 * c_0111_0 - c_0102_0 * c_1011_0,
            ...
        
        === If you have a magma installation ===
        
        Call p.compute_solutions() to automatically call magma on the above output
        and produce exact solutions!!!
        
        >>> try:
        ...     sols = p.compute_solutions()
        ... except:
        ...     # magma failed, use precomputed_solutions
        ...     sols = None

        Check solutions against manifold
        >>> if sols:
        ...     dummy = sols.check_against_manifold()
        
        === If you do not have a magma installation ===
        
        Load a precomputed example from magma which is provided with the package:
        
        >>> from snappy.ptolemy.processMagmaFile import _magma_output_for_4_1__sl3, solutions_from_magma
        >>> print(_magma_output_for_4_1__sl3)      #doctest: +ELLIPSIS
        <BLANKLINE>
        ==TRIANGULATION=BEGINS==
        % Triangulation
        4_1
        ...
        
        Parse the file and produce solutions:
        
        >>> sols = solutions_from_magma(_magma_output_for_4_1__sl3)

        >>> dummy = sols.check_against_manifold()
            
        === Continue here whether you have or do not have magma ===
        
        Pick the first solution of the three different solutions (up to Galois
        conjugates):
        
        >>> len(sols)
        3
        >>> solution = sols[0]
        
        Read the exact value for c_1020_0 (help(solution) for more information
        on how to compute cross ratios, volumes and other invariants):
        
        >>> solution['c_1020_0']
        Mod(-1/2*x - 3/2, x^2 + 3*x + 4)
        
        Example of simplified vs non-simplified variety for N = 4:
        
        >>> simplified = M.ptolemy_variety(4, obstruction_class = 1)
        >>> full = M.ptolemy_variety(4, obstruction_class = 1, simplify = False)
        >>> len(simplified.variables), len(full.variables)
        (21, 63)
        >>> len(simplified.equations), len(full.equations)
        (24, 72)
        """
        
        return ptolemyManifoldMethods.get_ptolemy_variety(
            self, N, obstruction_class,
            simplify = simplify,
            eliminate_fixed_ptolemys = eliminate_fixed_ptolemys)

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
                                        i, int(m), int(l), &num_rows)
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

    cdef big_homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.
        Preliminary simplification is done with arbitrary precision
        integers.  Smith form is then computed with PARI.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z
        """
        if not all_Dehn_coefficients_are_integers(self.c_triangulation):
            raise ValueError('All Dehn filling coefficients must be integers')
        choose_generators(self.c_triangulation, 0, 0)
        num_generators = self.c_triangulation.num_generators
        cdef EdgeClass* edge
        cdef PositionedTet ptet, ptet0
        cdef GeneratorStatus status
        cdef c_Tetrahedron* tet
        cdef VertexIndex vertex
        cdef FaceIndex side
        cdef Orientation orientation
        cdef int row, column, num_edges, m, l
        relation_matrix = PresentationMatrix(0, num_generators)
        edge = self.c_triangulation.edge_list_begin.next
        # find the edge relations
        while edge != &(self.c_triangulation.edge_list_end):
            row = relation_matrix.rows
            relation_matrix.add_rows(1)
            set_left_edge(edge, &ptet0)
            ptet = ptet0
            while True:
                column = ptet.tet.generator_index[ptet.near_face]
                status = ptet.tet.generator_status[ptet.near_face]
                if status == outbound_generator:
                    relation_matrix[row, column] += 1
                elif status == inbound_generator:
                    relation_matrix[row, column] -= 1
                elif status != not_a_generator:
                    raise RuntimeError('Invalid generator status')
                veer_left(&ptet)
                if same_positioned_tet(&ptet, &ptet0):
                    break
            row += 1
            edge = edge.next
        # find the cusp relations
        num_edges = relation_matrix.rows
        relation_matrix.add_rows(self.num_cusps())
        tet = self.c_triangulation.tet_list_begin.next
        while tet != &(self.c_triangulation.tet_list_end): 
            for vertex in range(4): 
                if tet.cusp[vertex].is_complete:
                    continue
                for side in range(4):
                    if side == vertex:
                        continue
                    if tet.generator_status[side] != inbound_generator:
                        continue
                    for orientation in range(2):
                        row = num_edges + tet.cusp[vertex].index
                        column = tet.generator_index[side]
                        m = <int>tet.cusp[vertex].m
                        l = <int>tet.cusp[vertex].l
                        relation_matrix[row, column] += (
                              m*tet.curve[0][orientation][vertex][side]
                            + l*tet.curve[1][orientation][vertex][side]
                        )
            tet = tet.next
        return AbelianGroup(relation_matrix.simplified_matrix())

    cdef csmall_homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.
        Preliminary simplification is done with 64 bit integers.
        Smith form is then computed with PARI.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z
        """
        cdef c_AbelianGroup *H
        cdef RelationMatrix R
        cdef int m, n

        if self.c_triangulation is NULL:
            return AbelianGroup()
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
            raise RuntimeError("The SnapPea kernel couldn't compute "
                             "the homology presentation matrix")
        result = AbelianGroup(relations)
        return result

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
            try:
                result = self.csmall_homology()
            except RuntimeError:
                result = self.big_homology()
        self._cache['homology'] = result
        return result

    def fundamental_group(self,
                          simplify_presentation = True,
                          fillings_may_affect_generators = True,
                          minimize_number_of_generators = True,
                          try_hard_to_shorten_relators = True):
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
        name_mangled = 'fundamental_group-%s-%s-%s-%s' %\
                       (simplify_presentation,
                        fillings_may_affect_generators,
                        minimize_number_of_generators,
                        try_hard_to_shorten_relators)
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = FundamentalGroup(
                self,
                simplify_presentation,
                fillings_may_affect_generators,
                minimize_number_of_generators,
                try_hard_to_shorten_relators)
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
        >>> N0.cover_info()['type']
        'irregular'
        >>> N0.cover_info()['base']
        'm004'
        >>> N0.cover_info()['degree']
        5
        
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
                permutation_rep = [ [x - 1 for x in perm.domain()]
                                   for perm in permutation_rep ]

        G = self.fundamental_group()
        c_representation = self.build_rep_into_Sn(permutation_rep)
        degree = len(permutation_rep[0])
        c_triangulation = construct_cover(self.c_triangulation,
                                          c_representation,
                                          degree)
        cover = self.__class__('empty')
        cover.set_c_triangulation(c_triangulation)
        cover.set_name(self.name() +'~')
        cover._cover_info = {
            'base'   : self.name(),
            'type'   : cover_types[c_representation.covering_type],
            'degree' : degree
            }
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
            T = self.__class__('empty')
            T.set_c_triangulation(cover)
            T._cover_info = info = {
                'base'   : self.name(),
                'type'   : cover_types[rep.covering_type],
                'degree' : degree
                }
            T.set_name(info['base'] + "~" + info['type'][:3] + '~%d' %
                       cover_count)
            covers.append(T)
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
                                             True, True, True, True)
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

        self._clear_cache(message='set_peripheral_curves')
        if peripheral_data == 'fillings':
            if which_cusp != None:
                raise ValueError("You must apply 'fillings' to all "
                                 "of the cusps.")
            install_current_curve_bases(self.c_triangulation)
            return
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

    def has_finite_vertices(self):
        """
        Returns True if and only if the triangulation has finite (non-ideal)
        vertices.

        >>> T = Triangulation("m004")
        >>> T.has_finite_vertices()
        False
        >>> T.dehn_fill((12,13))
        >>> S = T.filled_triangulation()
        >>> S.has_finite_vertices()
        True

        When trying to find a hyperbolic structure, SnapPea will eliminate
        finite vertices:

        >>> M = S.with_hyperbolic_structure()
        >>> M.has_finite_vertices()
        False
        """

        cdef c_Triangulation* copy_c_triangulation = NULL

        # Bail if empty
        if self.c_triangulation is NULL:
            return False

        # Copy so that we don't loose any reindexing of the cusps on the
        # original triangulation
        copy_triangulation(self.c_triangulation, &copy_c_triangulation)
        
        result = B2B(mark_fake_cusps(copy_c_triangulation))

        # Free the temporary copy
        free_triangulation(copy_c_triangulation)

        return result        

    def triangulation_isosig(self, decorated=False):
        """
        Returns a compact text representation of the triangulation, called an
        "isomorphism signature"::

        >>> T = Triangulation('m004')
        >>> T.triangulation_isosig()
        'cPcbbbiht'

        You can use this string to recreate an isomorphic triangulation later
        
        >>> A = Triangulation('y233')
        >>> A.triangulation_isosig()
        'hLMzMkbcdefggghhhqxqhx'
        >>> B = Triangulation('hLMzMkbcdefggghhhqxqhx')
        >>> A == B
        True

        *WARNING:* By default, the returned string does *not* encode
        the peripheral curves, but for orientable manifolds you can request
        a "decorated isosig" which is also a valid specifier for a Triangulation:

        >>> E = Triangulation('K3_1')   # the (-2, 3, 7) exterior
        >>> isosig = E.triangulation_isosig(); isosig
        'dLQacccjsnk'
        >>> F = Triangulation(isosig)
        >>> E.isomorphisms_to(F)[1]
        0 -> 0
        [1 18]
        [0  1]
        Extends to link
        >>> E.triangulation_isosig(True)
        'dLQacccjsnk_BaRsB'
        >>> F.triangulation_isosig(True)
        'dLQacccjsnk_BaaB'
        >>> G = Triangulation('dLQacccjsnk_BaRsB')
        >>> E.isomorphisms_to(G)[0]
        0 -> 0
        [1 0] 
        [0 1] 
        Extends to link

        The code has been copied from `Regina <http://regina.sf.org/>`_ where
        the corresponding method is called ``isoSig``.

        Unlike dehydrations for 3-manifold triangulations, an
        isomorphism signature uniquely determines a triangulation up
        to combinatorial isomorphism.  That is, two triangulations of
        3-dimensional manifolds are combinatorially isomorphic if and
        only if their isomorphism signatures are the same string.  For
        full details, see `Simplification paths in the Pachner graphs
        of closed orientable 3-manifold triangulations, Burton, 2011
        <http://arxiv.org/abs/1110.6080>`.

        For details about how the peripheral decorations work, see
        the SnapPy source code.
        """

        cdef char *c_string
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        name_mangled = 'triangulation_isosig-%s' % decorated
        if not name_mangled in self._cache.keys():
            if not decorated:
                try:
                    c_string = get_isomorphism_signature(self.c_triangulation)
                    self._cache[name_mangled] = to_str(c_string)
                finally:
                    free(c_string)
            else:
                if not self.is_orientable():
                    raise ValueError('Sorry, decorations are only implemented for orientable manifolds')
                self._cache[name_mangled] = decorated_isosig.decorated_isosig(
                    self, _triangulation_class)

        return self._cache[name_mangled]

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
    4.05686022
    >>> M.cusp_info('shape')
    [-4.27893632 + 1.95728680*I]

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
       + Dowker-Thistlethwaite code: e.g. 'DT:[(6,8,2,4)]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'Braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - From mapping class group data using Twister:

      'Bundle(S_{1,1}, [a0, B1])' or 'Splitting(S_{1,0}, [b1, A0], [a0,B1])'

      See the help for the 'twister' module for more.  

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.

    - A string containing the contents of a SnapPea triangulation or link
      projection file.
    """


    def __init__(self, spec=None):
        if self.c_triangulation != NULL:
            self.init_hyperbolic_structure()
            IF HIGH_PRECISION:
                self.c_triangulation.dilog = dilog_callback
            do_Dehn_filling(self.c_triangulation)

    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        """
        A class method for specifying a numerical conversion function.

        SnapPy includes its own number type, snappy.Number, which can
        represent floating point real or complex numbers of varying
        precision.  (In fact, Number is a wrapper for a pari number of
        type 't_INT', 't_FRAC', 't_REAL' or 't_COMPLEX', and the pari
        gen can be extracted as an attribute: x.gen .)  Methods of
        SnapPy objects which return numerical values will first compute
        the value as a Number, and then optionally convert the Number
        to a different numerical type which can be specified by calling
        this class method.

        By default SnapPy returns Numbers when loaded into python, and
        elements of a Sage RealField or ComplexField when loaded into
        Sage.  These will be 64 bit numbers for ordinary Manifolds and
        212 bit numbers for high precision manifolds.

        The func argument should be a function which accepts a number and
        returns a numerical type of your choosing.  Alternatively, the
        strings 'sage' or 'snappy' can be passed as arguments to select
        either of the two default behaviors.

        EXAMPLE::

            sage: M = Manifold('m004')
            sage: parent(M.volume())
            Real Field with 64 bits of precision
            sage: Manifold.use_field_conversion('snappy')
            sage: M = Manifold('m004')
            sage: parent(M.volume())
            SnapPy Numbers with 64 bits precision
            sage: Manifold.use_field_conversion('sage')
            sage: M = Manifold('m004')
            sage: parent(M.volume())
            Real Field with 64 bits of precision
        """
        if func == 'sage':
            cls._number_ = staticmethod(lambda n : n.sage())
        elif func == 'snappy':
            cls._number_ = staticmethod(lambda n : n)
        else:
            cls._number_ = staticmethod(func)

    def init_hyperbolic_structure(self):
        if not self.hyperbolic_structure_initialized:
            find_complete_hyperbolic_structure(self.c_triangulation)
            do_Dehn_filling(self.c_triangulation)
            self.hyperbolic_structure_initialized = True

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

        Note: due to rounding error, it is possible that this is not
        actually the canonical triangulation.          
        """
        cdef c_FuncResult result
        result = proto_canonize(self.c_triangulation)
        if FuncResult[result] != 'func_OK':
            raise RuntimeError('SnapPea failed to find the canonical '
                               'triangulation.')
        self._clear_cache(message='canonize')

    def _canonical_retriangulation(self, opacities = None):
        """
	If this triangulation is a subdivision of the canonical cell
        decomposition, return the canonical retriangulation as Triangulation.
	Warning: Many operations on a SnapPy Triangulation will remove the
        finite vertices or change the triangulation so it is no longer the
        canonical retriangulation.
        By default, the algorithm numerically checks that the tilts are close
        to zero to determine which faces are opaque or transparent.
        But it can also be passed an explicit list of 4 * num_tetrahedra bool's
        (one per face of each tet) that mark the opaque faces.

	For example, m412's canonical cell decomposition consists of a single
        cube. The canonical retriangulation thus has 12 simplices.

        >>> M = Manifold("m412")
        
        Some isometries of M are not visible in this triangulation (not
        even after calling M.canonize())

        >>> len(M.isomorphisms_to(M))
        4
        >>> T = M._canonical_retriangulation()
        >>> T.num_tetrahedra()
        12

        But all isometries of M are visible in its canonical retriangulation.

        >>> len(T.isomorphisms_to(T))
        8

        """

        cdef Boolean *c_opacities
        cdef c_Triangulation *c_retriangulated_triangulation
        cdef char *c_string
        cdef int n = get_num_tetrahedra(self.c_triangulation)
        cdef result

        if self.c_triangulation is NULL: return ""

        copy_triangulation(self.c_triangulation,
                           &c_retriangulated_triangulation)

        if opacities:
            if not len(opacities) == 4 * n:
                raise Exception("Number of opacities does not match.")
            c_opacities = <Boolean *>malloc(n * 4 * sizeof(Boolean))
            for i in range(4 * n):
                c_opacities[i] = 1 if opacities[i] else 0
        else:
            c_opacities = NULL
            result = proto_canonize(c_retriangulated_triangulation)
            if FuncResult[result] != 'func_OK':
                free_triangulation(c_retriangulated_triangulation)
                raise RuntimeError('SnapPea failed to find the canonical '
                                   'triangulation.')

        canonical_retriangulation_with_opacities(
                c_retriangulated_triangulation, c_opacities)

        free(c_opacities)

        new_tri = Triangulation('empty')
        new_tri.set_c_triangulation(c_retriangulated_triangulation)
        new_tri.set_name(self.name() + '_canonical')

        return new_tri

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

    def _from_string(self, string, initialize_structure=True):
        """
        WARNING: Users should not use this function directly.  To create a
        Manifold or ManifoldHP from a string containing the contents
        of a triangulation file, simply do:

        >>> M = Manifold('7_4')
        >>> seed = M._to_string()
        >>> N = Manifold(seed)
        >>> N.volume()
        5.13794120
        """
        Triangulation._from_string(self, string)
        if HIGH_PRECISION:
            self.c_triangulation.dilog = dilog_callback
        if initialize_structure:
            self.init_hyperbolic_structure()


    def _from_bytes(self, bytestring, initialize_structure=True):
        """
        Fill an empty manifold from a byte sequence generated by
        _to_bytes.
        >>> M = Manifold('m125')
        >>> N = Manifold('empty')
        >>> N._from_bytes(M._to_bytes())
        >>> N.is_isometric_to(M)
        True
        """
        Triangulation._from_bytes(self, bytestring)
        if HIGH_PRECISION:
            self.c_triangulation.dilog = dilog_callback
        if initialize_structure:
            self.init_hyperbolic_structure()

    def copy(self):
        """
        Returns a copy of the manifold

        >>> M = Manifold('m015')
        >>> N = M.copy()
        """
        empty = self.__class__('empty')
        if self.c_triangulation is NULL:
            return empty
        return Manifold_from_Triangulation(self, manifold_class=self.__class__)

    def cusp_neighborhood(self):
        """
        Returns information about the cusp neighborhoods of the
        manifold, in the form of data about the corresponding horoball
        diagrams in hyperbolic 3-space.
        
        >>> M = Manifold('s000')
        >>> CN = M.cusp_neighborhood()
        >>> CN.volume()
        0.32475953
        >>> len(CN.horoballs(0.01))
        178
        >>> CN.view()  # Opens picture of the horoballs  #doctest: +CYOPENGL
        """
        return CuspNeighborhood(self)

    def dirichlet_domain(self,
                         vertex_epsilon=default_vertex_epsilon,
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
        >>> D.view()   #Shows 3d-graphical view.  #doctest: +CYOPENGL
        
        Other options can be provided to customize the computation;
        the default choices are shown below:

        >>> M.dirichlet_domain(vertex_epsilon=10.0**-8,  displacement = [0.0, 0.0, 0.0],
        ... centroid_at_origin=True, maximize_injectivity_radius=True)
        32 finite vertices, 2 ideal vertices; 54 edges; 22 faces
        """
        name_mangled = 'dirichlet_domain-%s-%s-%s' % (
            displacement,
            centroid_at_origin,
            maximize_injectivity_radius)
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = DirichletDomain(
                self,
                vertex_epsilon,
                displacement,
                centroid_at_origin,
                maximize_injectivity_radius)
        return self._cache[name_mangled]

    def browse(self):
        """
        >>> M = Manifold('m125')
        >>> M.browse() # Opens browser window  #doctest: +CYOPENGL
        """
        if Browser is None:
            raise RuntimeError("Browser not imported, Tk or CyOpenGL is probably missing.")
        Browser(self)
        
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
        
        Filling cusps 0 and 2 :

        >>> M = Manifold('v3227(1,2)(3,4)(5,6)')
        >>> M.filled_triangulation([0,2])
        v3227_filled(3,4)
        """
        filled = _triangulation_class.filled_triangulation(self, cusps_to_fill)
        if filled.num_cusps() == 0:
            return Triangulation_from_Manifold(filled)
        return filled
       
    def fundamental_group(self,
                   simplify_presentation = True,
                   fillings_may_affect_generators = True,
                   minimize_number_of_generators = True,
                   try_hard_to_shorten_relators = True):
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
        matrix([[ 2.50000000 - 2.59807621*I, -6.06217783 - 0.50000000*I],
                [ 0.86602540 - 2.50000000*I, -4.00000000 + 1.73205081*I]])

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
        name_mangled = 'fundamental_group-%s-%s-%s-%s' %\
                       (simplify_presentation,
                        fillings_may_affect_generators,
                        minimize_number_of_generators,
                        try_hard_to_shorten_relators)
        if not name_mangled in self._cache.keys():
            result = HolonomyGroup(
               self,
               simplify_presentation,
               fillings_may_affect_generators,
               minimize_number_of_generators,
               try_hard_to_shorten_relators)
            result.use_field_conversion(self._number_)
            self._cache[name_mangled] = result
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

            symmetric_triangulation = self.__class__('empty')
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

          sage: M = Manifold('m004')

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
          8.00000000

        Or maybe we want larger cover coming from the kernel of this::

          sage: N3 = M.cover(f.Kernel())
          sage: N3.volume()/M.volume()
          168.00000000

        Check the homology against what Gap computes directly::
        
          sage: N3.homology().betti_number()
          32
          sage: len([ x for x in f.Kernel().AbelianInvariants().sage() if x == 0])
          32

        We can do the same for Magma::

          sage: G = magma(M.fundamental_group())             #doctest: +SKIP
          sage: Q, f = G.pQuotient(5, 1, nvals = 2)          #doctest: +SKIP
          sage: M.cover(f.Kernel()).volume()                 #doctest: +SKIP
          10.14941606
          sage: h = G.SimpleQuotients(1, 11, 2, 10000)[1,1]  #doctest: +SKIP
          sage: N4 = M.cover(h)                              #doctest: +SKIP
          sage: N2 == N4                                     #doctest: +SKIP
          True
        """
        cover = Triangulation.cover(self, permutation_rep)
        return Manifold_from_Triangulation(cover, recompute=False,
                                           manifold_class=self.__class__)

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
        return [Manifold_from_Triangulation(cover,
                                            recompute=False,
                                            manifold_class=self.__class__)
                for cover in covers]
    
    cdef _real_volume(self):
        """
        Return the (real) volume as a Number.
        """
        cdef int acc
        cdef solution_type
        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('Solution type is: %s'%solution_type)
        IF HIGH_PRECISION == True:
            # must provide a start value to get the correct precision
            result = sum(
                [z.volume() for z in self._get_tetrahedra_shapes('filled')],
                Number(0))
        ELSE:
            result = Number(volume(self.c_triangulation, &acc))
            result.accuracy = acc
        return result

    cdef _cusped_complex_volume(self, Complex *volume, int *accuracy):
        """
        Computes the complex volume of the manifold, computed using
        Goerner's implementation of Zickert's algorithm.  This only
        works for manifolds with at least one cusp.  A ValueError
        is raised if all cusps are filled.

        >>> M = Manifold('5_2')
        >>> M.cusped_complex_volume()
        2.828122088 - 3.024128377*I

        The return value has an extra attribute, accuracy, which is
        the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.cusped_complex_volume().accuracy in (11, 63) # Low, High
        True
        """
        cdef const_char_ptr err_msg=NULL
        cdef c_Triangulation* copy_c_triangulation
        cdef c_Triangulation
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        volume[0] = complex_volume(self.c_triangulation,
                                   &err_msg,
                                   accuracy)
        if err_msg is NULL:
            return
        raise ValueError(err_msg)

    def complex_volume(self):
        """
        Returns the complex volume, i.e.
            volume + i 2 pi^2 (chern simons)

        >>> M = Manifold('5_2')
        >>> M.complex_volume()
        2.82812209 - 3.02412838*I
        >>> M = Manifold("3_1")
	>>> M.complex_volume()
        0 - 1.64493407*I
        """
        cdef Complex volume
        cdef int accuracy
        if True in self.cusp_info('is_complete'):
            self._cusped_complex_volume(&volume, &accuracy)
            result = Complex2Number(volume)
            result.accuracy = accuracy
        else:
            result = self._real_volume() + self._chern_simons()*Number('I')
        return self._number_(result)

    def volume(self, accuracy=False):
        """
        Returns the volume of the current solution to the hyperbolic
        gluing equations; if the solution is sufficiently non-degenerate,
        this is the sum of the volumes of the hyperbolic pieces in
        the geometric decomposition of the manifold. 

        >>> M = Manifold('m004')
        >>> M.volume()
        2.02988321
        >>> M.solution_type()
        'all tetrahedra positively oriented'

        The return value has an extra attribute, accuracy, which is the
        number of digits of accuracy as *estimated* by SnapPea.  When
        printing the volume, the result is rounded to 1 more than this
        number of digits.

        >>> M.volume().accuracy in (11, 63) # Low precision, High precision
        True
        """
        vol = self._real_volume()
        if accuracy:
            return (self._number_(vol), vol.accuracy)
        else:
            return self._number_(vol)

    cdef _old_chern_simons(self):
        """
        Return the Chern-Simons invariant as a Number.  The value
        is computed using SnapPea's original algorithm, which is
        based on Meyerhoff-Hodgson-Neumann.
        """
        cdef Boolean is_known, requires_initialization
        cdef Real CS
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
        cs = Real2Number(CS)
        cs.accuracy = accuracy
        return cs

    cdef _chern_simons(self):
        """
        Return the Chern-Simons invariant as a Number.  The value
        is computed with Zickert's algorithm, for cusped manifolds,
        or by _old_chern_simons if all cusps are filled.
        """
        cdef Complex volume
        cdef Real cs_value
        cdef int accuracy
        if self.c_triangulation is NULL:
            return 0
        solution_type = self.solution_type()
        if solution_type in ('not attempted', 'no solution found'):
            raise ValueError('The solution type is: %s'%solution_type)
        if not True in self.cusp_info('is_complete'):
           result = self._old_chern_simons()
        else:
            self._cusped_complex_volume(&volume, &accuracy)
            cs_value = volume.imag / PI_SQUARED_BY_2
            result = Real2Number(cs_value)
            result.accuracy = accuracy - 1 if accuracy else None
            set_CS_value(self.c_triangulation, cs_value)
        return result

    def chern_simons(self):
        """
        Returns the Chern-Simons invariant of the manifold, if it is known.

        >>> M = Manifold('m015')
        >>> M.chern_simons()
        -0.15320413

        The return value has an extra attribute, accuracy, which
        is the number of digits of accuracy as *estimated* by SnapPea.

        >>> M.chern_simons().accuracy in (8, 9, 57) # Low and High precision
        True

        By default, when the manifold has at least one cusp, Zickert's
        algorithm is used; when the manifold is closed we use SnapPea's
        original algorithm, which is based on Meyerhoff-Hodgson-Neumann.
        
        Note: When computing the Chern-Simons invariant of a closed
        manifold, one must sometimes compute it first for the unfilled
        manifold so as to initialize SnapPea's internals.  For instance,

        >>> M = Manifold('5_2')
        >>> M.chern_simons()
        -0.15320413
        >>> M.dehn_fill( (1,2) )
        >>> M.chern_simons()
        0.07731787

        works, but will fail with 'Chern-Simons invariant not
        currently known' if the first call to chern_simons is not
        made.
        """
        return self._number_(self._chern_simons())
            
    def without_hyperbolic_structure(self):
        """
        Returns self as a Triangulation, forgetting the hyperbolic
        structure in the process.

        >>> M = Manifold('9_42')
        >>> T = M.without_hyperbolic_structure()
        >>> hasattr(T, 'volume')
        False
        """
        return Triangulation_from_Manifold(self)

    def _polish_hyperbolic_structures(self):
        polish_hyperbolic_structures(self.c_triangulation)

    def _two_to_three(self, tet_num, face_index):
        result = Triangulation._two_to_three(self, tet_num, face_index)
        polish_hyperbolic_structures(self.c_triangulation)
        return result

    def tetrahedra_shapes(self, part=None, fixed_alignment=True,
                          bits_prec=None, dec_prec=None,
                          intervals=False):
        """
        Gives the shapes of the tetrahedra in the current solution to
        the gluing equations.  Returns a list containing one info object
        for each tetrahedron.  The keys are:

        - rect : the shape of the tetrahedron, as a point in the
          complex plane.

        - log : the log of the shape

        - accuracies: a list of the approximate accuracies of the
          shapes, in order (rect re, rect im, log re, log im)

        If the optional variable 'part' is set to one of the above,
        then the function returns only that component of the data.
        
        If the flag 'fixed_alignment' is set to False, then the edges
        used to report the shape parameters are choosen so as to
        normalize the triangle.

        >>> M = Manifold('m015')
        >>> M.tetrahedra_shapes(part='rect')
        [0.66235898 + 0.56227951*I, 0.66235898 + 0.56227951*I, 0.66235898 + 0.56227951*I]
        >>> M.tetrahedra_shapes() #doctest:+SKIP
        [{'accuracies': (11, 11, 12, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I},
         {'accuracies': (11, 11, 11, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I},
         {'accuracies': (11, 11, 11, 11), 'log': -0.14059979 + 0.70385772*I, 'rect': 0.66235898 + 0.56227951*I}]
        """        
        cdef Real rect_re, rect_im, log_re, log_im
        cdef int acc_rec_re, acc_rec_im, acc_log_re, acc_log_im
        cdef Boolean is_geometric
        
        if self.c_triangulation is NULL: return []

        # The SnapPea kernel itself supports non-orientable manifolds,
        # but the extensions in snap and verify don't.
        # Thus work in double cover and use every other shape.
        if not self.is_orientable() and (bits_prec or dec_prec or intervals):
            result = self.orientation_cover().tetrahedra_shapes(
                                       part=part,
                                       fixed_alignment=fixed_alignment,
                                       bits_prec=bits_prec, dec_prec=dec_prec,
                                       intervals=intervals)[::2]
            if part != None:
                return result
            else:
                return ListOnePerLine(result)

        result = []
        if bits_prec or dec_prec:
            if fixed_alignment == False:
                raise ValueError(
                    'High precision shapes must be computed '
                    'in the fixed alignment')
            shapes = snap.polished_tetrahedra_shapes(
                self, dec_prec=dec_prec, bits_prec=bits_prec,
                ignore_solution_type=True)
            for z in shapes:
                result.append(ShapeInfo(rect=z,
                                        log=z.log(),
                                        accuracies=(None, None, None, None)))
        else:
            for i in range(self.num_tetrahedra()):
                get_tet_shape(self.c_triangulation, i, filled, fixed_alignment,
                              &rect_re, &rect_im, &log_re, &log_im,
                              &acc_rec_re, &acc_rec_im,
                              &acc_log_re, &acc_log_im,
                              &is_geometric)

                rect_shape=Number(RealImag2gen(rect_re, rect_im))
                rect_shape.accuracy=min(acc_rec_re, acc_rec_im)
                log_shape=Number(RealImag2gen(log_re, log_im))
                log_shape.accuracy=min(acc_log_re, acc_log_im)
                if part in ['rect', 'log', 'accuracies']:
                    # Don't compute volumes as they will be thrown away. 
                    result.append(ShapeInfo(
                            rect=self._number_(rect_shape),
                            log=self._number_(log_shape),
                            accuracies=(acc_rec_re, acc_rec_im,
                                        acc_log_re, acc_log_im)))
                else:
                    result.append(ShapeInfo(
                            rect=self._number_(rect_shape),
                            log=self._number_(log_shape),
                            volume=rect_shape.volume(),
                            accuracies=(acc_rec_re, acc_rec_im,
                                        acc_log_re, acc_log_im)))

        if intervals:
            if bits_prec or dec_prec:
                engine = verify.CertifiedShapesEngine(
                    self, [a['rect'] for a in result],
                    dec_prec=dec_prec, bits_prec=bits_prec)
            else:
                engine = verify.CertifiedShapesEngine(
                    self, [a['rect'] for a in result],
                    bits_prec = Number._default_precision)
            if not engine.expand_until_certified():
                raise RuntimeError('Could not certify shape intervals, either '
                                   'there are degenerate shapes or the '
                                   'precision must be increased.')
            result = [ ShapeInfo(rect=z,
                                 log=z.log(),
                                 accuracies=(None,None,None,None))
                       for z in engine.certified_shapes ]

        if part != None:
            try:
               return [a[part] for a in result]
            except KeyError:
                raise ValueError('A non-existent shape data type '
                                 'was specified.')
        else:
           return ListOnePerLine(result)

    def _get_tetrahedra_shapes(self, which_solution):
        """
        Return a list of the tetrahedra shapes from one of the two
        solutions maintained in the SnapPea triangulation structure.
        The which_solution argument must take one of the values
        'filled' or 'complete'.  This sets fixed_alignment = True.
        """
        cdef int i
        cdef Real rect_re, rect_im, log_re, log_im
        cdef int acc_rec_re, acc_rec_im, acc_log_re, acc_log_im
        cdef Boolean is_geometric
        cdef c_FillingStatus soln

        if which_solution == 'filled':
            soln = filled
        elif which_solution == 'complete':
            soln = complete
        else:
            raise ValueError("Please specify 'filled' or 'complete'.")
        result = []
        for i in range(self.num_tetrahedra()):
            get_tet_shape(self.c_triangulation, i, soln, True,
                          &rect_re, &rect_im, &log_re, &log_im,
                          &acc_rec_re, &acc_rec_im, &acc_log_re, &acc_log_im,
                          &is_geometric)
            result.append(Number(RealImag2gen(rect_re, rect_im)))
        return result

    def set_tetrahedra_shapes(self,
                              filled_shapes=None,
                              complete_shapes=None,
                              fillings=None, ):
        """
        Replaces the tetrahedron shapes with those in the given lists,
        and sets the Dehn filling coefficients as specified by the
        fillings argument.  The shapes will get double precision
        values; polishing will be needed for high precision shapes.
        """
        cdef int i, N
        cdef Complex *filled_shape_array = NULL
        cdef Complex *complete_shape_array = NULL

        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        N = get_num_tetrahedra(self.c_triangulation)
        if filled_shapes is not None:
            filled_shape_array = <Complex *>malloc(N*sizeof(Complex))
            for i, shape in enumerate(filled_shapes):
                filled_shape_array[i] = complex2Complex(complex(shape))
        if complete_shapes is not None:
            complete_shape_array = <Complex *>malloc(N*sizeof(Complex))
            for i, shape in enumerate(complete_shapes):
                complete_shape_array[i] = complex2Complex(complex(shape))
        if fillings:
            set_cusps(self.c_triangulation, fillings)
        set_tet_shapes(self.c_triangulation, 
                       filled_shape_array,
                       complete_shape_array)
        compute_holonomies(self.c_triangulation)
        compute_edge_angle_sums(self.c_triangulation)
        if filled_shape_array != NULL:
            free(filled_shape_array)
        if complete_shape_array != NULL:
            free(complete_shape_array)
        self._clear_cache(message='Manifold.set_tetrahedra_shapes')

    def solution_type(self, enum=False):
        """
        Returns the type of the current solution to the gluing
        equations, basically a summary of how degenerate the solution
        is.  If the flag enum=True is set, then an integer value is
        returned. The possible answers are:

        - 0: 'not attempted'
        
        - 1: 'all tetrahedra positively oriented' aka 'geometric_solution'
          Should correspond to a genuine hyperbolic structure.

        - 2: 'contains negatively oriented tetrahedra' aka 'nongeometric_solution'
          Probably correponds to a hyperbolic structure but some
          simplices have reversed orientiations.  
             
        - 3: 'contains flat tetrahedra' All tetrahedra have shape in R - {0, 1}.

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
        c_target = Object2Complex(target)
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
        0.11044502 + 0.94677098*I
        >>> c.modulus
        -0.12155872 + 1.04204128*I
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
        (-0.59883089 + 1.09812548*I, 0.89824633 + 1.49440443*I)

        The complex numbers returned for the shape and for the two
        holonomies have an extra attribute, accuracy, which is
        SnapPea's *estimate* of their accuracy.
        
        You can also get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : complete torus cusp of shape 0.11044502 + 0.94677098*I,
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('is_complete')
        [True, False, False]
        """
        cdef int cusp_index
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef Real m, l
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
        get_holonomy(self.c_triangulation, cusp_index,
                     &c_meridian, &c_longitude,
                     &meridian_accuracy, &longitude_accuracy)
        shape= Complex2Number(current_shape)
        shape.accuracy = current_shape_accuracy
        meridian = Complex2Number(c_meridian)
        meridian.accuracy = meridian_accuracy
        longitude = Complex2Number(c_longitude)
        longitude.accuracy = longitude_accuracy
        modulus = Complex2Number(current_modulus)
        info = {
            'index' : cusp_index,
            'topology' : CuspTopology[topology],
            'is_complete' : B2B(is_complete),
            'filling' : (Real2float(m), Real2float(l)),
            'shape': self._number_(shape),
            'shape_accuracy': current_shape_accuracy,
            'modulus': self._number_(modulus),
            'holonomies':(self._number_(meridian), self._number_(longitude)),
            'holonomy_accuracy':min(meridian_accuracy,longitude_accuracy)
        }

        core_geodesic(self.c_triangulation, cusp_index,
                      &singularity_index, &c_core_length, &accuracy)
            
        if singularity_index != 0:
            core_length = Complex2Number(c_core_length)
            core_length.accuracy = accuracy
            info.update({
                'core_length' : self._number_(core_length),
                'singularity_index':singularity_index
            })
                
        return CuspInfo(**info)
                
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
        self._clear_cache(message='Manifold.dehn_fill')

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
          [-2.49024467 + 2.97944707*I]
          >>> cob = M.set_peripheral_curves('shortest', return_matrices=True)
          >>> M.cusp_info('shape')
          [-0.49024467 + 2.97944707*I]
          >>> cob
          [[[1, 0], [-2, 1]]]
          >>> M.set_peripheral_curves(cob)
          >>> M.cusp_info('shape')
          [-2.49024467 + 2.97944707*I]

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

        self._clear_cache(message='Manifold.set_peripheral_curves')

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

        else:
            return Triangulation.set_peripheral_curves(
                self, peripheral_data, which_cusp,return_matrices)
        
    def dual_curves(self, max_segments=6):
        """
        Constructs a *reasonable* selection of simple closed curves in
        a manifold's dual 1-skeleton.  In particular, it returns thos e
        that appear to represent geodesics. The resulting curves can
        be drilled out.

        >>> M = Manifold('m015')
        >>> curves = M.dual_curves()
        >>> curves
        [  0: orientation-preserving curve of length 0.56239915 - 2.81543089*I,
           1: orientation-preserving curve of length 1.12479830 + 0.65232354*I,
           2: orientation-preserving curve of length 1.26080402 + 1.97804689*I,
           3: orientation-preserving curve of length 1.58826933 + 1.67347167*I,
           4: orientation-preserving curve of length 1.68719745 + 2.81543089*I]

        Each curve is returned as an info object with these keys
        
        >>> sorted(curves[0].keys())
        ['complete_length', 'filled_length', 'index', 'max_segments', 'parity']
        
        We can drill out any of these curves to get a new manifold
        with one more cusp.

        >>> N = M.drill(curves[0])
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        
        By default, this function only finds curves of length 6; this
        can be changed by specifying the optional argument
        max_segments

        >>> M.dual_curves(max_segments=2)
        [  0: orientation-preserving curve of length 0.56239915 - 2.81543089*I]
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
                    index=i,
                    parity=parity,
                    filled_length=self._number_(Complex2Number(filled_length)), 
                    complete_length=self._number_(Complex2Number(complete_length)),
                    max_segments=max_segments
                  )
               )
        free_dual_curves(num_curves, curve_list)
        return ListOnePerLine(result)
    
    def length_spectrum(self, cutoff=1.0, full_rigor=True):
        """
        M.length_spectrum(cutoff=1.0)

        Returns a list of geodesics (with multiplicities) of length
        up to the specified cutoff value. (The default cutoff is 1.0.)
        """
        D = self.dirichlet_domain()
        name_mangled = 'length_spectrum:%f'%cutoff
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = D.length_spectrum_dicts(
                cutoff_length=cutoff,
                full_rigor=full_rigor)
        return self._cache[name_mangled]
        
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
            result = self.__class__('empty')
            result.set_c_triangulation(c_triangulation)
            return result

    def is_isometric_to(self, Manifold other, return_isometries=False):
        """
        Returns True if M and N are isometric, False if they not.  A
        RuntimeError is raised in cases where the SnapPea kernel fails
        to determine either answer.  (This is fairly common for closed
        manifolds.)

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

        try:
            if self.homology() != other.homology():
                if return_isometries:
                    return []
                else:
                    return False
        except ValueError:
            pass
        

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
            raise RuntimeError('The SnapPea kernel was not able to '
                               'determine if the manifolds are isometric.')
        
        ans = bool(are_isometric)

        if return_isometries:
            if not ans:
                return []
            else:
                if isometries is NULL:  # means manifold is closed
                    raise ValueError("Can't get the list of isometries for closed manifolds")
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
        free_triangulation(c_canonized_triangulation)
        if is_two_bridge:
            return (p,q)
        else:
            return False

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
        cdef int neighbor0_idx, neighbor1_idx, neighbor2_idx, neighbor3_idx
        cdef int perm0, perm1, perm2, perm3	

        ans = []
        for i in range(self.num_tetrahedra()):
            choose_gen_tetrahedron_info(self.c_triangulation,
                                        i, &generator_path,
                                        &face0_gen, &face1_gen,
                                        &face2_gen, &face3_gen,
                                        &c0, &c1, &c2, &c3,
                                        &neighbor0_idx, &neighbor1_idx,
                                        &neighbor2_idx, &neighbor3_idx,
                                        &perm0, &perm1, &perm2, &perm3)
            ans.append(
                {'index':i,
                 'generators':(face0_gen, face1_gen, face2_gen, face3_gen),
                 'neighbors':(neighbor0_idx, neighbor1_idx,
                              neighbor2_idx, neighbor3_idx),
                 'gluings': ( tuple([ perm0>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm1>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm2>>(2 * i) & 3 for i in range(4)]),
                              tuple([ perm3>>(2 * i) & 3 for i in range(4)])),
                 'corners': ( self._number_(Complex2Number(c0)),
                              self._number_(Complex2Number(c1)),
                              self._number_(Complex2Number(c2)),
                              self._number_(Complex2Number(c3)) ),
                 'generator_path':generator_path}
                )
        return ans

    def splitting_surfaces(self):
        """
        Searches for connected closed normal surfaces of nonnegative Euler
        characteristic.  If spheres or projective planes are found, then
        tori and Klein bottles aren't reported.  There is no guarantee
        that all such normal surfaces will be found nor that any given
        surface is incompressible.  The search is confined to surfaces
        whose quads are in the tetrahedra that have degenerate shapes.
        
        You can split the manifold open along one of these surfaces
        using the method "split".
        
        A connect sum of two trefoils:

        >>> M1 = Manifold('DT: fafBCAEFD')
        >>> len(M1.splitting_surfaces())
        2

        First satellite knot in the table. 

        >>> M2 = Manifold('K13n4587')
        >>> M2.splitting_surfaces()
        [Orientable two-sided with euler = 0]
        """
        cdef NormalSurfaceList *surfaces
        cdef c_FuncResult result
                
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        result = find_normal_surfaces(self.c_triangulation, &surfaces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when computing normal surfaces')

        ans = []
        for i in range(number_of_normal_surfaces_on_list(surfaces)):
            ans.append( NormalSurfaceInfo(
                orientable=B2B(normal_surface_is_orientable(surfaces, i)),
                two_sided=B2B(normal_surface_is_two_sided(surfaces, i)),
                euler=normal_surface_Euler_characteristic(surfaces, i),
                index=i))
        
        free_normal_surfaces(surfaces)        
        return ListOnePerLine(ans)

    def split(self, which_surface):
        """
        Split the manifold open along a surface of positive characteristic found
        by the method "splitting_surfaces".  Returns a list of the pieces, with any 
        sphere boundary components filled in.
        
        Here's an example of a Whitehead double on the trefoil.        

        >>> M = Manifold('K14n26039')
        >>> S = M.splitting_surfaces()[0]
        >>> S
        Orientable two-sided with euler = 0

        >>> pieces = M.split(S); pieces
        [K14n26039.a(0,0)(0,0), K14n26039.b(0,0)]
        >>> pieces[0].volume()
        3.66386238
        >>> pieces[1].fundamental_group().relators()
        ['aabbb']
                
        You can also specify a surface by its index.

        >>> M = Manifold('L10n111') 
        >>> max( P.volume() for P in M.split(0) )
        5.33348957
        """
        cdef NormalSurfaceList *surfaces
        cdef c_FuncResult result
        cdef int num_surfaces, i
        cdef c_Triangulation  *pieces[2]
        cdef Manifold M0, M1
        
        if self.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')

        result = find_normal_surfaces(self.c_triangulation, &surfaces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when computing normal surfaces.')

        if isinstance(which_surface, NormalSurfaceInfo):
            which_surface = which_surface.index

        num_surfaces = number_of_normal_surfaces_on_list(surfaces)
        if not (0 <= which_surface < num_surfaces):
            raise ValueError('SnapPea only found %d surfaces, but you asked for surface number %s.' % (num_surfaces, which_surface))

        result =  split_along_normal_surface(surfaces, which_surface, pieces)
        if result != func_OK:
            raise RuntimeError('SnapPea kernel failed when splitting open along the given surface.')
        
        ans = []
        if pieces[0] != NULL:
            M0 = self.__class__('empty')
            M0.set_c_triangulation(pieces[0])
            M0.set_name(self.name() + '.a')
            ans.append(M0)
        if pieces[1] != NULL:
            M1 = self.__class__('empty')
            M1.set_c_triangulation(pieces[1])
            M1.set_name(self.name() + '.b')
            ans.append(M1)

        free_normal_surfaces(surfaces)
        return ans

    def identify(self, extends_to_link=False):
        """
        Look for the manifold in all of the SnapPy databases:

        >>> M = Manifold('m125')
        >>> M.identify()
        [m125(0,0)(0,0), L13n5885(0,0)(0,0)]
        
        One can require that there be an isometry taking merdians
        to meridians:

        >>> M.identify(extends_to_link=True)
        [m125(0,0)(0,0)]
        
        For closed manifolds, extends_to_link doesn't make sense because
        of how the kernel code works:
        
        >>> C = Manifold("m015(1,2)")
        >>> C.identify()
        [m006(-5,2)]
        >>> C.identify(True)
        []
        """
        ans = []
        for table in database.__all_tables__:
            match = table.identify(self, extends_to_link)
            if match:
                ans.append(match)
        return ans

    def _cusp_cross_section_info(self):
        cdef c_Tetrahedron *tet
        cdef Real temp
        allocate_cross_sections(self.c_triangulation)
        compute_cross_sections(self.c_triangulation)
        compute_tilts(self.c_triangulation)
        tilts, side_lengths = [], []
        tet = self.c_triangulation.tet_list_begin.next
        while tet != &self.c_triangulation.tet_list_end:
            one_tet_tilts, one_tet_lengths = [], []
            for v in range(4):
                one_tet_tilts.append(self._number_(Real2Number(<Real>tet.tilt[v])))
                one_vertex_lengths = []
                for f in range(4):
                    if v != f:
                         one_vertex_lengths.append(self._number_(
                             Real2Number(<Real>tet.cross_section.edge_length[v][f])))
                    else:
                        one_vertex_lengths.append(None)
                one_tet_lengths.append(one_vertex_lengths)
            tilts.append(one_tet_tilts)
            side_lengths.append(one_tet_lengths)
            tet = tet.next
                
        free_cross_sections(self.c_triangulation)
        return tilts, side_lengths

# PLink communication

def _plink_callback(LE):
    cdef Manifold manifold
    cdef c_Triangulation* c_triangulation = NULL
    if LE.manifold is None:
        LE.manifold = Manifold('empty')
    manifold = LE.manifold
    klp = LE.SnapPea_KLPProjection()
    if klp is not None:
        manifold._set_DTcode(spherogram.DTcodec(*LE.DT_code()))
        c_triangulation = get_triangulation_from_PythonKLP(klp)
        if c_triangulation is not NULL:
            find_complete_hyperbolic_structure(c_triangulation)
        if manifold.c_triangulation is not NULL:
            free_triangulation(manifold.c_triangulation)
        manifold.set_c_triangulation(c_triangulation)
        manifold._clear_cache(message='plink_callback')
        msg_stream.write('\nNew triangulation received from PLink!\n')
    else:
        raise RuntimeError('Communication with PLink failed.')


# Conversion functions Manifold <-> Triangulation

def Manifold_from_Triangulation(Triangulation T, recompute=True,
                                manifold_class=None):
    cdef c_Triangulation *c_triangulation
    cdef Manifold M

    M = _manifold_class('empty') if manifold_class is None else manifold_class('empty')
    if T.c_triangulation is NULL:
        return M
    copy_triangulation(T.c_triangulation, &c_triangulation)
    M.set_c_triangulation(c_triangulation)
    if recompute:
        if mark_fake_cusps(c_triangulation):
            remove_finite_vertices(c_triangulation)
        find_complete_hyperbolic_structure(c_triangulation)
        do_Dehn_filling(c_triangulation)
    M.set_name(T.name())
    M._cover_info = T._cover_info
    return M

def Triangulation_from_Manifold(Manifold M):
    cdef c_Triangulation *c_triangulation
    cdef Triangulation T

    if M.c_triangulation is NULL:
        return _triangulation_class('empty')

    copy_triangulation(M.c_triangulation, &c_triangulation)
    remove_hyperbolic_structures(c_triangulation)
    T = _triangulation_class('empty')
    T.set_c_triangulation(c_triangulation)
    T.set_name(M.name())
    T._cover_info = M._cover_info
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
    if not verbose_form:
        return word
    if re.search('\d', word):
        letters = re.findall('([xX]\d+)', word)
    else:
        letters = list(word)
    return '*'.join([a if a[0].islower() else a.lower() + '^-1' for a in letters])

cdef class CFundamentalGroup:
    cdef c_GroupPresentation *c_group_presentation
    cdef c_Triangulation *c_triangulation
    cdef readonly num_cusps

    cdef int_to_gen_string(self, int g):
        if self.num_generators() <=26:
            return Alphabet[g]
        else:
            ans = 'x' if g > 0 else 'X'
            return ans + '%d' % abs(g)

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
            word_list.append(self.int_to_gen_string(word[n]))
            n += 1
        return ''.join(word_list)

    cdef int *c_word_from_list(self, word_list):
        cdef int *c_word
        cdef int length, size, n
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
                  minimize_number_of_generators = True,
                  try_hard_to_shorten_relators = True):
        if triangulation.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        copy_triangulation(triangulation.c_triangulation,
                           &self.c_triangulation)
        self.c_group_presentation = fundamental_group(
            self.c_triangulation,
            simplify_presentation,
            fillings_may_affect_generators,
            minimize_number_of_generators,
            try_hard_to_shorten_relators)
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
        if len(generators) > 26:
            word = re.findall('[xX]\d+', word)
        for letter in word:
            try:
                if letter[0].islower():
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
        n = self.num_generators()
        return [ self.int_to_gen_string(i) for i in range(1, 1+n) ]

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

    def sage(self):
        """
        Returns the corresponding Sage FinitelyPresentedGroup
        """
        if not _within_sage:
            raise RuntimeError("Not within Sage")
        F = FreeGroup(self.generators())
        rels = [F(R) for R in self.relators(as_int_list=True)]
        return F/rels

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
        cdef int i, j
        word_list = self._word_as_list(word)
        c_word = self.c_word_from_list(word_list)
        result = fg_word_to_matrix(self.c_group_presentation, c_word, O, &M)
        if result == 0:
            sl2 = matrix(
                [[self._number_(Complex2Number(M.matrix[i][j]))
                  for j in range(2)] for i in range(2)] )
            o31 = matrix(
                [[self._number_(Real2Number(<Real>O[i][j]))
                  for j in range(4)] for i in range(4)] )
            L = self._number_(Complex2Number(complex_length_mt(&M)))
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
    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        cls._number_ = staticmethod(func)

if _within_sage:
    HolonomyGroup.__bases__ += (sage.structure.sage_object.SageObject,)
    HolonomyGroup.use_field_conversion(lambda n : n.sage())

# Dirichlet Domains

cdef WEPolyhedron* read_generators_from_file(
    file_name,
    double vertex_epsilon=default_vertex_epsilon):
    
    data = open(file_name).readlines()
    if data[0].strip() != '% Generators':
        raise ValueError('The generator file does not start with '
                         '"% Generators"')
    nums = []
    for line in data[1:]:
        nums +=  line.split()
    num_gens = int(nums[0])
    nums.pop(0)

    cdef O31Matrix *generators
    cdef MoebiusTransformation *temp_gens
    if len(nums) == 16 * num_gens:
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            for j in range(4):
                for k in range(4):
                    num_string = nums.pop(0) # save a reference
                    generators[i][j][k] =  <Real_struct>Real_from_string(
                        <char*>num_string)
    elif len(nums) == 8*num_gens:
        temp_gens = <MoebiusTransformation *>malloc(
            num_gens*sizeof(MoebiusTransformation))
        generators = <O31Matrix *>malloc(num_gens*sizeof(O31Matrix))
        for i in range(num_gens):
            temp_gens[i].parity = orientation_preserving
            for j in range(2):
                for k in range(2):
                    num_string = nums.pop(0) # save a reference
                    temp_gens[i].matrix[j][k].real = Real_from_string(
                        <char*>num_string)
                    num_string = nums.pop(0) # save a reference
                    temp_gens[i].matrix[j][k].imag = Real_from_string(
                        <char*>num_string)
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
 
    if not O31_determinants_OK(generators, num_gens, det_error_epsilon):
        raise ValueError('The data given do not have the '
                         'right determinants.')
        
    cdef WEPolyhedron *dirichlet_domain
    dirichlet_domain = Dirichlet_from_generators(generators,
                                                 num_gens,
                                                 vertex_epsilon,
                                                 Dirichlet_keep_going,
                                                 True);
    free(generators)
    return dirichlet_domain

cdef WEPolyhedron* dirichlet_from_O31_matrix_list(
    matrices,
    double vertex_epsilon=default_vertex_epsilon)except*:

    cdef WEPolyhedron* c_dirichlet_domain
    cdef O31Matrix* generators
    cdef int i, j, k, num_gens
    num_gens = len(matrices)
    generators = <O31Matrix*>malloc(num_gens*sizeof(O31Matrix))
    for i, A in enumerate(matrices):
        for j in range(4):
            for k in range(4):
                generators[i][j][k] = <Real_struct>Object2Real(A[j,k])
    if not O31_determinants_OK(generators, num_gens, det_error_epsilon):
        raise ValueError('The data given do not have the '
                         'right determinants.')
    c_dirichlet_domain = Dirichlet_from_generators(generators,
                                                   num_gens,
                                                   vertex_epsilon,
                                                   Dirichlet_keep_going,
                                                   True);
    free(generators)
    return c_dirichlet_domain

cdef class CDirichletDomain:
    cdef WEPolyhedron *c_dirichlet_domain
    cdef c_Triangulation *c_triangulation

    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        cls._number_ = staticmethod(func)

    def __cinit__(self, 
                  Manifold manifold=None,
                  vertex_epsilon=default_vertex_epsilon,
                  displacement = [0.0, 0.0, 0.0],
                  centroid_at_origin=True,
                  maximize_injectivity_radius=True,
                  generator_file = None,
                  O31_generators = None):
        cdef double c_displacement[3]
        self.c_dirichlet_domain = NULL
        if generator_file != None:
            self.c_dirichlet_domain = read_generators_from_file(
                generator_file)
            self.manifold_name = generator_file
        elif O31_generators != None:
            self.c_dirichlet_domain = dirichlet_from_O31_matrix_list(
                O31_generators)
            self.manifold_name = 'unnamed'
        else:
            if manifold is None:
                raise ValueError('Supply a manifold, or generators.')
            if manifold.c_triangulation is NULL:
                raise ValueError('The Triangulation is empty.')
            for n from 0 <= n < 3:
                c_displacement[n] = <double>displacement[n] 
            copy_triangulation(manifold.c_triangulation,
                               &self.c_triangulation)
            self.c_dirichlet_domain = Dirichlet_with_displacement(
                self.c_triangulation,
                c_displacement,
                vertex_epsilon,
                centroid_at_origin,
                Dirichlet_keep_going,
                maximize_injectivity_radius )
            self.manifold_name = manifold.name()
        if self.c_dirichlet_domain == NULL:
            raise RuntimeError('The Dirichlet construction failed.')
        if manifold:
            self._number_ = manifold._number_

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
        radius = Real2Number(self.c_dirichlet_domain.inradius)
        return self._number_(radius)

    def out_radius(self):
        """
        Return the radius of the smallest circubscribed sphere.
        """
        radius = Real2Number(self.c_dirichlet_domain.outradius)
        return self._number_(radius)

    def spine_radius(self):
        """
	Return the infimum of the radii (measured from the origin) of all 
	spines dual to the Dirichlet domain.
	"""
        radius = Real2Number(self.c_dirichlet_domain.spine_radius)
        return self._number_(radius)

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
                        Object2Real(cutoff_length),
                        full_rigor,
                        multiplicities,
                        Object2Real(user_radius),
                        &geodesics,
                        &num_lengths)
        spectrum = []
        for n from 0 <= n < num_lengths:
            length = Complex2Number(geodesics[n].length)
            spectrum.append(
               LengthSpectrumInfo(
                  length=self._number_(length),
                  parity=MatrixParity[geodesics[n].parity],
                  topology=Orbifold1[geodesics[n].topology],
                  multiplicity=geodesics[n].multiplicity
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
          vertices.append(
              ( self._number_(Real2Number(<Real>vertex.x[1])),
                self._number_(Real2Number(<Real>vertex.x[2])),
                self._number_(Real2Number(<Real>vertex.x[3])) ) )
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
                vertices.append(tuple( 
                    self._number_(Real2Number(<Real>vertex.x[i]))
                    for i in range(1,4) ))
                # get the next edge
                if edge.f[left] == face:
                    edge = edge.e[tip][left];
                else:
                    edge = edge.e[tail][right];
                if edge == face.some_edge:
                    break
            faces.append(
                {'vertices' : vertices,
                 'distance' : self._number_(Real2Number(<Real>face.dist)),
                 'closest'  : [
                     self._number_(Real2Number(<Real>face.closest_point[i]))
                     for i in range(1,4) ],
                 'hue'      : Real2double(face.f_class.hue) })
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

    def manifold(self):
        """
        Returns a Manifold computed directly from the Dirichlet
        domain, regarded as polyhedron with faces identified in pairs.
        Only works if this gives a manifold not an orbifold.

        >>> M = Manifold('7_3')
        >>> D = M.dirichlet_domain()
        >>> A = D.manifold()
        >>> M.is_isometric_to(A)
        True

        """
        return self.triangulate(_manifold_class)

    def triangulation(self):
        """
        Returns a Triangulation computed directly from the Dirichlet
        domain, regarded as polyhedron with faces identified in pairs.
        Only works if this gives a manifold not an orbifold.

        >>> M = Manifold('7_3')
        >>> D = M.dirichlet_domain()
        >>> B = D.triangulation()
        >>> M.is_isometric_to(B.with_hyperbolic_structure())
        True

        """
        return self.triangulate(_triangulation_class)

    cdef triangulate(self, return_class):
        cdef c_Triangulation *c_triangulation
        cdef Triangulation M
        c_triangulation = Dirichlet_to_triangulation(self.c_dirichlet_domain)
        if c_triangulation is NULL:
            raise ValueError('The Dirichlet domain could not be '
                             'triangulated; perhaps this is an '
                             'orbifold group?')
        M = return_class('empty')
        M.set_c_triangulation(c_triangulation)
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
        volume = Real2Number(self.c_dirichlet_domain.approximate_volume)
        return self._number_(volume)

    def pairing_matrices(self):
        """
        Returns a list of the O31Matrices which pair the faces of
        this DirichletDomain.

        >>> M = Manifold('s345')
        >>> D = M.dirichlet_domain()
        >>> matrices = D.pairing_matrices()
        >>> D1 = DirichletDomain(O31_generators=matrices)
        >>> N = D1.manifold()
        >>> M.is_isometric_to(N)
        True
        """
        cdef O31Matrix* M
        cdef int i, j
        cdef WEFace* face

        if self.c_dirichlet_domain == NULL:
            raise ValueError('The Dirichlet Domain was not computed.')
        matrices = []
        face = self.c_dirichlet_domain.face_list_begin.next
        while face != &self.c_dirichlet_domain.face_list_end:
            M = face.group_element
            matrices.append(matrix(
                [[self._number_(Real2Number(<Real>M[0][i][j]))
                  for j in range(4)] for i in range(4)] ))
            face = face.next
        return matrices

    def save(self, filename):
        """
        Save the Dirichlet domain as a text file in "% Generators" format.
        """
        matrices = self.pairing_matrices()
        with open(filename, 'wb') as output:
            output.write('% Generators\n')
            output.write('%s\n'%len(matrices))
            for matrix in matrices:
                for row in matrix:
                    output.write(' %s\n'%' '.join([str(x) for x in row]))
                output.write('\n')

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

# Cusp Neighborhoods

cdef class CCuspNeighborhood:
    cdef c_CuspNeighborhoods *c_cusp_neighborhood
    cdef c_Triangulation *c_triangulation
    cdef int _num_cusps
    cdef original_indices

    @staticmethod
    def _number_(number):
        return number

    @classmethod
    def use_field_conversion(cls, func):
        cls._number_ = staticmethod(func)

    def __cinit__(self, Manifold manifold):
        if manifold.c_triangulation is NULL:
            raise ValueError('The Triangulation is empty.')
        is_complete = manifold.cusp_info('is_complete')
        self.original_indices = [n for n, c in enumerate(is_complete) if c]
        copy_triangulation(manifold.c_triangulation,
                           &self.c_triangulation)
        self.c_cusp_neighborhood = initialize_cusp_neighborhoods(
            self.c_triangulation)
        if self.c_cusp_neighborhood == NULL:
            raise RuntimeError('The cusp neighborhood construction failed.')
        self.manifold_name = manifold.name()
        self._num_cusps = get_num_cusp_neighborhoods(self.c_cusp_neighborhood)
        self._number_ = manifold._number_

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_cusp_neighborhood != NULL:
            free_cusp_neighborhoods(self.c_cusp_neighborhood)

    def __repr__(self):
        N = self._num_cusps
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

    def original_index(self, which_cusp):
        """
        Returns the index by which the Manifold identifies this cusp.
        """
        return self.original_indices[which_cusp]

    def check_index(self, which_cusp):
        """
        Raises an IndexError if the cusp index is invalid.
        """
        N = int(which_cusp)
        if 0 <= N < self._num_cusps:
            return N
        else:
            raise IndexError('The specified cusp (%s) does not '
                             'exist.'%which_cusp)
        
    def num_cusps(self):
        """
        Return the number of cusps.
        """
        return self._num_cusps

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
        disp =  Number(Real2gen(get_cusp_neighborhood_displacement(
                self.c_cusp_neighborhood, N)))
        return self._number_(disp)

    def set_displacement(self, new_displacement, which_cusp=0):
        """
        Set the displacement of the specified cusp.
        """
        N = self.check_index(which_cusp)
        set_cusp_neighborhood_displacement(self.c_cusp_neighborhood,
                                           N, Object2Real(new_displacement))

    def stopping_displacement(self, which_cusp=0):
        """
        Return the displacement at which the specified cusp
        neighborhood bumps into itself or another cusp neighborhood.
        (Assumes the other displacements are fixed.)
        """
        disp = Number(Real2gen(get_cusp_neighborhood_stopping_displacement(
            self.c_cusp_neighborhood, which_cusp)))
        return self._number_(disp)

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
        N = self.check_index(which_cusp)
        reach = Number(Real2gen(get_cusp_neighborhood_reach(
            self.c_cusp_neighborhood, N)))
        return self._number_(reach)

    def max_reach(self):
        """
        Return the maximum reach over all cusps.
        """
        reach = Number(Real2gen(get_cusp_neighborhood_max_reach(
            self.c_cusp_neighborhood)))
        return self._number_(reach)

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
        volume = Number(Real2gen(get_cusp_neighborhood_cusp_volume(
                self.c_cusp_neighborhood, N)))
        return self._number_(volume)

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
        M, L = Complex2Number(meridian), Complex2Number(longitude)
        return self._number_(M), self._number_(L)

    def horoballs(self, cutoff=0.1, which_cusp=0, full_list=True,
                  high_precision=False):
        """
        Return a list of dictionaries describing the horoballs with
        height at least cutoff.  The keys are 'center', 'radius', 'index'.
        
        If the high_precision flag is set to the default value False, these
        are Python complexes and floats.  Otherwise they are SnapPy Numbers.
        """
        cdef CuspNbhdHoroballList* horoball_list
        cdef CuspNbhdHoroball ball
        horoball_list = get_cusp_neighborhood_horoballs(
            self.c_cusp_neighborhood,
            which_cusp,
            full_list,
            Object2Real(cutoff))
        if horoball_list == NULL:
            raise RuntimeError('The horoball construction failed.')
        result = []
        for n from 0 <= n < horoball_list.num_horoballs:
            ball = horoball_list.horoball[n]
            if high_precision:
                dict = {
                    'center' : self._number_(Complex2Number(ball.center)),
                    'radius' : self._number_(Real2Number(ball.radius)),
                    'index'  : ball.cusp_index}
            else:
                dict = {'center' : Complex2complex(ball.center),
                        'radius' : Real2float(ball.radius),
                        'index'  : ball.cusp_index}
            result.append(dict)
        free_cusp_neighborhood_horoball_list(horoball_list)
        return result

    def Ford_domain(self, which_cusp=0, high_precision=False):
        """
        Return a list of pairs of complex numbers describing the
        endpoints of the segments obtained by projecting the edges of
        the Ford domain to the xy-plane in the upper half space model.

        If the high_precision flag is set to False (the default), the
        coordinates are Python complex numbers.  Otherwise they are
        SnapPy Numbers.
        """
        cdef CuspNbhdSegmentList* segment_list
        cdef CuspNbhdSegment segment
        segment_list = get_cusp_neighborhood_Ford_domain(
            self.c_cusp_neighborhood,
            which_cusp)
        if segment_list == NULL:
            raise RuntimeError('The Ford domain construction failed.')
        result = []
        for n from 0 <= n < segment_list.num_segments:
            segment = segment_list.segment[n]
            if high_precision:
                pair = ( 
                    self._number_(Complex2Number(segment.endpoint[0])),
                    self._number_(Complex2Number(segment.endpoint[1])) )
            else:
                pair = ( Complex2complex(segment.endpoint[0]),
                         Complex2complex(segment.endpoint[1]) )
            result.append(pair)
        free_cusp_neighborhood_segment_list(segment_list)
        return result

    def triangulation(self, which_cusp=0, high_precision=False):
        """
        Return a list of dictionaries describing the endpoints of the
        segments obtained by projecting the edges of the triangulation
        dual to the Ford domain into the xy-plane in the upper half
        space model.  The keys are 'endpoints' and 'indices'.
        """
        cdef CuspNbhdSegmentList* segment_list
        cdef CuspNbhdSegment segment
        segment_list = get_cusp_neighborhood_triangulation(
            self.c_cusp_neighborhood,
            which_cusp)
        if segment_list == NULL:
            raise RuntimeError('The triangulation construction failed.')
        result = []
        for n from 0 <= n < segment_list.num_segments:
            segment = segment_list.segment[n]
            if high_precision:
                endpoints = (
                    self._number_(Complex2Number(segment.endpoint[0])),
                    self._number_(Complex2Number(segment.endpoint[1])) )
            else:
                endpoints = (Complex2complex(segment.endpoint[0]),
                             Complex2complex(segment.endpoint[1]))
            indices = (segment.start_index,
                       segment.middle_index,
                       segment.end_index)
            result.append({'endpoints' : endpoints, 'indices' : indices})
        free_cusp_neighborhood_segment_list(segment_list)
        return result

    def view(self, which_cusp=0, cutoff=None):
        """
        Create a 3D picture of the horoball packing.  One can specify
        which cusp to put at infinity and how large of horoballs to
        look at, e.g.

        >>> M = Manifold('m125')
        >>> C = M.cusp_neighborhood()
        >>> C.view(which_cusp = 1, cutoff=0.2)   #doctest: +CYOPENGL
        """
        if HoroballViewer:
            self.viewer = HoroballViewer(
                self, which_cusp=which_cusp, cutoff=cutoff,
                title='Cusp neighborhood%s of %s'%(
                    's' if self.num_cusps > 1 else '',
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
        elif self.order() == 1:
            return 'unknown'
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

split_filling_info = re.compile('(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))*$)')
is_census_manifold = re.compile('([msvtxy])([0-9]+)$|o9_\d\d\d\d\d$')
is_torus_bundle = re.compile('b([+-no])([+-])([lLrR]+)$')
is_knot_complement = re.compile('([0-9]+_[0-9]+)$')
is_link_complement1_pat = '(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$'
is_link_complement2_pat = '(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$'
is_link_complement3_pat = '[lL](?P<components>[0-9]{1})(?P<crossings>[0-9]{2})(?P<index>[0-9]+)$'
is_link_complement1 = re.compile(is_link_complement1_pat)
is_link_complement2 = re.compile(is_link_complement2_pat)
is_link_complement3 = re.compile(is_link_complement3_pat)
rolfsen_link_regexs = [is_link_complement1, is_link_complement2, is_link_complement3]
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)$')
is_HT_link = re.compile('[KL][0-9]+[an]([0-9]+)$')
is_braid_complement = re.compile('[Bb]raid[:]?(\[[0-9, \-]+\])$')
is_int_DT_exterior = re.compile('DT[:]? *(\[[0-9, \-\(\)]+\](?: *, *\[[01, ]+\])?)$')
is_alpha_DT_exterior = re.compile('DT[:\[] *([a-zA-Z]+(?:\.[01]+)?)[\]]?$')
is_census_knot = re.compile('[kK][2-8]_([0-9]+)$')
twister_word = '\[*([abcABC_\d!*]*|[abcABC_\d!,]*)\]*'
is_twister_bundle = re.compile('Bundle\(S_\{(\d+),(\d+)\},'+twister_word+'\)')
is_twister_splitting = re.compile('Splitting\(S_\{(\d+),(\d+)\},'+twister_word+','+twister_word+'\)')
is_isosig = re.compile('([a-zA-Z0-9\+\-]+)$')
is_decorated_isosig = decorated_isosig.isosig_pattern

# Hooks so that global module can monkey patch in modified versions
# of the Triangulation and Manifold classes.

_triangulation_class = Triangulation
_manifold_class = Manifold

def bundle_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_bundle.match(desc)
    if m:
        g, n, monodromy = m.groups()
        g, n = int(g), int(n)
        monodromy = monodromy.replace(',', '*').replace('_', '')
        surface = twister.Surface( (g, n) )
        return surface.bundle(monodromy, return_type='triangulation')
    
def splitting_from_string(desc):
    desc = desc.replace(' ', '')
    m = is_twister_splitting.match(desc)
    if m:
        g, n, gluing, handles = m.groups()
        g, n = int(g), int(n)
        gluing = gluing.replace(',', '*').replace('_', '')
        handles = handles.replace(',', '*').replace('_', '')
        surface = twister.Surface( (g, n) )        
        return surface.splitting(gluing, handles, return_type='triangulation')

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
                          is_complete,
                          Object2Real(meridian),
                          Object2Real(longitude))
    return 0

# Testing code for get_triangulation
def get_triangulation_tester():
    """
    >>> get_triangulation_tester()
    L13n9331(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    m003(0,0) 2.02988321 Z/5 + Z
    m004(0,0) 2.02988321 Z
    v1205(2,3) 4.70744340 Z/40
    x012(0,0)(0,0) 3.54972978 Z/2 + Z
    y123(0,0) 5.02755480 Z
    L13n9331(3,4)(2,3)(2,1) 14.60215339 Z/53
    K7_1(0,0) 3.57388254 Z
    6_1(0,0) 3.16396323 Z
    5^2_1(3,4)(1,-2) 2.73300075 Z/3
    8^3_3(0,0)(0,0)(0,0) 8.96736085 Z + Z + Z
    4_1(0,0) 2.02988321 Z
    12n123(0,0) 18.15036328 Z
    16n1235(0,0) 21.29383093 Z
    b++RL(0,0) 2.02988321 Z
    b-+RRL(0,0) 2.40690959 Z/3 + Z
    b+-RL(0,0) 2.02988321 Z/5 + Z
    b--RRL(0,0) 2.40690959 Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) 4.05976643 Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT[4,6,2](0,0) 0.0 Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) 16.64369585 Z + Z + Z
    a0*B1(0,0) 2.02988321 Z
    b1*A0 a0*B1(1,0) 0.0 Z/2
    dLQbcccdxwb(0,0) 2.56897060 Z/3 + Z
    L13n9331(0,0)(0,0)(0,0) Z + Z + Z
    m003(0,0) Z/5 + Z
    m004(0,0) Z
    v1205(2,3) Z/40
    x012(0,0)(0,0) Z/2 + Z
    y123(0,0) Z
    L13n9331(3,4)(2,3)(2,1) Z/53
    K7_1(0,0) Z
    6_1(0,0) Z
    5^2_1(3,4)(1,-2) Z/3
    8^3_3(0,0)(0,0)(0,0) Z + Z + Z
    4_1(0,0) Z
    12n123(0,0) Z
    16n1235(0,0) Z
    b++RL(0,0) Z
    b-+RRL(0,0) Z/3 + Z
    b+-RL(0,0) Z/5 + Z
    b--RRL(0,0) Z/3 + Z
    Braid:[1, 2, -1, -2](0,0)(0,0) Z + Z
    DT:[(8, 10, -14), (2, 6, 20), (-4, 22, 24, 12, 26, 18, 16)](0,0)(0,0)(0,0) Z + Z + Z
    DT[4,6,2](0,0) Z
    DT[mcccgdeGacjBklfmih](0,0)(0,0)(0,0) Z + Z + Z
    DT:mcccgdeGacjBklfmih(0,0)(0,0)(0,0) Z + Z + Z
    a0*B1(0,0) Z
    b1*A0 a0*B1(1,0) Z/2
    dLQbcccdxwb(0,0) Z/3 + Z
    """

    M = database.HTLinkExteriors['L13n9331']
    specs = [M._to_string(), 'm003', 'm004', 'v1205(2,3)', 'x012', 'y123',
         'L13n9331(3,4)(2,3)(2,1)', 'K7_1', '6_1',
         '5^2_1(3,4)(1,-2)', '8_3^3', 'L104001', '12n123', '16n1235',
         'b++RL', 'b-+RRL', 'b+-RL', 'b--RRL',
         'Braid[1,2,-1,-2]', 'DT:'+repr(M.DT_code()), 'DT[4,6,2]',
         'DT['+M.DT_code(alpha=True) + ']',
         'DT:'+M.DT_code(alpha=True),
         'Bundle(S_{1,1}, [a0, B1])', 'Splitting(S_{1,0}, [b1, A0], [a0,B1])',
         'dLQbcccdxwb',
         ]

    for spec in specs:
        M = Manifold(spec)
        vol = M.volume()
        if abs(vol) < 0.1:
            vol = 0.0
        print M, vol, M.homology()

    for spec in specs:
        M = Triangulation(spec)
        print M, M.homology()

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

def get_HT_knot_by_index(alternation, index, manifold_class):
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
    return manifold_class(name)

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
    path = str(manifold_path)

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
            self.path, num_tet, self.orientability, census_index)
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
    path = str(manifold_path)

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
            self.path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Triangulation(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill(( int(m),int(l)) )
        return result.with_hyperbolic_structure()

class ObsNonorientableClosedCensus(Census):
    """
    Obsolete.
    """
    data = None
    orientability = Orientability.index('nonorientable')
    path = str(manifold_path)

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
            self.path, int(num_tet), self.orientability, int(index))
        if c_triangulation == NULL:
            print(num_tet, index)
            raise RuntimeError('SnapPea failed to read the census manifold.')
        result = Triangulation(spec='empty')
        result.set_c_triangulation(c_triangulation)
        result.dehn_fill( (int(m),int(l)) )
        return result.with_hyperbolic_structure()

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
            return get_HT_knot_by_index(self.alternation, n, _manifold_class)

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
                result =  Triangulation('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result.with_hyperbolic_structure()
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
                result =  Triangulation('empty')
                result.set_c_triangulation(c_triangulation)
                result.set_name(name)
                return result.with_hyperbolic_structure()              

#----------------------------------------------------------------
#
#  Morwen's link table (data taken from Snap).  Now obsolete.
#
#----------------------------------------------------------------

class MorwenLinks:
    def __init__(self, num_components, num_crossings = None):
        raise ValueError(
            'The MorwenLinks class has been replaced by HTLinkExteriors. '
            'Please refer to the documentation for HTLinkExteriors.'
            )

def left_pad_string(str,  length, pad=' '):
    return pad*(length - len(str)) + str

class ObsMorwenLinks(Census):
    """
    Obsolete.
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
            name = str('DT:%s'%self.DT_codes[n])
            result = Manifold(name)
            result.set_name(name)
            return result

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
    name = to_byte_str('Braid:' + repr(braid_word))
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
