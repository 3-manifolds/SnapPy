# Python modules
import os, sys, operator, types, re, gzip, struct, tempfile
import tarfile, atexit, math, string, time
python_major_version = sys.version_info[0]

# Sage interaction
from snappy.sage_helper import _within_sage
from snappy.pari import pari as pari
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
except:
    pass

## SnapPy components
import spherogram
from .manifolds import __path__ as manifold_paths
from . import database
from . import twister
from . import snap
from . import verify
from . import decorated_isosig
from .ptolemy import manifoldMethods as ptolemyManifoldMethods
from .export_stl import stl
from .exceptions import SnapPeaFatalError
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

# Integer types for Python 2 vs 3.
if python_major_version < 3:
    def is_int(x):
        return isinstance(x, int) or isinstance(x, long)
else:
    def is_int(x):
        return isinstance(x, int)

# Helper function to emulate the behavior of range(n)[i]
def extract_index(i, n, formatStr):
    try:
        #
        # See https://docs.python.org/2.5/whatsnew/pep-357.html
        #
        # Note that this is not backwards compatible with python < 2.5.
        index = i.__index__()
    except:
        raise TypeError("object (%r) cannot be interpreted as index" % i)
        
    if not is_int(index):
        raise TypeError("__index__ returned non-(int, long)")
    
    if index < -n or index >= n:
        raise IndexError(formatStr % index)
    
    if index < 0:
        return index + n
    
    return index

# A stream for asynchronous messages
class MsgIO(object):
    def __init__(self):
        self.write = sys.stdout.write
        self.flush = sys.stdout.flush()

msg_stream = MsgIO()

# A very basic matrix class
class SimpleMatrix(object):
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


show_uAcknowledge = False

cdef public void uAcknowledge(const_char_ptr message):
    if show_uAcknowledge:
        sys.stderr.write(to_str(<char *> message))
        sys.stderr.write('\n')
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

cdef class AbelianGroup(object):
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

cdef class PresentationMatrix(object):
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
    
class Isometry(object):
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


# PLink communication

def _plink_callback(LE):
    cdef Manifold manifold
    cdef c_Triangulation* c_triangulation = NULL
    if LE.manifold is None:
        LE.manifold = _manifold_class('empty')
    manifold = LE.manifold
    klp = LE.SnapPea_KLPProjection()
    if klp is not None:
        manifold._set_DTcode(spherogram.DTcodec(*LE.DT_code()))
        manifold._set_PDcode(LE.PD_code())
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
