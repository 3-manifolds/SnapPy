# Python modules
import os
import sys
import operator
import types
import re
import gzip
import struct
import tempfile
import tarfile
import atexit
import math
import string
import time
import typing
python_major_version = sys.version_info.major

# Sage interaction
from .sage_helper import _within_sage, SageNotAvailable
from .pari import pari as pari
try:
    import sage.structure.sage_object
    from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.groups.free_group import FreeGroup
    from sage.interfaces.gap import gap
    from sage.interfaces.magma import magma
    try:
        from sage.interfaces.abc import GapElement, MagmaElement

        def is_GapElement(elt):
            return isinstance(elt, GapElement)

        def is_MagmaElement(elt):
            return isinstance(elt, MagmaElement)
    except:
        from sage.interfaces.gap import is_GapElement
        from sage.interfaces.magma import is_MagmaElement
    # for testing:
    from sage.matrix.constructor import matrix as sage_matrix
except ImportError:
    pass

from .matrix import matrix
from . import number

# SnapPy components
import spherogram
from .manifolds import __path__ as manifold_paths
from . import database
from . import twister
from . import snap
from . import verify
from .verify import complex_volume as verify_complex_volume
from .cusps import compute_cusp_shapes as cusps_compute_cusp_shapes
from . import decorated_isosig
from .ptolemy import manifoldMethods as ptolemyManifoldMethods
from .export_stl import stl
from .exceptions import SnapPeaFatalError
try:
    from plink import LinkEditor, LinkManager
except ImportError:
    LinkEditor, LinkManager = None, None
try:
    from .gui import ViewerWindow
except ImportError:
    ViewerWindow = None
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


def valid_index(i, n, format_str):
    """
    Return range(n)[i] or raises a nicely formatted IndexError
    using format_str.

    This does several things for us::

        * avoid that cython happily converts a float to an int, so a call
          such as M.dehn_fill((1,0), 0.6) would succeed.
        * converts a negative index such as -1 to n - 1
        * checks that i is in the right range
        * supports Sage and numpy Integers: they are not python int's but
          have __index__ (see PEP 357)

    It is faster than reimplementing these behaviors.
    """
    try:
        return range(n)[i]
    except IndexError:
        raise IndexError(format_str % i)


# A stream for asynchronous messages
class MsgIO():
    def __init__(self):
        if sys.stdout is None:
            self.write = None
            self.flush = None
        else:
            self.write = sys.stdout.write
            self.flush = sys.stdout.flush


msg_stream = MsgIO()


def to_str(s):
    return s.decode()


def bytearray_to_bytes(x):
    return bytes(x)


def to_byte_str(s):
    return s.encode('utf-8') if type(s) != bytes else s


# Types of covering spaces
cover_types = {1: "irregular", 2: "regular", 3: "cyclic"}

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
    # Only raise exception the first time so that we see the first
    # uFatalError which is usually the root cause of the problem.
    if not PyErr_Occurred():
        raise SnapPeaFatalError('SnapPea crashed in function %s(), '
                                'defined in %s.c.' % (function, file))

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
    sys.stderr.write('Q: %s\nA:  %s\n' % (<char *> message, default))
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


def set_rand_seed(seed):
    """
    Set seed of the random number generator used for choosing
    triangulation moves.
    """
    srand(seed)


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

    python_matrix = [[c_matrix.entries[i][j]
                      for j in range(c_matrix.num_cols)]
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


class MatrixWithExplanations():

    def __init__(self, mat, explain_rows, explain_columns):

        self.matrix = mat
        self.explain_rows = explain_rows
        self.explain_columns = explain_columns

    def __repr__(self, _type_str = "MatrixWithExplanations"):

        def format_explain_list(l):
            if len(l) <= 6:
                return repr(l)

            return "[ %s, ..., %s]" % (
                ', '.join(map(repr, l[:3])),
                ', '.join(map(repr, l[-3:])))

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

# C type for a function of Real returning an object
ctypedef object (*func_real_to_obj)(Real)

# Convert Real to gen in an appropriate manner for this environment
cdef func_real_to_obj Real2gen

if hasattr(pari, '_real_coerced_to_bits_prec'):  # Cypari
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
        # Pari idiosyncratically formats small and large numbers as,
        # e.g., "1.0 E-10" (note the space before "E").
        # Remove it - otherwise it cannot be parsed.
        string = string.replace(' ', '')
        float(string)
    except:
        raise ValueError('Cannot convert %s to a Real.' % type(obj))
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

cdef B2B(Boolean B):
    return B != 0

# Infos are immutable containers which hold information about SnapPy
# objects.  The base class for these is Info. Subclasses should
# override __repr__ to appropriately display the information they
# contain.


class Info(dict):
    """
    Immutable dictionary whose data can be accessed either as
    attributes or by mapping.
    Initialize with keyword arguments, or **dict.
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
            return ('Cusp %-2d: %s with Dehn filling coefficients (M, L) = %s'%
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


def _StrLongestLen(l):
    return str(max(len(e) for e in l))


LengthSpectrumFormatStringBase = (
    '%-4s '                                   # Multiplicity
    '%-40s '                                  # Length
    '%-' + _StrLongestLen(Orbifold1) + 's ')  # Topology

LengthSpectrumFormatString = (
    LengthSpectrumFormatStringBase +
    '%s')                                     # Parity

LengthSpectrumFormatStringWithWord = (
    LengthSpectrumFormatStringBase +
    '%-6s '                                   # Parity
    '%s')                                     # Word

ShortMatrixParity = { MatrixParity[0] : '-', MatrixParity[1] : '+' }


class LengthSpectrumInfo(Info):
    def __repr__(self):
        lenStr = "%17.14f" % self.length.real()
        absImag = abs(self.length.imag())
        if absImag > 1e-9:
            if self.length.imag() > 0:
                lenStr += " + "
            else:
                lenStr += " - "
            lenStr += "%16.14f*I" % absImag

        if 'word' in self:
            return LengthSpectrumFormatStringWithWord % (
                self.multiplicity,
                lenStr,
                self.topology,
                ShortMatrixParity[self.parity],
                self.word)
        else:
            return LengthSpectrumFormatString % (
                self.multiplicity,
                lenStr,
                self.topology,
                ShortMatrixParity[self.parity] )


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
        if len(self) > 0 and 'word' in self[0]:
            base = LengthSpectrumFormatStringWithWord % (
                'mult', ' length', 'topology', 'parity', 'word')
        else:
            base = LengthSpectrumFormatString % (
                'mult', ' length', 'topology', 'parity')
        return '\n'.join([base] + [repr(s) for s in self])


class ListOnePerLine(list):
    def __repr__(self):
        return '[' + ',\n '.join([repr(s) for s in self]) + ']'


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
            line0.append(l0.ljust(L))
            line1.append(l1.ljust(L))
            line2.append(l2.ljust(L))
        line3 = 'Extends to link' if self.extends_to_link() else 'Does not extend to link'
        return '\n'.join(['  '.join(l).strip()
                          for l in [line0, line1, line2]] + [line3])


cdef IsometryListToIsometries(IsometryList *isometries):
    cdef int n, c, i = 0, j = 0, c_cusp_image
    cdef MatrixInt22  c_cusp_map
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
    if klp is None:
        raise RuntimeError('Communication with PLink failed.')
    elif klp == (0, 0, 0, []):
        msg_stream.write('\nGot the empty link from PLink and so not changing the triangulation.\n')
    else:
        manifold._set_DTcode(spherogram.DTcodec(*LE.DT_code()))
        manifold._set_PDcode(LE.PD_code())
        c_triangulation = get_triangulation_from_PythonKLP(klp)
        if c_triangulation is not NULL:
            find_complete_hyperbolic_structure(c_triangulation)
        if manifold.c_triangulation is not NULL:
            free_triangulation(manifold.c_triangulation)
        manifold.set_c_triangulation(c_triangulation)
        manifold._cache.clear(message='plink_callback')
        msg_stream.write('\nNew triangulation received from PLink!\n')


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
        # Need to remove any finite vertices to be able to find a
        # hyperbolic structure.
        count_cusps(c_triangulation)
        if get_num_fake_cusps(c_triangulation) > 0:
            remove_finite_vertices(c_triangulation)
            count_cusps(c_triangulation)

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
