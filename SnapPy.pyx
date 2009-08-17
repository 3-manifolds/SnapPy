import os, sys, operator, types, re, gzip, struct, tempfile, tarfile, atexit
from signal import signal, SIGINT, SIG_DFL
from manifolds import __path__ as manifold_paths

include "SnapPy.pxi"

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
    warning = 'To do matrix algebra, please install numpy or run SnapPy in sage.'

    def __init__(self, list_of_lists):
        self.data = list_of_lists
        self.type = type(self.data[0][0])

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
        if type(key) == types.SliceType:
            raise TypeError, "Simple matrices don't slice."
        if type(key) == types.IntType:
            raise TypeError, "Simple matrices need 2 indices."
        i, j = key
        if i < 0 or j < 0:
            raise TypeError, "Simple matrices don't have negative indices." 
        return key

    def __getitem__(self, key):
        i, j = self._check_indices(key)
        return self.data[i][j]

    def __setitem__(self, key, value):
        i, j = self._check_indices(key)
        self.data[i][j] = value

    def __add__(self, other):
        raise TypeError, self.warning

    def __sub__(self, other):
        raise TypeError, self.warning

    def __mul__(self, other):
        raise TypeError, self.warning

    def __div__(self, other):
        raise TypeError, self.warning

    def __inv__(self, other):
        raise TypeError, self.warning

try:
    from sage.matrix.constructor import matrix
except ImportError:
    try:
        from numpy import matrix
    except ImportError:
        matrix = SimpleMatrix

# Sage interaction
try:
    import sage.structure.sage_object
    from sage.groups.perm_gps.permgroup_element import is_PermutationGroupElement
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from sage.interfaces.gap import gap
    from sage.interfaces.gap import is_GapElement
    from sage.interfaces.magma import magma
    from sage.interfaces.magma import is_MagmaElement
    
    _within_sage = True
except ImportError:
    _within_sage = False

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

# Paths
manifold_path = manifold_paths[0] + os.sep
closed_census_directory = os.path.join(manifold_path, 'ClosedCensusData')
link_directory = os.path.join(manifold_path, 'ChristyLinks')
link_archive = os.path.join(manifold_path, 'ChristyLinks.tgz')
census_knot_archive = os.path.join(manifold_path, 'CensusKnots.tgz')
table_directory = os.path.join(manifold_path, 'HTWKnots')

# These are the gzipped files holding the knot tables.
Alternating_table = gzip.open(os.path.join(table_directory, 'alternating.gz') )
Nonalternating_table = gzip.open(os.path.join(table_directory, 'nonalternating.gz') )

# This is the gzipped tarball of Joe Christy's link complements.
Christy_links = tarfile.open(link_archive, 'r:*')

# This is the gzipped tarball of the knots in the SnapPea census,
# as classified by Callahan, Dean, Weeks, Champanerkar, Kofman, Patterson.
Census_Knots = tarfile.open(census_knot_archive, 'r:*')

# Implementation of the SnapPea UI functions and their global variables.
cdef extern from *:
    ctypedef char* const_char_ptr "const char*"
    ctypedef int const_int "const int"

class SnapPeaFatalError(Exception):
    """
    This exception is raised by SnapPy when the SnapPea kernel encounters
    a fatal error.
    """

cdef public void uFatalError(char *function, char *file) except *:
    raise SnapPeaFatalError,  """
SnapPea crashed in function %s(), defined in %s.c."""%(function, file)

cdef public Boolean gLongComputationInProgress
cdef public Boolean gLongComputationCancelled
    
def SnapPea_handler(signal, stackframe):
    """
    A Python signal handler which cancels the SnapPea computation.
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
    # Set SnapPea's flags
    gLongComputationCancelled = False
    gLongComputationInProgress = True
    # Install the SnapPea handler on SIGINT and SIGALRM
    signal(SIGINT, SnapPea_handler)
#    signal(SIGALRM, SnapPea_handler)

cdef public c_FuncResult uLongComputationContinues():
    global gLongComputationCancelled
    global gLongComputationInProgress
    # While a SnapPea function is executing, Python saves all of its
    # calls to interrupt handlers on a list of "Pending Calls".
    # This forces any waiting handlers to be called.
    Py_MakePendingCalls()
    if gLongComputationCancelled:
        gLongComputationCancelled = False
        gLongComputationInProgress = False
        return func_cancelled
    else:
        return func_OK

cdef public void uLongComputationEnds():
    global gLongComputationCancelled
    global gLongComputationInProgress
    # Restore Python's default signal handler on SIGINT and SIGALRM
    signal(SIGINT, python_handler)
#    signal(SIGALRM, python_handler)
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

# PARI support for Smith normal form

# We do this to keep PARI from stealing our keyboard interrupts.
python_handler = signal(SIGINT, SIG_DFL)
pari_init_opts(1000000,500000,0)
signal(SIGINT, python_handler)

def smith_form(M):
    cdef GEN pari_matrix
    cdef GEN pari_vector
    cdef GEN pari_int
    cdef int i, j
    try:
        m, n = M.shape
    except AttributeError:
        # probably means we're within Sage
        m, n = M.nrows(), M.ncols()
        
    pari_matrix = cgetg(n+1, t_MAT)
    for j from 1 <= j <= n:
        pari_matrix[j] = <long>cgetg(m+1, t_COL)
    for i from 1 <= i <= m:
        for j from 1 <= j <= n:
            (<GEN*>pari_matrix)[j][i] =  <long>stoi(M[i-1,j-1])
    pari_vector = matsnf0(pari_matrix, 4)
    result = []
    for i from 1 <= i < lg(pari_vector):
        pari_int = (<GEN*>pari_vector)[i]
        result.append(itos(pari_int))
    if m < n:
        result = result + [0]*(n-m)

    cgiv(pari_vector)
    cgiv(pari_matrix)
    # PARI views the input to matsnf0 as square.
    if m > n:
        for i in range(m - n):
            result.remove(0)

    # For consistency with SnapPea, need to switch the order of the factors.
    zeros = [x for x in result if x == 0]
    nonzeros = [x for x in result if x != 0]
    nonzeros.sort()
    return nonzeros + zeros

# Enum conversions
CuspTopology = ['torus cusp', 'Klein bottle cusp', 'unknown']
MatrixParity = ['orientation-reversing', 'orientation-preserving']
Orientability = ['orientable', 'nonorientable', 'unknown']
Orbifold1 = ['unknown', 'circle', 'mirrored arc']
FuncResult = ['func_OK', 'func_cancelled', 'func_failed', 'func_bad_input']
SolutionType = ["not attempted", "all tetrahedra positively oriented",
                "contains negatively oriented tetrahedra", "contains flat tetrahedra",
                "contains degenerate tetrahedra", "unrecognized solution type",
                "no solution found"]

# global functions
def check_SnapPea_memory():
    verify_my_malloc_usage()

# SnapPea Classes

# Abelian Groups

cdef class AbelianGroup:
    """
    An AbelianGroup object represents a finitely generated abelian group,
    usually the first homology group of a snappy Manifold.

    Instantiate as AbelianGroup([n_1, n_2, ... ]) where the n_i are the
    orders of the cyclic factors (or 0, in the case of an infinite cyclic
    factor).

    >>> A = AbelianGroup([5,15,0,0])
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

    cdef readonly coefficients

    def __init__(self, coefficients):
        try:
            self.coefficients = list(coefficients)
        except:
            raise RuntimeError, "Argument is not a sequence."
        for c in self.coefficients:
            assert type(c) == types.IntType and c >= 0,\
                'Coefficients must be non-negative integers.\n'
        for i in range(len(coefficients) - 1):
            n,m = coefficients[i:i+2]
            assert (n == m == 0) or (m % n == 0), 'Each coefficient must divide the subsequent one, but was given %s' % coefficients

    def __repr__(self):
        if len(self.coefficients) == 0:
            return '0'
        factors = ( ['Z/%d'%n for n in self.coefficients if n > 1] +  ['Z' for n in self.coefficients if n == 0])
        return ' + '.join(factors)

    def __len__(self):
        return len(self.coefficients)
    
    def __getitem__(self, i):
        return self.coefficients[i]

    def __cmp__(self, other):
        return cmp(self.coefficients, other.coefficients)
            
    def rank(self):
        """
        The rank of the group.
        """
        return len(self.coefficients)

    def betti_number(self):
        """
        The rank of the maximal free abelian subgroup.
        """
        return len([n for n in self.coefficients if n == 0])
    
    def order(self):
        """
        The order of the group.  Returns the string 'infinite' if the
        group is infinite.        
        """
        det = reduce(operator.mul, [1] + self.coefficients)
        return 'infinite' if det == 0 else det

# Helper class for cusp info

class CuspInfoDict(dict):
    def __repr__(self):
        if self['complete?']:
            if self.has_key('shape'):
                return 'Cusp %-2d: complete %s of shape %s' % \
                      (self['index'], self['topology'], self['shape'])
            return 'Cusp %-2d: %s, not filled'% (self['index'], self['topology'])
        else:
            return 'Cusp %-2d: %s with Dehn filling coeffients (M, L) = %s'%\
                   (self['index'], self['topology'], self['filling'])

class DualCurveDict(dict):
    def __repr__(self):
        return '%3d: %s curve of length %s'% \
               (self['index'],MatrixParity[self['parity']],
                self['filled length'])

class ListOnePerLine(list):
    def __repr__(self):
        return "[" + ",\n ".join([repr(s) for s in self]) + "]"

# Isometry


def format_two_by_two(mat):
    a,b,c,d = ["%d" % x for x in [mat[0,0], mat[0,1], mat[1,0], mat[1,1]]]
    w0 = max(len(a), len(c))
    w1 = max(len(b), len(d))
    return "[" + a.rjust(w0) + " " + b.rjust(w1) + "]", "[" + c.rjust(w0) + " " + d.rjust(w1) + "]"
    
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
            l0 = "%d -> %d" % (i, images[i])
            l1, l2 = format_two_by_two(maps[i])
            L = max(len(l0), len(l1), len(l2))
            line0.append( l0.ljust(L) )
            line1.append( l1.ljust(L) )
            line2.append( l2.ljust(L) )
        line3 = "Extends to link" if self.extends_to_link() else "Does not extend to link"
        return "\n".join(["  ".join(line0), "  ".join(line1), "  ".join(line2), line3])

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
        
    

# Triangulations

cdef class Triangulation:
    """
    A Triangulation object represents a compact 3-manifold with
    boundary a union of tori by an ideal triangulation of the
    manifold's interior.  A Dehn-filling can be specified for each
    boundary component, allowing the description of closed 3-manifolds
    and some orbifolds.  For non-orientable 3-manifolds, the boundary
    components can also be Klein bottles. Two Triangulations are equal
    ('==') if they represent combinatorially isomorphic
    triangulations.  A Triangulation does *not* have any geometric
    structure, and usually one works with the subclass Manifold which
    adds this.

    A Triangulation can be specified in a number of ways, e.g.

    - Triangulation('9_42') : The complement of the knot 9_42 in S^3.
    - Triangulation('125(1,2)(4,5)') : The SnapPea census manifold m125
       where the first cusp has Dehn filling (1,2) and the second cusp has
       filling (4,5).
    - Triangulation() : Opens a link editor window where can you
       specify a link complement.
    
    In general, the specification can be from among the below, with
    information on Dehn fillings added.

    - SnapPea cusped census manifolds: e.g. 'm123', 's123', 'v123'.

    - Link complements:
       + Rolfsen's table: e.g. '4_1', '04_1', '5^2_6', '6_4^7', 'L20935', 'l104001'.
       + Hoste-Thistlethwaite Knotscape table:  e.g. '11a17' or '12n345'
       + Dowker-Thistlethwaite code: e.g. 'DT[6,8,2,4]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.
      
    """

    cdef c_Triangulation* c_triangulation
    cdef readonly _cache
    cdef readonly LE

    def __new__(self, spec=None):
        cdef c_Triangulation *c_triangulation = NULL
        # Answers to potentially hard computations are cached
        self._cache = {}
        self.LE = None
        if spec is not None and spec != 'empty':
            if type(spec) != types.StringType:
                raise TypeError, triangulation_help%self.__class__.__name__
            c_triangulation = get_triangulation(spec)
            if c_triangulation == NULL:
                raise RuntimeError, "An empty triangulation was generated."
        if spec is None:
            if LinkEditor:
                print 'Starting the link editor.\n'\
                      'Select PLink->Send to SnapPy to load the link complement.'
                self.LE = LinkEditor(no_arcs=True,
                                     callback=self._plink_callback,
                                     cb_menu='Send to SnapPy')
            else:
                raise RuntimeError, "PLink was not imported."

        if c_triangulation != NULL:    
            self.set_c_triangulation(c_triangulation)
            remove_hyperbolic_structures(c_triangulation)

    def _plink_callback(self):
        cdef c_Triangulation* c_triangulation = NULL
        if self.LE is not None:
            klp = self.LE.SnapPea_KLPProjection() 
            if klp is not None:
                c_triangulation = get_triangulation_from_PythonKLP(klp)
                self.set_c_triangulation(c_triangulation)
                msg_stream.write('\nNew triangulation received from PLink!\n')
                return
        else:
            raise RuntimeError, "Communication with PLink failed." 

    def plink(self):
        """
        Brings the corresponding link editor window to the front, if
        there is one.
        """
        if self.LE is not None:
            self.LE.reopen()
        else:
            raise ValueError, "This manifold does not have a PLink window."

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
            raise ValueError, "Acceptable cusp types are ['all', 'orientable', 'nonorientable']"

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
            raise ValueError, "Triangulation is already orientable"

        cdef c_Triangulation* cover_c_triangulation = NULL
        cdef Triangulation new_tri
        
        cover_c_triangulation = double_cover(self.c_triangulation)
        new_tri = Triangulation('empty')
        new_tri.set_c_triangulation(cover_c_triangulation)
        new_tri.set_name(self.name() + "~")
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

    def with_hyperbolic_structure(self):
        """
        Add a (possibly degenerate) hyperbolic structure, turning the
        Triangulation into a Manifold.

        >>> M = Triangulation('m004')
        >>> N = M.with_hyperbolic_structure()
        >>> N.volume()
        2.0298832128193069
        """
        return Manifold_from_Triangulation(self)

    def save(self, file_name):
        """
        Save the triangulation as a SnapPea triangulation file.

        >>> M = Triangulation('m004')
        >>> M.save('fig-eight.tri')
        """
        if self.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'
        write_triangulation(self.c_triangulation, file_name)

    def __dealloc__(self):
        if self.c_triangulation is not NULL:
            free_triangulation(self.c_triangulation)

    def __richcmp__(Triangulation self, Triangulation other, case):
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
        ValueError: Can't compare triangulations of manifolds with Dehn fillings
        """
        cdef c_Triangulation *c_triangulation1
        cdef c_Triangulation *c_triangulation2
        cdef Boolean answer
        if case != 2:
            return NotImplemented
        if type(self) != type(other):
            return False
        if False in self.cusp_info('complete?') + other.cusp_info('complete?'):
            raise ValueError, "Can't compare triangulations of manifolds with Dehn fillings"
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
                if info['complete?']:
                    repr += '(0,0)'
                else:
                    repr += '(%g,%g)'% info['filling']
            return repr
 
    def set_name(self, new_name):
        """
        Give the triangulation a new name.

        >>> M = Triangulation('4_1')
        >>> M.set_name('figure-eight-comp')
        >>> M
        figure-eight-comp(0,0)
        """
        cdef char* c_new_name = new_name
        if self.c_triangulation is NULL:
            raise ValueError, 'The empty triangulation has no name.'
        set_triangulation_name(self.c_triangulation, c_new_name)

    def name(self):
        """
        Return the name of the triangulation.

        >>> M = Triangulation('4_1')
        >>> M.name()
        'L104001'
        """
        if self.c_triangulation is NULL: return
        return get_triangulation_name(self.c_triangulation)

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
          L408001(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          L408001(0,0)(0,0)(2,3)(1,5)
        
        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          L408001(3,0)(1,-4)(2,3)(1,5)

        - When there is only one cusp, there's a shortcut

          >>> N = Triangulation('m004')
          >>> N.dehn_fill( (-3,4) )
          >>> N
          m004(-3,4)
        
        Does not return a new Triangulation.
        """
        if self.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty'

        if which_cusp != None:
            try:
                which_cusp = range(self.num_cusps())[which_cusp]
            except IndexError:
                raise IndexError, 'The specified cusp (%s) does not exist'%which_cusp

            meridian, longitude = filling_data
            complete = ( meridian == 0 and longitude == 0)
            set_cusp_info(self.c_triangulation,
                          which_cusp, complete, meridian, longitude)
            self._cache = {}
        else:
            if self.num_cusps() == 1 and len(filling_data) == 2:
                self.dehn_fill(filling_data, 0)
                return 
            if len(filling_data) > self.num_cusps():
                raise IndexError, 'Provided more filling data than there are cusps.'
            for i, fill in enumerate(filling_data):
                self.dehn_fill(fill, i)
                    

    def cusp_info(self, data_spec=None):
        """
        Returns a dictionary containing information about the given
        cusp.   Usage:

        >>> M = Triangulation('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0)
        >>> c = M.cusp_info(1)
        >>> c['complete?']
        False
        >>> c.keys()
        ['index', 'filling', 'topology', 'complete?']

        You can get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : torus cusp, not filled,
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('complete?')
        [True, False, False]
        """
        cdef int cusp_index
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef double m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_precision, current_shape_precision,
        cdef Complex initial_modulus, current_modulus
        cdef int meridian_precision, longitude_precision
        cdef Complex c_meridian, c_longitude

        if self.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'

        if data_spec == None:
            return ListOnePerLine([self.cusp_info(i) for i in range(self.num_cusps())])

        if type(data_spec) == type(''):
            return [c[data_spec] for c in self.cusp_info()]

        try:
            cusp_index = range(self.num_cusps())[data_spec]
        except IndexError:
            raise IndexError, 'The specified cusp (%s) does not exist'%data_spec

        get_cusp_info(self.c_triangulation, cusp_index,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_precision, &current_shape_precision,
                      &initial_modulus, &current_modulus)
        ans = {'index' : cusp_index,
               'topology' : CuspTopology[topology],
               'complete?' : B2B(is_complete),
               'filling' : (m, l)}

        #If there's a hyperbolic structure, there more information to
        #pass on.
        if hasattr(self, 'tetrahedra_shapes'):
            get_holonomy(self.c_triangulation, cusp_index,
                         &c_meridian, &c_longitude,
                         &meridian_precision, &longitude_precision)


            ans = CuspInfoDict({'index' : cusp_index,
                                'topology' : CuspTopology[topology],
                                'complete?' : B2B(is_complete),
                                'filling' : (m, l),
                                'shape' : C2C(current_shape),
                                'shape precision' : current_shape_precision,
                                'modulus' : C2C(current_modulus),
                                'holonomies' : (C2C(c_meridian), C2C(c_longitude)),
                                'holonomy precision' : min(meridian_precision, longitude_precision)
                                })

        return CuspInfoDict(ans)

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
            raise ValueError, "Manifold not orientable so can't reverse orientation."
        reorient(self.c_triangulation)
        self._cache = {}
            
    def filled_triangulation(self, cusps_to_fill="all"):
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
            raise ValueError, 'Triangulation is empty.'
        n = self.num_cusps()
        if cusps_to_fill == "all":
            cusps_to_fill = [c for c in range(n) if cusp_is_fillable(self.c_triangulation, c)]
                
        if False in [(c in range(n)) for c in cusps_to_fill]:
            raise IndexError, "Specified indices to be filled are beyond the actual number of cusps"
        if 0 in [cusp_is_fillable(self.c_triangulation, c) for c in cusps_to_fill]:
            raise IndexError, "To permanently fill a cusp, the Dehn filling coefficients must be relatively prime integers."

        cdef c_Triangulation* c_filled_tri = NULL
        cdef Triangulation filled_tri
        cdef Boolean *fill_cusp_spec = NULL
        
        copy_triangulation(self.c_triangulation, &c_filled_tri)

        fill_cusp_spec = <Boolean*>malloc(n*sizeof(Boolean))
        for i in range(n):
            fill_cusp_spec[i] = 1 if i in cusps_to_fill else 0

        fill_all = 1 if not False in [i in cusps_to_fill for i in range(n)] else 0

        
        c_filled_tri = fill_cusps(self.c_triangulation, fill_cusp_spec, "", fill_all)

        free(fill_cusp_spec)

        filled_tri = Triangulation('empty')
        filled_tri.set_c_triangulation(c_filled_tri)
        filled_tri.set_name(self.name() + "_filled")

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
            raise ValueError, 'Triangulation is empty.'
        while get_num_edge_classes(self.c_triangulation, v, 1) > 0:
            c = get_num_edge_classes(self.c_triangulation, v, 0)
            if c > 0:
                ans[v] = c
            v += 1
        return ans
        

    def gluing_equations(self,form='log'):
        """
        In the default mode, this function returns a matrix with rows of the form

                  a b c  d e f  ...

        which means

            a*log(z0) + b*log(1/(1-z0)) + c*log((z0-1)/z0) + d*log(z1) +... = 2 pi i

        for an edge equation, and (same) = 1 for a cusp equation.
        Here, the cusp equations come at the bottom of the matrix.  

        In terms of the tetrahedra, a is the invariant of the edge
        (2,3), b the invariant of the edge (0,2) and c is the
        invariant of the edge (1,2).  See kernel_code/edge_classes.c
        for a detailed account of the convention used.  

        If the optional argument form='rect' is given, then this
        function returns a list of tuples of the form:

           ( [a0, a1,..,a_n], [b_0, b_1,...,b_n], c)

        where this corresponds to the equation

           z0^a0 (1 - z0)^b0 + z1^a1(1 - z1)^b1 + ...  = c

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
            raise ValueError, 'Triangulation is empty.'
        c_eqns = get_gluing_equations(self.c_triangulation, &num_rows, &num_cols)
        eqns = [ [c_eqns[i][j] for j in range(num_cols)] for i in range(num_rows)]
        free_gluing_equations(c_eqns, num_rows)

        for i in range(self.num_cusps()):
            cusp_info = self.cusp_info(i)
            to_do = [(1,0), (0,1)] if cusp_info["complete?"] else [cusp_info["filling"]]
            for (m, l) in to_do:
                eqn = get_cusp_equation(self.c_triangulation, i, m, l, &num_rows)
                eqns.append([eqn[j] for j in range(num_rows)])
                free_cusp_equation(eqn)

        if form == "log":
            return matrix(eqns)

        if form != "rect":
            raise ValueError, "Equations available in 'log' and 'rect' forms only."

        ans = []
        for row in eqns:
            n = self.num_tetrahedra()
            a, b = [0,]*n, [0,]*n
            c = 1
            for j in range(n):
                r = row[3*j + 2]
                a[j] = row[3*j] - r
                b[j] = -row[3*j + 1] + r
                c *= -1 if r % 2 else 1
            ans.append( (a, b, c) )
        return ans
                                                     
    def homology(self):
        """
        Returns an AbelianGroup representing the first integral
        homology group of the underlying (Dehn filled) manifold.

        >>> M = Triangulation('m003')
        >>> M.homology()
        Z/5 + Z
        """
        if "homology" in self._cache.keys():
            return self._cache["homology"]
        
        cdef c_AbelianGroup *H
        cdef RelationMatrix R
        cdef int m, n

        if self.c_triangulation is NULL:
            return AbelianGroup([])
        coefficient_list = []
        H = homology(self.c_triangulation)
        if H != NULL:
            compress_abelian_group(H)
            for n from 0 <= n < H.num_torsion_coefficients:
                coefficient_list.append(H.torsion_coefficients[n])
            free_abelian_group(H)
        else:
            homology_presentation(self.c_triangulation, &R)
            relations = []
            if R.relations != NULL:
                for m from 0 <= m < R.num_rows:
                    row = []
                    for n from 0 <= n < R.num_columns:
                        row.append(R.relations[m][n])
                    relations.append(row)
                coefficient_list = smith_form(matrix(relations))
                free_relations(&R)

        self._cache["homology"] = AbelianGroup(coefficient_list)
        return self._cache["homology"]
    

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
            raise ValueError, 'Triangulation is empty.'
        name_mangled = "fundamental_group-%s-%s-%s" %\
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
            raise ValueError, 'Triangulation is empty.'
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
                if input_type == "GrpFP":
                    GG = magma(self.fundamental_group())
                    f = GG.CosetAction(permutation_rep)
                elif input_type == "HomGrp":
                    f = permutation_rep
                    if not repr(f.Image().Type()) == "GrpPerm":
                        raise TypeError, "Homorphism image is not a permutation group."
                else:
                    raise TypeError, "Magma type not recogonized."
                
                magma.eval("""\
                     FormatHomForSnapPea := function(f)
                         subone := function(L)   return [x - 1 : x in L]; end function;
                         return [subone(Eltseq(f(g))) : g in Generators(Domain(f))];
                       end function;""")
                permutation_rep = f.FormatHomForSnapPea().sage()

            # Not a useful GAP or MAGMA object, so let's try.  
            elif not False in [is_PermutationGroupElement(p) for p in permutation_rep]:
                permutation_rep = [[x - 1 for x in perm.list()] for perm in permutation_rep]

        G = self.fundamental_group()
        c_representation = self.build_rep_into_Sn(permutation_rep)
        degree = len(permutation_rep[0])
        c_triangulation = construct_cover(self.c_triangulation,
                                          c_representation,
                                          degree)
        cover = Triangulation('empty')
        cover.set_c_triangulation(c_triangulation)
        cover.set_name(self.name()+'~')
        free_representation(c_representation,
                            G.num_orig_gens(),
                            self.num_cusps())
        return cover

    def covers(self, degree, method=None):
        """
        Returns a list of Triangulations corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Triangulation('m003')
        >>> covers = M.covers(4)
        >>> [(N, N.homology()) for N in covers]
        [(m003~0(0,0)(0,0), Z/5 + Z + Z), (m003~1(0,0), Z/3 + Z/15 + Z)]

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
            raise ValueError, 'Triangulation is empty.'
        if method:
            if not _within_sage:
                raise RuntimeError, "Only the default method of finding subgroups is available, as you are not using Sage"
            if method == "gap":
                G = gap(self.fundamental_group())
                return [self.cover(H) for H in G.LowIndexSubgroupsFpGroup(degree) if G.Index(H) == degree]
            if method == "magma":
                G = magma(self.fundamental_group())
                return [self.cover(H) for H in G.LowIndexSubgroups("<%d, %d>" % (degree, degree))]

        
        reps = find_representations(self.c_triangulation,
                                        degree,
                                        permutation_subgroup_Sn)
        covers = []
        rep = reps.list
        while rep != NULL:
            cover = construct_cover(self.c_triangulation,
                                    rep,
                                    reps.num_sheets)
            T = Triangulation('empty')
            T.set_c_triangulation(cover)
            covers.append(T)
            rep = rep.next
        free_representation_list(reps)
        for i in range(len(covers)):
            covers[i].set_name(self.name() + '~%d'%i)
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
                raise ValueError, "Not a valid permutation list."

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
            message = """
        Invalid permutation data."""
            failed = True
        if c_repn_in_original_gens == NULL:
            message = """
        Failed to construct permutation representation."""
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
            raise RuntimeError, message
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
    and some orbifolds.

    A Manifold can be specified in a number of ways, e.g.

    - Manifold('9_42') : The complement of the knot 9_42 in S^3.
    - Manifold('125(1,2)(4,5)') : The SnapPea census manifold m125
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
       + Callahan-Dean-Weeks-Champanerkar-Kofman-Patterson knots: e.g. 'K_6_21'.
       + Dowker-Thistlethwaite code: e.g. 'DT[6,8,2,4]'

    - Once-punctured torus bundles: e.g. 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'

    - Fibered manifold associated to a braid: 'braid[1,2,-3,4]'
    
      Here, the braid is thought of as a mapping class of the
      punctured disc, and this manifold is the corresponding
      mapping torus.  If you want the braid closure, do (1,0) filling
      of the last cusp.

    - A SnapPea triangulation or link projection file: 'filename'

      The file will be loaded if found in the current directory or the
      path given by the shell variable SNAPPEA_MANIFOLD_DIRECTORY.
      


    """

    def __init__(self, spec=None):
        if self.c_triangulation != NULL:
            find_complete_hyperbolic_structure(self.c_triangulation)
            do_Dehn_filling(self.c_triangulation)

    
    def canonize(self):
        """
        Change the triangulation to the canonical retriangulation of
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
            raise RuntimeError, "SnapPea failed to find the canonical triangulation"
        

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

    def filled_triangulation(self, cusps_to_fill="all"):
        """
        Return a new manifold where the specified cusps have been
        permanently filled in.  Examples:

        Filling all the cusps, wich this results in a Tiangulation rather
        than a manifold, since SnapPea can't deal with hyperbolic
        structures in that case.  
        
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
        >>> G.SL2C("baaBA")
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
            raise ValueError, 'Triangulation is empty.'
        name_mangled = "fundamental_group-%s-%s-%s" %\
                       (simplify_presentation,
                        fillings_may_affect_generators,
                        minimize_number_of_generators)
        if not name_mangled in self._cache.keys():
            self._cache[name_mangled] = HolonomyGroup(self, simplify_presentation, fillings_may_affect_generators, minimize_number_of_generators)
        return self._cache[name_mangled]

    def symmetry_group(self, of_link=False):
        """
        Returns the symmetry group of the Manifold.
        If the flag"of_link" is set, then it only returns symmetries that preserves the meridians.
        """

        cdef c_SymmetryGroup* symmetries_of_manifold = NULL
        cdef c_SymmetryGroup* symmetries_of_link = NULL
        cdef c_Triangulation* c_symmetric_triangulation = NULL
        cdef Manifold symmetric_triangulation
        cdef Boolean is_full_group
        cdef c_FuncResult result
        cdef SymmetryGroup symmetry_group

        if self.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'

        name_mangled = 'symmetry_group-%s' % of_link
        if not name_mangled in self._cache.keys():
            result = compute_symmetry_group(self.c_triangulation, &symmetries_of_manifold,
                                            &symmetries_of_link, &c_symmetric_triangulation, &is_full_group)

            if result != func_OK:
                raise ValueError, "SnapPea failed to compute any part of the symmetry group."

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
        if not self._cache.has_key('symmetric_triangulation'):
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

    def covers(self, degree, method=None):
        """
        M.covers(degree, method=None)

        Returns a list of Manifolds corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        >>> M = Manifold('m003')
        >>> covers = M.covers(4)
        >>> [(N, N.homology()) for N in covers]
        [(m003~0(0,0)(0,0), Z/5 + Z + Z), (m003~1(0,0), Z/3 + Z/15 + Z)]
        
        If you are using Sage, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional
        argument method = 'gap' If you have Magma installed, you can
        used it to do the heavy lifting by specifying method =
        'magma'.
        """
        covers = Triangulation.covers(self, degree, method)
        return [Manifold_from_Triangulation(cover, False) for cover in covers]
    
    def volume(self, accuracy=False):
        """
        Returns the volume of the manifold.

        >>> M = Manifold('m004')
        >>> M.volume()
        2.0298832128193069

        If the flag accuracy is set to True, then it returns the
        volume of the manifold together with the number of digits of
        accuracy as *estimated* by SnapPea.

        >>> M.volume(True)
        (2.0298832128193069, 10)
        """
        cdef int acc
        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ("not attempted", "no solution found"):
            raise ValueError, 'Solution type is: %s'%solution_type
        vol = volume(self.c_triangulation, &acc)
        if accuracy:
            return (vol, acc)
        else:
            return vol

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

    def tetrahedra_shapes(self, part=None, fixed_alignment=True):
        """
        Gives the shapes of the tetrahedra in the current solution to
        the gluing equations.  Returns a list containing one dictionary
        for each tetrahedron.  The dictionary keys are:

        - rect : the shape of the tetrahedron, as a point in the
          complex plane.

        - log : the log of the shape

        - precision: a list of the approximate precisions of the
          shapes, in order (rect re, rect im, log re, log im)

        If the optional variable 'part'is set to one of the above,
        then the function returns only that component of the data.
        
        If the flag 'fixed_alignment' is set to False, then the edges
        used to report the shape parameters are choosen so as to
        normalize the triangle.

        >>> M = Manifold('m015')
        >>> M.tetrahedra_shapes(part='rect')
        [(0.66235897862237314+0.56227951206230109j), (0.66235897862237303+0.56227951206230109j), (0.66235897862237292+0.56227951206230098j)]
        >>> M.tetrahedra_shapes()
        [{'log': (-0.14059978716148094+0.70385772130147628j), 'rect': (0.66235897862237314+0.56227951206230109j), 'precisions': (11, 11, 12, 11)},
         {'log': (-0.14059978716148103+0.70385772130147639j), 'rect': (0.66235897862237303+0.56227951206230109j), 'precisions': (11, 11, 11, 11)},
         {'log': (-0.14059978716148125+0.70385772130147639j), 'rect': (0.66235897862237292+0.56227951206230098j), 'precisions': (11, 11, 11, 11)}]
        """        
        cdef double rect_re, rect_im, log_re, log_im
        cdef int prec_rec_re, prec_rec_im, prec_log_re, prec_log_im
        cdef Boolean is_geometric
        
        if self.c_triangulation is NULL: return []
        ans = []
        for i in range(self.num_tetrahedra()):
            get_tet_shape(self.c_triangulation, i,  fixed_alignment,
                          &rect_re, &rect_im, &log_re, &log_im,
                          &prec_rec_re, &prec_rec_im, &prec_log_re, &prec_log_im,
                          &is_geometric)
            ans.append({"rect":(rect_re + rect_im*(1J)),"log":(log_re + log_im*(1J)),
                        "precisions":(prec_rec_re, prec_rec_im, prec_log_re, prec_log_im)})

        if part != None:
            if part not in ["rect", "log", "precisions"]:
                raise ValueError, "A non-existent shape data type was specified."
            return [a[part] for a in ans]

        return ans if part else ListOnePerLine(ans)


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
            raise ValueError, 'Triangulation is empty.'
        N = get_num_tetrahedra(self.c_triangulation)
        shape_array = <Complex *>malloc(N*sizeof(Complex))
        set_cusps(self.c_triangulation, fillings)
        for i from 0 <= i < N:
            shape = complex(shapes[i]) 
            shape_array[i].real = shape.real
            shape_array[i].imag = shape.imag
        set_tet_shapes(self.c_triangulation, shape_array)
        free(shape_array)

    def solution_type(self):
        """
        Returns the type of the current solution to the gluing
        equations, basically a summary of how degenerate the solution
        is.  The possible answers are:

        - 'not attempted'
        
        - 'all tetrahedra positively oriented' aka 'geometric_solution'
          Should correspond to a genuine hyperbolic structure

        - 'contains negatively oriented tetrahedra' aka 'nongeometric_solution'
          Probably correponds to a hyperbolic structure but some
          simplices have reversed orientiations.  
             
        - 'contains flat tetrahedra' Contains some tetrahedra with
          shapes in R - {0, 1}.

        - 'contains degenerate tetrahedra' Some shapes are close to
          {0,1, or infinity}.  
        
        - 'unrecognized solution type'
        
        - 'no solution found'

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
        solution_type = get_filled_solution_type(self.c_triangulation)

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
        Returns a dictionary containing information about the given
        cusp.   Usage:

        >>> M = Manifold('v3227(0,0)(1,2)(3,2)')
        >>> M.cusp_info(1)
        Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0)

        To get more detailed information about the cusp, we do

        >>> c = M.cusp_info(0)
        >>> c['shape']
        (0.11044501762139303+0.94677097849790615j)
        >>> c['modulus']
        (-0.12155871955249957+1.0420412829322609j)
        >>> c.keys()
        ['index', 'holonomies', 'shape', 'complete?', 'filling', 'shape precision', 'holonomy precision', 'modulus', 'topology']

        Here 'shape' is the shape of the cusp,
        i.e. (longitude/meridian) and 'modulus' is its shape in the
        geometrically preferred basis, i.e.  ( (second shortest
        translation)/(shortest translation)).  For cusps that are
        filled, one instead cares about the holonomies:
        
        >>> M.cusp_info(-1)['holonomies']
        ((-0.59883088859413069+1.0981254817102275j), (0.89824633289119604+1.494404431024452j))
        
        You can also get information about multiple cusps at once:

        >>> M.cusp_info()
        [Cusp 0 : complete torus cusp of shape (0.110445017621+0.946770978498j),
         Cusp 1 : torus cusp with Dehn filling coeffients (M, L) = (1.0, 2.0),
         Cusp 2 : torus cusp with Dehn filling coeffients (M, L) = (3.0, 2.0)]
        >>> M.cusp_info('complete?')
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
          L408001(0,0)(0,0)(2,3)(0,0)

        - Fill the last cusp:

          >>> M.dehn_fill((1,5), -1)
          >>> M
          L408001(0,0)(0,0)(2,3)(1,5)
        
        - Fill the first two cusps:

          >>> M.dehn_fill( [ (3,0), (1, -4) ])
          >>> M
          L408001(3,0)(1,-4)(2,3)(1,5)

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

    def set_peripheral_curves(self, peripheral_data, which_cusp=None):
        """
        Each cusp has a preferred marking. In the case of torus
        boundary, this is pair of essential simple curves meeting in
        one point; equivalently, a basis of the first homology. These
        curves are called the meridian and the longitude.

        This method changes these markings.

        - Make the shortest curves the meridians, and the second
          shortest curves the longitudes.  

          >>> M = Manifold('5_2')
          >>> M.cusp_info('shape')
          [(-2.4902446675066177+2.9794470664789769j)]
          >>> M.set_peripheral_curves('shortest')
          >>> M.cusp_info('shape')
          [(-0.49024466750661766+2.9794470664789769j)]
          
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
            raise ValueError, 'Triangulation is empty'

        if peripheral_data == 'shortest':
            if which_cusp != None:
                raise ValueError, "Must apply 'shortest' to all the cusps"
            install_shortest_bases(self.c_triangulation)

        elif peripheral_data == 'fillings':
            if which_cusp != None:
                raise ValueError, "Must apply 'fillings' to all the cusps"
            install_current_curve_bases(self.c_triangulation)
            return
        
        elif which_cusp != None:
            try:
                which_cusp = range(self.num_cusps())[which_cusp]
            except IndexError:
                raise IndexError, 'The specified cusp (%s) does not exist'%which_cusp

            meridian, longitude = peripheral_data
            a, b = meridian
            c, d = longitude
            if a*d - b*c != 1:
                raise ValueError, 'Data provided does not give a (pos. oriented) basis'


            matrices = <MatrixInt22 *>malloc(self.num_cusps()*sizeof(MatrixInt22))

            for n in range(self.num_cusps()):
                for i,j in [(0,0),(0,1),(1,0),(1,1)]:
                    matrices[n][i][j] = 1 if i == j else 0

            matrices[which_cusp][0][0] = a
            matrices[which_cusp][0][1] = b
            matrices[which_cusp][1][0] = c
            matrices[which_cusp][1][1] = d
            

            result = change_peripheral_curves(self.c_triangulation, matrices)
            
            if result == func_bad_input:
                raise ValueError, 'Peripheral data ((%d, %d), (%d,%d)) not acceptable' % (a,b,c,d)

            free(matrices)
            
        else:
            if self.num_cusps() == 1 and len(peripheral_data) == 2:
                self.set_peripheral_curves(peripheral_data, 0)
                return 
            if len(peripheral_data) > self.num_cusps():
                raise IndexError, 'Provided more peripheral data than there are cusps.'
            for i, basis in enumerate(peripheral_data):
                self.set_peripheral_curves(basis, i)

        self._cache = {}


    def dual_curves(self, max_segments=6):
        """
        Constructs a *reasonable* selection of simple closed curves in
        a manifold's dual 1-skeleton.  In particular, it returns those
        that appear to represent geodesics. The resulting curves can
        be drilled out.

        >>> M = Manifold('m015')
        >>> curves = M.dual_curves()
        >>> curves
        [  0: orientation-preserving curve of length (0.562399148646-2.81543088521j),
           1: orientation-preserving curve of length (1.12479829729+0.652323536768j),
           2: orientation-preserving curve of length (1.26080401747+1.97804689023j),
           3: orientation-preserving curve of length (1.58826932598+1.67347167369j),
           4: orientation-preserving curve of length (1.68719744594+2.81543088521j)]

        Each curve is returned as a dictionary with these keys
        
        >>> curves[0].keys()
        ['complete length', 'index', 'parity', 'filled length']

        We can drill out any of these curves to get a new manifold
        with one more cusp.

        >>> N = M.drill(curves[0])
        >>> (M.num_cusps(), N.num_cusps())
        (1, 2)
        
        By default, this function only finds curves of length 6; this
        can be changed by specifying the optional argument
        max_segments

        >>> M.dual_curves(max_segments=2)
        [  0: orientation-preserving curve of length (0.562399148646-2.81543088521j)]
        """
        cdef int i, num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_MatrixParity parity
        cdef Complex complete_length, filled_length

        if self.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'
        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)
        result = []
        for i from 0 <= i < num_curves:
            info_dict = {}
            get_dual_curve_info(curve_list[i], 
                           &complete_length,
                           &filled_length,
                           &parity)
            info_dict['index'] = i
            info_dict['parity'] = parity
            info_dict['filled length'] = C2C(filled_length)
            info_dict['complete length'] = C2C(complete_length)
            D = DualCurveDict(info_dict)
            D.max_segments = max_segments
            result.append(D)
        free_dual_curves(num_curves, curve_list)
        return ListOnePerLine(result)
    
    def length_spectrum(self, cutoff=1.0):
        """
        M.length_spectrum(cutoff=1.0)

        Print a list of geodesics (with multiplicities) of length
        up to the specified cutoff value. (The default cutoff is 1.0.)
        """
        try:
            D = DirichletDomain(self)
        except:
            raise RuntimeError, 'Length spectrum not available: '\
                                'no Dirichlet Domain.'
        spectrum = D.length_spectrum_dicts(cutoff_length=cutoff)
        print '%-4s %-32s %-12s  %s'%('mult', 'length',
                                      'topology', 'parity')
        for curve in spectrum:
            print '%-4d %-32s %-14s%s'%(
                curve['multiplicity'],
                curve['length'],
                curve['topology'],
                curve['parity'] )

    def chern_simons(self, accuracy=False):
        """
        Returns the Chern-Simons of the manifold, if it is known

        >>> M = Manifold('m015')
        >>> M.chern_simons()
        -0.15320413329715188

        If the flag accuracy is set to True, then it returns the
        volume of the manifold together with the number of digits of
        accuracy as *estimated* by SnapPea.

        >>> M.chern_simons(True)
        (-0.15320413329715188, 12)
        """

        cdef Boolean is_known, requires_initialization
        cdef double CS
        cdef int precision

        if self.c_triangulation is NULL: return 0
        solution_type = self.solution_type()
        if solution_type in ("not attempted", "no solution found"):
            raise ValueError, 'Solution type is: %s'%solution_type

        get_CS_value(self.c_triangulation, &is_known, &CS, &precision, &requires_initialization)

        if not is_known:
            raise ValueError, "Chern-Simons invariant isn't currently known"        
        if accuracy:
            return (CS, precision)
        else:
            return CS

        
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

        if isinstance(which_curve,DualCurveDict):
            max_segments = which_curve.max_segments
            which_curve = which_curve['index']


        new_name = self.name()+'-%d'%which_curve
        c_new_name = new_name

        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)

        if which_curve not in range(num_curves):
            raise IndexError, "Drilling curve requested is not in range(%d)." % num_curves
        
        c_triangulation = drill_cusp(self.c_triangulation,
                                     curve_list[which_curve],
                                     c_new_name)
        free_dual_curves(num_curves, curve_list)

        if c_triangulation == NULL:
            raise RuntimeError, "Curve is not isotopic to a geodesic."
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
            raise ValueError, 'Manifolds must be non-empty.'

        if return_isometries:
            result = compute_isometries(self.c_triangulation, other.c_triangulation, 
                                        &are_isometric, &isometries, NULL)
        else:
            result = compute_isometries(self.c_triangulation, other.c_triangulation, 
                                        &are_isometric, NULL, NULL)
            
        if FuncResult[result] == 'func_bad_input':
            raise ValueError, "Dehn filling coefficients must be relatively prime integers."

        if FuncResult[result] == 'func_failed':
            raise RuntimeError, "SnapPea failed to determine whether the manifolds are isometric."

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
        if self.c_triangulation is NULL: return False
        two_bridge(self.c_triangulation, &is_two_bridge, &p, &q)        
        return (p,q) if  is_two_bridge else False

        

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

cdef class CFundamentalGroup:
    cdef c_GroupPresentation *c_group_presentation
    cdef c_Triangulation *c_triangulation
    cdef readonly num_cusps
        
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

    def __new__(self, Triangulation triangulation,
                      simplify_presentation = True,
                      fillings_may_affect_generators = True,
                      minimize_number_of_generators = True):
        if triangulation.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'
        copy_triangulation(triangulation.c_triangulation, &self.c_triangulation)
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
        if type(word) != types.StringType:
            raise TypeError, "Words must be represented as Python strings."
        word_list = []
        generators = self.generators()
        for letter in word:
            try:
                if letter.islower():
                    word_list.append(1 + generators.index(letter))
                else:
                    word_list.append(-1 - generators.index(letter.lower()))
            except ValueError:
                raise RuntimeError, "Word contains a non-generator."
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
                            
    def num_orig_gens(self):
        """
        Return the number of geometric generators (before simplification).
        """
        return fg_get_num_orig_gens(self.c_group_presentation)

    def generators(self):
        """
        Return the letters representing the generators in the presentation.
        """
        return [ Alphabet[i] for i in range(1, 1 + self.num_generators()) ]

    def relators(self, verbose_form = False):
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
            word = self.c_word_as_string(relation)
            if verbose_form:
                word = "*".join([a if a.islower() else a.lower() + "^-1" for a in list(word)])
            relation_list.append(word)
            fg_free_relation(relation)
        return relation_list

    def meridian(self, int which_cusp=0):
        """
        Returns a word representing a conjugate of the current meridian for
        the given cusp.  Guaranteed to commute with the longitude for the same
        cusp.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.meridian(0)
        'aaba'
        >>> G.meridian(-1)  # The last cusp
        'baaba'
        """
        try:
            which_cusp = range(self.num_cusps)[which_cusp]
        except IndexError:
            raise IndexError, "Specified cusp (%s) does not exist."%which_cusp
        return self.c_word_as_string(
            fg_get_meridian(self.c_group_presentation, which_cusp))

    def longitude(self, int which_cusp=0):
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
            raise IndexError, "Specified cusp (%s) does not exist."%which_cusp

        return self.c_word_as_string(
            fg_get_longitude(self.c_group_presentation, which_cusp))

    def peripheral_curves(self):
        """
        Returns a list of meridian-longitude pairs for all cusps.

        >>> G = Manifold('m125').fundamental_group()
        >>> G.peripheral_curves()
        [('aaba', 'abb'), ('baaba', 'Ba')]
        """
        return [ (self.meridian(n), self.longitude(n))
                 for n in range(self.num_cusps) ]

    def magma_string(self):
        """
        Returns a string which will define this group within MAGMA.
        """
        return "Group<" + ",".join(self.generators()) + "|" + ", ".join(self.relators(verbose_form = True)) + ">"

    def gap_string(self):
        """
        Returns a string which will define this group within GAP.
        """
        gens = ", ".join(self.generators())
        gen_names = ", ".join(['"' + x + '"' for x in self.generators()])
        relators = ", ".join(self.relators(verbose_form = True))
        assignments = "".join(["%s := F.%d; " % (x, i+1) for (i, x) in enumerate(self.generators())])
        return "CallFuncList(function() local F, %s; F := FreeGroup(%s); %s  return F/[%s]; end,[])"  % (gens, gen_names, assignments, relators)

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
    operator is named"+", according to Python conventions).

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
            o31 = matrix([[O[0][0], O[0][1], O[0][2]],
                          [O[1][0], O[1][1], O[2][2]],
                          [O[2][0], O[2][1], O[2][2]]])
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

cdef class CDirichletDomain:
    cdef WEPolyhedron *c_dirichlet_domain
    cdef c_Triangulation *c_triangulation

    def __new__(self, Manifold manifold,
                      vertex_epsilon=10.0**-8,
                      displacement = [0.0, 0.0, 0.0],
                      centroid_at_origin=True,
                      maximize_injectivity_radius=True):
        cdef double c_displacement[3]
        if manifold.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'
        for n from 0 <= n < 3:
            c_displacement[n] = displacement[n] 
        copy_triangulation(manifold.c_triangulation, &self.c_triangulation)
        self.c_dirichlet_domain = Dirichlet_with_displacement(
            self.c_triangulation,
            c_displacement, 
            vertex_epsilon,
            centroid_at_origin,
            Dirichlet_keep_going,
            maximize_injectivity_radius )
        if self.c_dirichlet_domain == NULL:
            raise RuntimeError, "Dirichet construction failed."
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
        Return a list of dictionaries describing the short geodesics
        up to the specified cutoff length.  The keys are 'length',
        'parity', 'topology', and 'multiplicity'.  The length is the
        complex length; the parity specifies whether orientation is
        preserved; and topology distinguishes between circles and
        mirrored intervals.
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
            spectrum.append({
                    'length' : C2C(geodesics[n].length),
                    'parity' : MatrixParity[geodesics[n].parity],
                    'topology' : Orbifold1[geodesics[n].topology],
                    'multiplicity': geodesics[n].multiplicity })
        free_length_spectrum(geodesics)
        return spectrum

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
            raise RuntimeError, "PolyhedronViewer was not imported."


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
    """
    pass

# Cusp Neighborhoods

cdef class CCuspNeighborhood:
    cdef CuspNeighborhoods *c_cusp_neighborhood
    cdef c_Triangulation *c_triangulation

    def __new__(self, Manifold manifold):
        if manifold.c_triangulation is NULL:
            raise ValueError, 'Triangulation is empty.'
        copy_triangulation(manifold.c_triangulation, &self.c_triangulation)
        self.c_cusp_neighborhood = initialize_cusp_neighborhoods(
            self.c_triangulation)
        if self.c_cusp_neighborhood == NULL:
            raise RuntimeError, "Cusp Neighborhood construction failed."
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

    def check_index(self, which_cusp):
        N = int(which_cusp)
        if 0 <= N < self.num_cusps():
            return N
        else:
            raise IndexError, "Specified cusp (%s) does not exist."%which_cusp
        
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
            raise RuntimeError, "Horoball construction failed."
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
            raise RuntimeError, "Ford domain construction failed."
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
            raise RuntimeError, "Triangulation construction failed."
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

    def view(self, which_cusp=0, cutoff=0.1):
        if HoroballViewer:
            self.viewer = HoroballViewer(
                self, cutoff, which_cusp,
                title='Cusp neighborhood #%s of %s'%(
                    which_cusp,
                    self.manifold_name
                    ))
        else:
            raise RuntimeError, "HoroballViewer was not imported."
        
class CuspNeighborhood(CCuspNeighborhood):
    """
    A CuspNeighborhood object represents an equivariant collection of disjoint
    horoballs that project to cusp neighborhoods.

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
    
    def __new__(self, is_full_group, owns_c_symmetry_group):
        self.c_symmetry_group = NULL 
        self._is_full_group = is_full_group
        self._owns_c_symmetry_group = owns_c_symmetry_group

    def __dealloc__(self):
        #if self._owns_c_symmetry_group:
        #    free_symmetry_group(self.c_symmetry_group)
        pass

    cdef _set_c_symmetry_group(self, c_SymmetryGroup * c_symmetry_group):
        if c_symmetry_group is NULL:
            raise ValueError, "Tried to create an *empty* SymmetryGroup"
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
            theText = 'D%d'%(self.order()/2)
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
            raise ValueError, "Symmetry group is not abelian"

        coeffs = []
        for n from 0 <= n < A.num_torsion_coefficients:
                coeffs.append(A.torsion_coefficients[n])

        free_abelian_group(A)
        return AbelianGroup(coeffs)
            
    
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
            raise ValueError, "Symmetry group is not polyhedral"

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
            raise ValueError, "Symmetry group is not a nontrivial, nonabelian direct product"

        cdef c_SymmetryGroup* c_factor_0
        cdef c_SymmetryGroup* c_factor_1
        cdef SymmetryGroup factor_0
        cdef SymmetryGroup factor_1
        
        c_factor_0 = get_symmetry_group_factor(self.c_symmetry_group, 0)
        c_factor_1 = get_symmetry_group_factor(self.c_symmetry_group, 1)
        
        factor_0, factor_1 = SymmetryGroup(True, False), SymmetryGroup(True, False)
        factor_0._set_c_symmetry_group(c_factor_0), factor_1._set_c_symmetry_group(c_factor_1)
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
            raise ValueError, "Full symmetry group not known"

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
            raise ValueError, "Full symmetry group not known"

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
                raise ValueError, "Symmetry group has only %d elements" % order

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
        free_isometry_list(isometries)
        return ans


# get_triangulation

split_filling_info = re.compile("(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))+)")
is_census_manifold = re.compile("([msvxy])([0-9]+)$")
is_torus_bundle = re.compile("b([+-no])([+-])([lLrR]+)$")
is_knot_complement = re.compile("(?P<crossings>[0-9]+)_(?P<index>[0-9]+)$")
is_link_complement1 = re.compile("(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$")
is_link_complement2 = re.compile("(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$")
is_link_complement3 = re.compile("[lL]([0-9]+)$")
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)$')
is_braid_complement = re.compile("braid(\[[0-9, -]+\])$")
is_DT_exterior = re.compile("DT(\[[0-9, -]+\])$")
is_census_knot = re.compile("K[2-7]_([0-9]+)$")

#Orientability.orientable = 0
spec_dict = {'m' : (5, 0),
             's' : (6, 0),
             'v' : (7, 0),
             'x' : (6, 1),
             'y' : (7, 1)}

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

cdef c_Triangulation* get_triangulation(spec) except ? NULL:
    cdef c_Triangulation* c_triangulation = NULL
    cdef LRFactorization* gluing
    cdef int LRlength

    # get filling info, if any
    m = split_filling_info.match(spec)
    if m:
        real_name = m.group(1)
        fillings = re.subn("\)\(", "),(", m.group(2))[0]
        fillings = eval( "[" + fillings + "]" )
    else:
        real_name = spec
        fillings = ()

    # Step 1. Check for a census manifold
    m = is_census_manifold.match(real_name)
    if m:
        num_tet, orientable = spec_dict[m.group(1)]
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, num_tet, orientable, int(m.group(2)))
        set_cusps(c_triangulation, fillings)
        return c_triangulation

     # Step 2. Check for a punctured torus bundle 
    m = is_torus_bundle.match(real_name)
    if m:
        LRstring = m.group(3).upper()
        LRlength = len(LRstring)
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
        strncpy(gluing.LR_factors, LRstring, 1+LRlength)
        c_triangulation =  triangulate_punctured_torus_bundle(gluing)
        free_LR_factorization(gluing)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 3. Check for a Rolfsen link complement
    filename = None
    m = is_knot_complement.match(real_name)
    if m:
        filename = "L1%.2d%.3d" % (int(m.group("crossings")),
                                   int(m.group("index")))
    m = is_link_complement1.match(real_name)
    if m:
        filename = "L%.1d%.2d%.3d" % (int(m.group("components")),
                                      int(m.group("crossings")),
                                      int(m.group("index")))
    m = is_link_complement2.match(real_name)
    if m:
        filename = "L%.1d%.2d%.3d" % (int(m.group("components")),
                                      int(m.group("crossings")),
                                      int(m.group("index")))
    m = is_link_complement3.match(real_name)
    if m:
        filename = "L" + m.group(1)
    if filename:
        tarpath =  'ChristyLinks/%s'%filename
        try:
            filedata = Christy_links.extractfile(tarpath).read()
            c_triangulation = read_triangulation_from_string(filedata)
        except: 
            raise IOError, "The link complement %s was not found."%real_name
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 4. Check for a Hoste-Thistlethwaite knot.
    m = is_HT_knot.match(real_name)
    if m:
        c_triangulation = get_HT_knot(int(m.group("crossings")),
                             m.group("alternation"),
                             int(m.group("index")))
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 5. Check for a census knot.
    m = is_census_knot.match(real_name)
    if m:
        tarpath =  'CensusKnots/%s'%real_name
        try:
            filedata = Census_Knots.extractfile(tarpath).read()
            c_triangulation = read_triangulation_from_string(filedata)
        except: 
            raise IOError, "The census knot %s was not found."%real_name
        set_cusps(c_triangulation, fillings)
        return c_triangulation
        
    # Step 6. See if a (fibered) braid complement is requested

    m = is_braid_complement.match(real_name)
    if m:
        word = eval(m.group(1))
        num_strands = max([abs(x) for x in word]) + 1
        c_triangulation = get_fibered_manifold_associated_to_braid(num_strands, word)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 7.  See if a knot exterior is requested via its
    # Dowker-Thistlethwaite code:

    m = is_DT_exterior.match(real_name)
    if m:
        word = eval(m.group(1))
        c_triangulation = get_link_exterior_from_DT(word)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 8. If all else fails, try to load a manifold from a file.
    try:
        locations = [os.curdir, os.environ["SNAPPEA_MANIFOLD_DIRECTORY"]]
    except KeyError:
        locations = [os.curdir]
    found = 0
    for location in locations:
        pathname = os.path.join(location, real_name)
        if os.path.isfile(pathname):
            file = open(pathname, "r")
            first_line = file.readline()[:-1]
            file.close()
            if first_line.find("% Link Projection") > -1:
                c_triangulation = triangulate_link_complement_from_file(pathname, "")
            else:
                c_triangulation = read_triangulation(pathname)
            set_cusps(c_triangulation, fillings)
            return c_triangulation

    # Step 9. Give up.
    raise IOError, "The manifold file %s was not found.\n%s"%(
        real_name, triangulation_help%'Triangulation or Manifold')
        
cdef int set_cusps(c_Triangulation* c_triangulation, fillings) except -1:
    if c_triangulation == NULL:
        return 0
    if len(fillings) > 0:
        num_cusps = get_num_cusps(c_triangulation) 
        if len(fillings) > num_cusps:
            raise ValueError, "The number of fillings specified exceeds the number of cusps."
        for i in range(len(fillings)):
            meridian, longitude = fillings[i]
            is_complete = (meridian == 0 and longitude == 0)
            set_cusp_info(c_triangulation, i, is_complete, meridian, longitude)
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
    size = (1+crossings)/2
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
    size = (1 + crossings)/2
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
        raise ValueError, """
        You have specified a Hoste-Thistlethwaite knot with an
        inappropriate index or number of crossings."""

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
    name = "DT"+repr(DT)
    set_triangulation_name(c_triangulation, name)
    free(DT_array)
    return c_triangulation

cdef c_Triangulation* get_HT_knot(crossings, alternation, index) except ? NULL:
    cdef c_Triangulation* c_triangulation
    DT = get_HT_knot_DT(crossings, alternation, index)
    c_triangulation = get_link_exterior_from_DT(DT)
    name = "%d" % crossings + alternation + "%d" % index
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
    name = "%d" % crossings + alternation + "%d" % (index_within_crossings + 1)
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

    def next(self):
        if self.index >= self.stop:
            raise StopIteration
        self.index = self.index + self.step
        return self[self.index-self.step]

    def __len__(self):
        return self.length
    
    # subclasses override this
    def __getitem__(self, n):
        pass

#  Cusped Census

Orientable_lengths = (301, 962, 3552, 301+962+3552)
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
    five_length, six_length, seven_length, length = Orientable_lengths
    orientability = Orientability.index('orientable')

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)

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
        else:
            raise IndexError, "Index out of range."
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, num_tet, self.orientability, census_index)
        if c_triangulation == NULL:
            print num_tet, census_index
            raise RuntimeError, "SnapPea failed to read census manifold."
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        return result

class OrientableCuspedCensus(CuspedCensus):
    """
    Iterator/Sequence for orientable manifolds in the SnapPea
    Cusped Census.

    >>> C = OrientableCuspedCensus()
    >>> for M in C[:5]:   # Just the first 5 manifolds
    ...     print M, M.volume()
    m003(0,0) 2.02988321282
    m004(0,0) 2.02988321282
    m006(0,0) 2.56897060094
    m007(0,0) 2.56897060094
    m009(0,0) 2.66674478345
    """

class NonorientableCuspedCensus(CuspedCensus):
    """
    Iterator/Sequence for nonorientable manifolds in the SnapPea
    Cusped Census.
    """
    five_length, six_length, seven_length, length = Nonorientable_lengths
    orientability = Orientability.index('nonorientable')

    def __init__(self, indices=(0, length, 1)):
        Census.__init__(self, indices)

    def lookup(self, n):
        return five_tet_nonorientable[n]

# Closed Census

class OrientableClosedCensus(Census):
    """
    Iterator/Sequence for orientable closed manifolds in the SnapPea
    Closed Census.

    >>> C = OrientableClosedCensus()
    >>> M = C[0]
    >>> M.volume() # The smallest hyperbolic manifold!
    0.94270736277692779
    """
    data = None
    def __init__(self, indices=(0,11031,1)):
        if OrientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedOrientableDistinct.txt')
            closed_orientable = open(datafile)
            OrientableClosedCensus.data = closed_orientable.readlines()
            closed_orientable.close()
        self.length = len(OrientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = OrientableClosedCensus.data[n].split()
        code = rev_spec_dict[(int(num_tet), 0)]
        spec = '%s%s(%s,%s)'%(code,index,m,l)
        return Manifold(spec)

class NonorientableClosedCensus(Census):
    """
    Iterator/Sequence for orientable closed manifolds in the SnapPea
    Closed Census.
    """
    data = None
    def __init__(self, indices=(0,17,1)):
        if NonorientableClosedCensus.data is None:
            datafile = os.path.join(closed_census_directory,
                                    'ClosedNonorientableDistinct.txt')
            closed_nonorientable = open(datafile)
            NonorientableClosedCensus.data = closed_nonorientable.readlines()
            closed_nonorientable.close()
        self.length = len(NonorientableClosedCensus.data)
        Census.__init__(self, indices)

    def __getitem__(self,n):
        if isinstance(n, slice):
            return self.__class__(n.indices(self.length))
        volume, num_tet, index, m, l = NonorientableClosedCensus.data[n].split()
        code = rev_spec_dict[(int(num_tet), 1)]
        spec = '%s%s(%s,%s)'%(code,index,m,l)
        return Manifold(spec)

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
    Hoste-Thistlethwaite tables.
    """

class NonalternatingKnotExteriors(KnotExteriors):
    """
    Iterator/Sequence for nonAlternating knot exteriors from the
    Hoste-Thistlethwaite tables.
    """
    length = sum(Nonalternating_numbers.values())
    alternation = 'n'

    def __init__(self, indices=(0, sum(Nonalternating_numbers.values()), 1)):
        Census.__init__(self, indices)

census_knot_numbers = [0, 0, 1, 2, 4, 22, 43, 129]

class CensusKnots(Census):
    """
    Iterator/Sequence for knot exteriors in the SnapPea Census as
    tabulated by Callahan, Dean, Weeks, Champanerkar, Kofman and
    Patterson.

    >>> K = CensusKnots()
    >>> M = K[75]
    >>> M
    K7_4(0,0)
    >>> M.volume()
    3.6352511866719941
    >>> Manifold('v0114').volume()
    3.6352511866719941
    """
    length = sum(census_knot_numbers)

    def __init__(self, indices=(0, sum(census_knot_numbers), 1)):
        Census.__init__(self, indices)

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
                return  Manifold(name)
            else:
                raise IndexError, 'There are only 201 census knots.'
                
class LinkExteriors(Census):
    """
    Census of links/knots using the classical numbering system of
    Tait/Conway/Rolfsen/Christy.  Includes knots through 11 crossings,
    and links through 10 crossings.  Mostly useful just for links as
    the Hoste-Thistlethwaite table of knots is much more extensive.
    Takes as argument the number of components.

    >>> C = LinkExteriors(2)    # 2 component links
    >>> len(C)
    273
    >>> C[20]
    8^2_8(0,0)(0,0)
    >>> for M in LinkExteriors(5):
    ...     print M, M.volume()
    10^5_1(0,0)(0,0)(0,0)(0,0)(0,0) 14.6030607534
    10^5_2(0,0)(0,0)(0,0)(0,0)(0,0) 12.8448530047
    10^5_3(0,0)(0,0)(0,0)(0,0)(0,0) 10.1494160641
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
         if not (1 <= components < len(self.num_links) ):
            raise IndexError, "No data on links with %s many components." % components

         self.components = components

         self.length = sum(self.num_links[components])

         Census.__init__(self, (indices[0], min(self.length, indices[1]), indices[2]))

    def __repr__(self):
        return 'Christy census of link complements in S^3 with %s components' % self.components

    def __getitem__(self,j):
        if isinstance(j, slice):
            return self.__class__(n.indices(self.length))
        so_far = 0
        for k in range(self.max_crossings + 1):
            n =  self.num_links[self.components][k]
            so_far = so_far + n
            if so_far > j:
                l = j - so_far + n + 1
                name = "%d^%d_%d" % (k, self.components, l) if self.components > 1 \
                       else "%d_%d" % (k,  l)
                M =  Manifold(name)
                M.set_name(name)
                return M
                
# Creating fibered manifolds from braids

cdef c_Triangulation*  get_fibered_manifold_associated_to_braid(num_strands, braid_word):
    if num_strands < 2:
        raise ValueError, "Must have at least 2 strands."
    allowed_letters = range(1,num_strands) + range(-num_strands+1, 0)
    if False in [b in allowed_letters for b in braid_word]:
        raise ValueError, "Invalid braid word."

    cdef int* word
    cdef c_Triangulation* c_triangulation

    n = len(braid_word)
    word = <int*>malloc(n*sizeof(int))
    for  i, w in enumerate(braid_word):
        word[i] = w
    for i in range(n):
        c_triangulation = fibered_manifold_associated_to_braid(num_strands, n, word)
    free(word)
    name = "braid" + repr(braid_word)
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
        raise RuntimeError, "Couldn't allocate crossing table."

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
        raise RuntimeError, 'Could not create triangulation.'

    # The triangulation must have a name or SnapPea will segfault when
    #   trying to copy the triangulation.
    set_triangulation_name(c_triangulation, 'unnamed link');
    return c_triangulation

# Code for interacting with the OS X GUI SnapPea.

SnapPeaX_name = "SnapPea"

def execute_applescript( data ):
    script_name = tempfile.mktemp() + ".scpt"
    open(script_name, "w").write(data)
    script = os.popen("osascript " + script_name)
    ans = script.read()
    os.remove(script_name)
    return ans 

def activate_SnapPeaX():
    script = """\
    tell application "%s"
           activate
    end tell""" % SnapPeaX_name
    execute_applescript(script)

# Commands for interacting with SnapPeaX
             
def get_from_SnapPeaX():
    file_name = tempfile.mktemp() + ".tri"
    script = """\
    on suffix(s, i)
       set len to length of s
       if len < i then
            return s
       else
             return text -i thru -1 of s
       end if
     end suffix

     set f to POSIX file "%s"
     tell application "%s"
          set winName to (name of window 1)
      end tell
      set tail to suffix(winName, 11)
      if tail = "Link Editor" then
          return 0
      end if
      tell application "%s"
           try
              set doc to (document of window 1)
              save doc in f
                  return 1
              on error
                  return 0
              end try
      end tell"""  % (file_name, SnapPeaX_name, SnapPeaX_name)
    ans = int(execute_applescript(script))
    if not ans:
        return None
    manifold =  Manifold(file_name)
    os.remove(file_name)
    return manifold

def get_all_from_SnapPeaX():
    num_docs = int(execute_applescript(
        'tell application "%s" to count of the documents' % SnapPeaX_name))
    file_names = [tempfile.mktemp() + ".tri" for i in range(num_docs)]
    script = """\
    set savenames to {%s}
    set usednames to {}

    tell application "SnapPea"
           repeat with win in windows
                 set winname to (name of win)
                     set tail to (text -4 thru -1 of winname)
                           if tail = "Main" then
                                set savename to item 1 of savenames
                                set savenames to rest of savenames
                                save (document of win) in POSIX file savename
                                set usednames to usednames & savename
                            end if
	   end repeat
    end tell
    usednames""" % repr(file_names)[1:-1].replace("'", '"')

    used_names = [name.strip() for name in execute_applescript(script).split(",")]
    manifolds =  [Triangulation(file_name) for file_name in used_names]
    [os.remove(file_name) for file_name in used_names]
    return manifolds

# So that SnapPeaX can raise the correct Terminal Window, we signal
# which window this is by changing the title.  

def connect_to_SnapPeaX():
    script = """\
    tell application "Terminal"
    set custom title of front window to "SnapPeaPython"
    end tell
    
    tell application "%s" to activate
    tell application "Terminal" to activate 
    """ % SnapPeaX_name
    execute_applescript(script)
    
def detach_from_SnapPeaX():
    script = """\
    tell application "Terminal"
         set custom title of front window to ""
    end tell"""
    execute_applescript(script)
