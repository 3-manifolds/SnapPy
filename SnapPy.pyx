import os, sys, operator, types, re, gzip, struct, tempfile
from signal import signal, SIGINT, SIG_DFL, SIGALRM
from manifolds import __path__ as manifold_paths

include "SnapPy.pxi"

# We need a matrix class

try:
    from sage.matrix.constructor import matrix
except ImportError:
    from numpy import matrix

# Enable graphical link input, if plink is available.

try:
    import plink
except:
    pass

# Enable OpenGL display of DirichletDomains

try:
    from polyviewer import PolyhedronViewer
except ImportError:
    PolyhedronViewer = None

# SAGE interaction

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

# Paths
manifold_path = manifold_paths[0] + os.sep
closed_census_directory = os.path.join(manifold_path, 'ClosedCensusData')
link_directory = os.path.join(manifold_path, 'ChristyLinks')
table_directory = os.path.join(manifold_path, 'HTWKnots')

# Implementation of the SnapPea UI functions and their global variables

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

cdef public void uLongComputationBegins(char *message, Boolean is_abortable):
    global gLongComputationCancelled
    global gLongComputationInProgress
    # Set SnapPea's flags
    gLongComputationCancelled = False
    gLongComputationInProgress = True
    # Install the SnapPea handler on SIGINT and SIGALRM
    signal(SIGINT, SnapPea_handler)
    signal(SIGALRM, SnapPea_handler)

cdef public c_FuncResult uLongComputationContinues():
    global gLongComputationCancelled
    # While a SnapPea function is executing, Python saves all of its
    # calls to interrupt handlers on a list of "Pending Calls".
    # This forces any waiting handlers to be called.
    Py_MakePendingCalls()
    if gLongComputationCancelled:
        return func_cancelled
    else:
        return func_OK

cdef public void uLongComputationEnds():
    global gLongComputationCancelled
    global gLongComputationInProgress
    # Restore Python's default signal handler on SIGINT and SIGALRM
    signal(SIGINT, python_handler)
    signal(SIGALRM, python_handler)
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
        # probably means we're within SAGE
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

    # PARI views the input to matsnf0 as square.
    if m > n:
        for i in range(m - n):
            result.remove(0)
    cgiv(pari_vector)
    cgiv(pari_matrix)
    return result

# Enum conversions
CuspTopology = ['torus cusp', 'Klein bottle cusp', 'unknown']
MatrixParity = ['orientation-reversing', 'orientation-preserving']
Orientability = ['orientable', 'nonorientable', 'unknown']
FuncResult = ['func_OK', 'func_cancelled', 'func_failed', 'func_bad_input']

# global functions

def check_SnapPea_memory():
    verify_my_malloc_usage()

def doc(X=None):
    if X is None:
        print __doc__
    else:
        print X.__doc__

# SnapPea Classes

# Abelian Groups

cdef class AbelianGroup:
    """
    An AbelianGroup object represents a finitely generated abelian group,
    usually the first homology group of a SnapPeaX Manifold.
bg
    Instantiate as AbelianGroup([n_1, n_2, ... ]) where the n_i are the
    orders of the cyclic factors (or 0, in the case of an infinite cyclic
    factor).

    Methods:
      Betti_number() --> rank of maximal free abelian subgroup
      order()        --> the order of the group, or the string 'infinite'
      G[n] is the order of the nth cyclic factor
    """

    cdef readonly coefficients

    def __init__(self, coefficients):
        try:
            self.coefficients = list(coefficients)
        except:
            raise RuntimeError, """
        Argument is not a sequence"""
        for c in self.coefficients:
            assert type(c) == types.IntType and c >= 0,\
                'Coefficients must be non-negative integers.\n'
        self.coefficients.sort()

    def __repr__(self):
        factors = ( ['Z' for n in self.coefficients if n == 0] +
                    ['Z/%d'%n for n in self.coefficients if n > 1] )
        return ' + '.join(factors)

    def __len__(self):
        return len(self.coefficients)
    
    def __getitem__(self, i):
        return self.coefficients[i]

    def Betti_number(self):
        return len([n for n in self.coefficients if n == 0])

    def rank(self):
        return len(self.coefficients)
    
    def order(self):
        det = reduce(operator.mul, [1] + self.coefficients)
        if det == 0:
            return 'infinite'
        else:
            return det

# Triangulations

cdef class Triangulation:
    """
    A Triangulation object represents the interior of a 3-manifold
    with non-empty boundary, each component of which is a torus.  The
    3-manifold comes equipped with an ideal triangulation.  Two
    Triangulations are equal ("==") if they represent combinatorially
    isomorphic triangulations.

    A Triangulation does NOT have any geometric structure.  The
    subclass Manifold adds the geometric structure to a Triangulation,
    and is the object one usually wants to work with.

    Convention: methods which change the triangulation always return
    a new Triangulation.

    A triangulation can be specified by a string, according to the
    following conventions:

    1. Numbers in parens at the end specify Dehn fillings.  For example
    'm125(1,2)(4,5)' means do (1,2) filling on the first cusp and (4,5)
    filling on the second cusp of the census manifold m125.

    2. Strings of the form 'm123', 's123', 'v123', and so on refer to the
    SnapPea Cusped Census manifolds.

    3. Strings of the form '4_1', '04_1', '4_01', '5^2_6', '6_4^7',
    etc, refer to complements of links in Rolfsen's table.  Similarly
    for 'L20935', 'l104001', etc.

    4. Strings of the form 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'
    refer to the correponding punctured torus bundle.

    5. Strings of the form '11a17' or '12n345' refer to complements of
    knots in the Hoste-Thistlethwaite tables.

    6. Strings of the form 'braid[1,2,-3,4]' creates fibered manifold
    corresponding the given braid.  In other words, think of the braid
    as givening an element of the mapping class group of the
    numStrands-punctured disc.  This function returns the
    corresponding mapping torus.  If you want the braid closure, you
    have to do (1,0) filling of the last cusp.

    If the string is not in any of the above forms it is assumed to be
    the name of a SnapPea manifold file.  The file will be loaded
    if found in the current directory or the path given by the user
    variable SNAPPEA_MANIFOLD_DIRECTORY.       

    """

    cdef c_Triangulation* c_triangulation
    cdef readonly num_cusps, num_or_cusps, num_nonor_cusps, is_orientable
    cdef readonly _cache

    def __new__(self, spec=None):
        cdef c_Triangulation *c_triangulation = NULL
        if spec is not None and spec != 'empty':
            if type(spec) != types.StringType:
                raise TypeError, triangulation_help%self.__class__.__name__
            c_triangulation = get_triangulation(spec)
            if c_triangulation == NULL:
                raise TypeError, "Specified empty manifold"
        if spec is None:
            try:
                dialog = plink.LinkEditor(no_arcs=True)
                dialog.window.mainloop()
            except:
                raise RuntimeError, "Please install PLink to use this feature"
            klp = dialog.SnapPea_KLPProjection() 
            if klp is not None:
                c_triangulation = get_triangulation_from_PythonKLP(klp)

        if c_triangulation != NULL:    
            self.set_c_triangulation(c_triangulation)
            remove_hyperbolic_structures(c_triangulation)

            
    cdef set_c_triangulation(self, c_Triangulation* c_triangulation):
        self.c_triangulation = c_triangulation
        self.num_cusps = get_num_cusps(self.c_triangulation)
        self.num_or_cusps = get_num_or_cusps(self.c_triangulation)
        self.num_nonor_cusps = get_num_nonor_cusps(self.c_triangulation)
        orientability = Orientability[get_orientability(self.c_triangulation)]
        if orientability == 'orientable': self.is_orientable = True
        elif orientability == 'nonorientable': self.is_orientable = False
        else: self.is_orientable = None

        # Answers to potentially hard computations are cached

        self._cache = {}

    def copy(self):
        """
        Returns a copy of this triangulation.
        """
        cdef c_Triangulation* copy_c_triangulation = NULL
        cdef Triangulation new_tri
        
        copy_triangulation(self.c_triangulation, &copy_c_triangulation)
        new_tri = Triangulation()
        new_tri.set_c_triangulation(copy_c_triangulation)
        return new_tri

    def randomize(self):
        """
        Perform random Pachner moves on the triangulation.
        """
        randomize_triangulation(self.c_triangulation)
        self._cache = {}

    def simplify(self):
        """
        Try to simplify the triangulation by doing Pachner moves.
        """
        basic_simplification(self.c_triangulation)
        self._cache = {}

    def save(self, file_name):
        write_triangulation(self.c_triangulation, file_name)

    def __dealloc__(self):
        if self.c_triangulation is not NULL:
            free_triangulation(self.c_triangulation)

    def __richcmp__(Triangulation self, Triangulation other, case):
        cdef c_Triangulation *c_triangulation1
        cdef c_Triangulation *c_triangulation2
        cdef Boolean answer
        if case != 2:
            return NotImplemented
        if type(self) != type(other):
            return False
        if same_triangulation(self.c_triangulation, other.c_triangulation):
            return True
        else:
            return False

    def __repr__(self):
        if self.c_triangulation is NULL:
            return 'Empty Triangulation'
        else:
            repr = self.name()
            for i in range(self.num_cusps):
                info = self.cusp_info_dict(i)
                if info['complete?']:
                    repr += '(0,0)'
                else:
                    repr += '(%g,%g)'%(info['m'],info['l'])
            return repr
 
    def set_name(self, new_name):
        cdef char* c_new_name = new_name
        if self.c_triangulation is not NULL:
            set_triangulation_name(self.c_triangulation, c_new_name)

    def name(self):
        if self.c_triangulation is not NULL:
            return get_triangulation_name(self.c_triangulation)

    def num_tetrahedra(self):
        return get_num_tetrahedra(self.c_triangulation)
    
    def dehn_fill(self, meridian, longitude, which_cusp=0):
        """
        dehn_fill(self, meridian, longitude, which_cusp=0)
        
        Assigns the specified Dehn filling coefficients.
        Does not return a new Triangulation.
        """
        if not 0 <= which_cusp < self.num_cusps:
            raise IndexError, "Specified cusp does not exist"
        complete = ( meridian == 0 and longitude == 0)
        set_cusp_info(self.c_triangulation,
                      which_cusp, complete, meridian, longitude)
        self._cache = {}

    def cusp_info_dict(self, int which_cusp = 0):
        cdef c_CuspTopology topology
        cdef Boolean is_complete,
        cdef double m, l
        cdef Complex initial_shape, current_shape
        cdef int initial_shape_precision, current_shape_precision,
        cdef Complex initial_modulus, current_modulus
        if which_cusp >= self.num_cusps or which_cusp < 0:
            raise IndexError, """
        There are %d cusps!"""%self.num_cusps
        get_cusp_info(self.c_triangulation, which_cusp,
                      &topology, &is_complete, &m, &l,
                      &initial_shape, &current_shape,
                      &initial_shape_precision, &current_shape_precision,
                      &initial_modulus, &current_modulus)
        return {'topology' : CuspTopology[topology],
                'complete?' : is_complete,
                'm' : m, 'l' : l,
                'initial shape' : C2C(initial_shape),
                'current shape' : C2C(current_shape),
                'initial shape precision' : initial_shape_precision,
                'current shape precision' : current_shape_precision,
                'initial modulus' : C2C(initial_modulus),
                'current modulus' : C2C(current_modulus)}

    def cusp_info(self):
        """
        Print a list of the cusps, showing the topological type, the
        complex modulus and the current Dehn filling coefficients.
        """

        for i in range(self.num_cusps):
            info_dict = self.cusp_info_dict(i)
            if info_dict['complete?']:
                print 'Cusp %-2d: %s, not filled'%\
                    (i, info_dict['topology'])
            else:
                print 'Cusp %-2d: %s with Dehn filling coeffients M = %g, L = %g'%\
                    (i, info_dict['topology'], info_dict['m'], info_dict['l'])

    def filled_triangulation(self, cusps_to_fill="all"):
        """
        Return a new triangulation where the specified cusps have
        been permently filled in.
        """
        n = self.num_cusps
        if cusps_to_fill == "all":
            cusps_to_fill = [c for c in range(n) if cusp_is_fillable(self.c_triangulation, c)]
                
        if False in [(c in range(n)) for c in cusps_to_fill]:
            raise ValueError, "Specified indices to be filled are beyond the actual number of cusps"
        if 0 in [cusp_is_fillable(self.c_triangulation, c) for c in cusps_to_fill]:
            raise ValueError, "To permanently fill a cusp, the Dehn filling coefficients must be relatively prime integers."

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

        filled_tri = Triangulation()
        filled_tri.set_c_triangulation(c_filled_tri)
        filled_tri.set_name(self.name() + "_filled")

        return filled_tri
        
        
    def edge_valences(self):
        """
        Returns a dictionary whose keys are the valences of the edges
        in the triangulation, and the value associated to a key is the
        number of edges of that valence.
        """
        cdef int c, v = 1
        ans = {}
        while get_num_edge_classes(self.c_triangulation, v, 1) > 0:
            c = get_num_edge_classes(self.c_triangulation, v, 0)
            if c > 0:
                ans[v] = c
            v += 1
        return ans
        

    def gluing_equations(self,form="log"):
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
        for a detailed account of the convention.

        If the optional argument form="rect" is given, then this
        function returns a list of tuples of the form:

           ( [a0, a1,..,a_n], [b_0, b_1,...,b_n], c)

        where this corresponds the the equation

           z0^a0 (1 - z0)^b0 + z1^a1(1 - z1)^b1 + ...  = c

        where c = 1 or -1.  
        """
        
        cdef int **c_eqns
        cdef int num_rows, num_cols
        c_eqns = get_gluing_equations(self.c_triangulation, &num_rows, &num_cols)
        eqns = [ [c_eqns[i][j] for j in range(num_cols)] for i in range(num_rows)]
        free_gluing_equations(c_eqns, num_rows)

        cdef int* eqn
        for i in range(self.num_cusps):
            cusp_info = self.cusp_info_dict(i)
            to_do = [(1,0), (0,1)] if cusp_info["complete?"] else [(cusp_info["m"], cusp_info["l"])]
            for (m, l) in to_do:
                eqn = get_cusp_equation(self.c_triangulation, i, m, l, &num_rows)
                eqns.append([eqn[j] for j in range(num_rows)])
                free_cusp_equation(eqn)

        if form == "log":
            return matrix(eqns)

        if form != "rect":
            raise ValueError, "Equations available in 'log' and 'rect' forms, not what you specified"

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
        homology group of the (Dehn filled) manifold.
        """
        if "homology" in self._cache.keys():
            return self._cache["homology"]
        
        cdef c_AbelianGroup *H
        cdef RelationMatrix R
        cdef int m, n
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
        Returns a FundamentalGroup representing the fundamental group
        of the manifold.  If integer Dehn surgery parameters have been
        set, then the corresponding peripheral element is killed.

        There are three optional arguments all of which default to True.

             simplify_presentation
             fillings_may_affect_generators
             minimize_number_of_generators


        """
        if not "fundamental_group" in self._cache.keys():
            self._cache["fundamental_group"] = FundamentalGroup(self, simplify_presentation, fillings_may_affect_generators, minimize_number_of_generators)
        return self._cache["fundamental_group"]

    def cover(self, permutation_rep):
        """
        Returns a Triangulation representing the finite cover
        specified by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the simplified presentation of the
        fundamental group.  Each permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.

        If within SAGE the permutations can also be of type
        PermutationGroupElement, in which case they act on the set
        range(1, d + 1).  Or, you can specify a GAP or MAGMA subgroup
        of the fundamental group.  For examples, see the docstring for
        Manifold.cover
        """
        cdef RepresentationIntoSn* c_representation
        cdef c_Triangulation* c_triangulation
        cdef Triangulation cover

        # For SAGE, we need to check if we have been given some
        # alternate inputs
        
        if _within_sage:
            if is_GapElement(permutation_rep):
                if permutation_rep.IsSubgroupFpGroup():
                    GG = gap(self.fundamental_group())
                    coset_action = GG.FactorCosetAction(permutation_rep)
                    perms = PermutationGroup(coset_action.Image().GeneratorsOfGroup()).gens()
                    return self.cover(perms)
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
                        raise TypeError, "Homorphism image is not a permutation group"
                else:
                    raise TypeError, "Magma type not recogonized"
                
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
        cover = Triangulation()
        cover.set_c_triangulation(c_triangulation)
        cover.set_name(self.name()+'~')
        free_representation(c_representation,
                            G.num_orig_gens(),
                            self.num_cusps)
        return cover

    def all_covers(self, degree, method=None):
        """
        Returns a list of Triangulations corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        If you are using SAGE, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional argument

        method = "gap"

        If in addition you have Magma installed, you can use it to do
        the heavy-lifting by specifying method = "magma".
        
        """
        cdef RepresentationList* reps
        cdef RepresentationIntoSn* rep
        cdef c_Triangulation* cover
        cdef Triangulation T

        if method:
            if not _within_sage:
                raise RuntimeError, "Only default method of finding subgroups available as your are not using SAGE"
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
            T = Triangulation()
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
        space.  (This awful mess, like, totally belongs in the kernel!)
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
                raise ValueError, """
        Not a valid permutation list"""

        # Initialize
        num_cusps = self.num_cusps
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

    def copy_to_SnapPeaX(manifold):
        """
        Copies the manifold over into the OS X SnapPea GUI.
        """
        file_name = tempfile.mktemp() + ".tri"
        manifold.save(file_name)
        script = """\
        set f to POSIX file "%s"
        tell application "%s"
        open f
        end tell""" % (file_name, SnapPeaX_name)
        execute_applescript(script)
        activate_SnapPeaX()

# Manifolds

cdef class Manifold(Triangulation):
    """
    A Manifold is a Triangulation together with a geometric structure
    defined by assigning shapes to the tetrahedra.

    A manifold can be specified by a string, according to the
    following conventions:

    1. Numbers in parens at the end specify Dehn fillings.  For example
    'm125(1,2)(4,5)' means do (1,2) filling on the first cusp and (4,5)
    filling on the second cusp of the census manifold m125.

    2. Strings of the form 'm123', 's123', 'v123', and so on refer to the
    SnapPea Cusped Census manifolds.

    3. Strings of the form '4_1', '04_1', '4_01', '5^2_6', '6_4^7',
    etc, refer to complements of links in Rolfsen's table.  Similarly
    for 'L20935', 'l104001', etc.

    4. Strings of the form 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'
    refer to the correponding punctured torus bundle.

    5. Strings of the form '11a17' or '12n345' refer to complements of
    knots in the Hoste-Thistlethwaite tables.

    6. Strings of the form 'braid[1,2,-3,4]' create a fibered manifold
    corresponding the given braid.  In other words, think of the braid
    as givening an element of the mapping class group of the
    numStrands-punctured disc.  This function returns the
    corresponding mapping torus.  If you want the braid closure, you
    have to do (1,0) filling of the last cusp.

    If the string is not in any of the above forms it is assumed to be
    the name of a SnapPea manifold file.  The file will be loaded
    if found in the current directory or the path given by the user
    variable SNAPPEA_MANIFOLD_DIRECTORY.
    """

    def __init__(self, spec=None):
        if self.c_triangulation != NULL:
            find_complete_hyperbolic_structure(self.c_triangulation)
            do_Dehn_filling(self.c_triangulation)

    def copy(self):
        """
        Returns a copy of this manifold.
        """
        return Manifold_from_Triangulation(Triangulation.copy(self))

    def filled_triangulation(self, cusps_to_fill="all"):
        """
        Return a new manifold where the specified cusps have been
        permently filled in.  If every cusp is filled, then it returns
        a triangulation rather than manifold since SnapPea can't deal
        with hyperbolic structures in that case.
        """
        filled = Triangulation.filled_triangulation(self, cusps_to_fill)
        if filled.num_cusps == 0:
            return filled
        return Manifold_from_Triangulation(filled)
       
    def fundamental_group(self,
                   simplify_presentation = True,
                   fillings_may_affect_generators = True,
                   minimize_number_of_generators = True):
        """
        Return a HolonomyGroup representing the fundamental group of
        the manifold, together with its holonomy representation.

        There are three optional arguments all of which default to True.

             simplify_presentation
             fillings_may_affect_generators
             minimize_number_of_generators
        
        """

        if not "fundamental_group" in self._cache.keys():
            self._cache["fundamental_group"] = HolonomyGroup(self, simplify_presentation, fillings_may_affect_generators, minimize_number_of_generators)
        return self._cache["fundamental_group"]

    def cover(self, permutation_rep):
        """
        Returns a Manifold representing the finite cover
        specified by a transitive permutation representation.  The
        representation is specified by a list of permutations, one for
        each generator of the simplified presentation of the
        fundamental group.  Each permutation is specified as a list P
        such that set(P) == set(range(d)) where d is the degree of the
        cover.

        If within SAGE the permutations can also be of type
        PermutationGroupElement, in which case they act on the set
        range(1, d + 1).  Or, you can specify a GAP or MAGMA subgroup
        of the fundamental group.     Some examples:

        sage: M = SnapPy.Manifold("m004")

        # The basic method
        sage: N0 = M.cover([[1, 3, 0, 4, 2], [0, 2, 1, 4, 3]])

        # From a Gap subgroup
        
        sage: G = gap(M.fundamental_group())
        sage: H = G.LowIndexSubgroupsFpGroup(5)[9]
        sage: N1 = M.cover(H)
        sage: N0 == N1
        True

        # Or a homomorphism to a permutation group

        sage: f = G.GQuotients(PSL(2,7))[1]
        sage: N2 = M.cover(f)
        sage: N2.volume()/M.volume()
        7.9999999999999947

        # Or maybe we want larger cover coming from the kernel of this.  

        sage: N3 = M.cover(f.Kernel())
        sage: N3.volume()/M.volume()
        167.99999999999858

        # Check the homology against what Gap computes directly
        
        sage: N3.homology().Betti_number()
        32
        sage: len([ x for x in f.Kernel().AbelianInvariants().sage() if x == 0])
        32

        # We can do the same for Magma

        sage: Q, f = G.pQuotient(5, 1, nvals = 2)
        sage: M.cover(f.Kernel()).volume()
        sage: h = G.SimpleQuotients(1, 11, 2, 10^4)[1,1]
        sage: N4 = M.cover(h)
        sage: N2 == N4
        True
        
        """

        cover = Triangulation.cover(self, permutation_rep)
        return Manifold_from_Triangulation(cover, False)

    def all_covers(self, degree, method=None):
        """
        Returns a list of Manifolds corresponding to all of the
        finite covers of the given degree.

        WARNING: If the degree is large this might take a very, very,
        very long time.

        If you are using SAGE, you can use GAP to find the subgroups,
        which is often much faster, by specifying the optional argument

        method = "gap"

        If you have Magma installed, you can used it to do the heavy
        lifting by specifying method = "magma".
        """
        covers = Triangulation.all_covers(self, degree, method)
        return [Manifold_from_Triangulation(cover, False) for cover in covers]

    def volume(self, precision=False):
        """
        Returns the volume of the manifold.  If the flag precision is
        set to True, then it returns the pair (V,p), where V is the
        computed volume of the manifold, and p is the number of digits
        of accuracy estimated by SnapPea.
        """
        cdef int prec
        vol = volume(self.c_triangulation, &prec)
        if precision:
            return (vol, prec)
        else:
            return vol

    def without_hyperbolic_structure(self):
        """
        Returns self as a Triangulation, forgetting the hyperbolic
        structure in the process.
        """
        return Triangulation_from_Manifold(self)

    def tetrahedra_shapes(self, part=None, fixed_alignment=True):
        """
        Gives the shapes of the tetrahedra in the current solution to
        the gluing equations.  Returns a list with one dictionary for each
        tetrahedra with keys

            rect : the shape of the tetrahedra as a point in the complex plane.
            log : the log of shape
            precision: a list of the approximate precisions of the shapes, in order
            (rect re, rect im, log re, log im)

        If the optional variable part is set to one of the above, then
        the function returns only that component of the data.
        
        The if the flag fixed_alignment is set to False, then the
        edges used to report the shape parameters are choosen so as to
        normalize the triangle.
        """        
        cdef double rect_re, rect_im, log_re, log_im
        cdef int prec_rec_re, prec_rec_im, prec_log_re, prec_log_im
        cdef Boolean is_geometric
        
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
                raise ValueError, "Specified non-existent shape data type"
            return [a[part] for a in ans]
            
        return ans

    def cusp_info(self):
        """
        Print a list of the cusps, showing the topological type, the
        complex modulus and the current Dehn filling coefficients.
        """
        for i in range(self.num_cusps):
            info_dict = self.cusp_info_dict(i)
            if info_dict['complete?']:
                print 'Cusp %-2d: complete %s of modulus %s'%\
                    (i, info_dict['topology'],info_dict['current modulus'])
            else:
                print 'Cusp %-2d: %s with Dehn surgery coeffients M = %g, L = %g'%\
                    (i, info_dict['topology'], info_dict['m'], info_dict['l'])

    
    def dehn_fill(self, meridian, longitude, which_cusp=0):
        """
        dehn_fill(self, meridian, longitude, which_cusp=0)
        
        Assigns the specified Dehn filling coefficients and computes
        the associated hyperbolic structure.  Does not return a new
        Manifold.
        """
        Triangulation.dehn_fill(self, meridian, longitude, which_cusp)
        do_Dehn_filling(self.c_triangulation)
        self._cache = {}

    def curve_info(self, max_segments=6):
        """
        curve_info(max_segments=6)

        Lists the dual curves with at most max_segments, showing
        their length and parity.
        """
        dicts = self.curve_info_dicts(max_segments)
        i = 0
        for dict in dicts:
            print '%3d: %s curve of length %s'%(
                i,
                MatrixParity[dict['parity']],
                dict['filled length'])
            i = i+1

    def curve_info_dicts(self, max_segments=6):
        cdef int i, num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_MatrixParity parity
        cdef Complex complete_length, filled_length

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
            info_dict['parity'] = parity
            info_dict['filled length'] = C2C(filled_length)
            info_dict['complete length'] = C2C(complete_length)
            result.append(info_dict)
        free_dual_curves(num_curves, curve_list)
        return result

    def drill(self, which_curve, max_segments=6):
        """
        Drills out the specified dual curve from among all dual
        curves with at most max_segments.
        """
        
        cdef int num_curves
        cdef DualOneSkeletonCurve **curve_list
        cdef c_Triangulation *c_triangulation
        cdef Triangulation result
        cdef char* c_new_name

        new_name = self.name()+'-%d'%which_curve
        c_new_name = new_name

        dual_curves(self.c_triangulation,
                    max_segments,
                    &num_curves,
                    &curve_list)

        if which_curve not in range(num_curves):
            raise IndexError, "Drilling curve requested in not in range(%d)." % num_curves
        
        c_triangulation = drill_cusp(self.c_triangulation,
                                     curve_list[which_curve],
                                     c_new_name)
        free_dual_curves(num_curves, curve_list)

        if c_triangulation == NULL:
            raise RuntimeError, """
        Curve is not isotopic to a geodesic."""
        else:
            result = Manifold()
            result.set_c_triangulation(c_triangulation)
            return result

    def is_isometric_to(self, Manifold other):
        cdef Boolean are_isometric
        cdef c_FuncResult result

        result = compute_isometries(self.c_triangulation, other.c_triangulation, 
                                       &are_isometric, NULL, NULL)
        if FuncResult[result] == 'func_bad_input':
            raise ValueError, """
        Dehn filling coefficients must be relatively prime integers."""

        if FuncResult[result] == 'func_failed':
            raise RuntimeError, """
        SnapPea failed to determine whether the manifolds are isometric."""

        return bool(are_isometric)

# Conversion functions Manifold <-> Triangulation

def Manifold_from_Triangulation(Triangulation T, recompute=True):
    cdef c_Triangulation *c_triangulation
    cdef Manifold M

    copy_triangulation(T.c_triangulation, &c_triangulation)
    M = Manifold()
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
    T = Triangulation()
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
        copy_triangulation(triangulation.c_triangulation, &self.c_triangulation)
        self.c_group_presentation = fundamental_group(
            self.c_triangulation,
            simplify_presentation,
            fillings_may_affect_generators,
            minimize_number_of_generators)
        self.num_cusps = triangulation.num_cusps

    def __dealloc__(self):
        free_triangulation(self.c_triangulation)
        free_group_presentation(self.c_group_presentation)

    def __repr__(self):
        return 'Generators:\n   %s\nRelators:\n   %s'%(
            ','.join(self.generators()),
            '\n   '.join(self.relators()))

    def _word_as_list(self, word):
        if type(word) != types.StringType:
            raise TypeError, """
        Words are represented as Python strings."""
        word_list = []
        generators = self.generators()
        for letter in word:
            try:
                if letter.islower():
                    word_list.append(1 + generators.index(letter))
                else:
                    word_list.append(-1 - generators.index(letter.lower()))
            except ValueError:
                raise RuntimeError, """
        Word contains a non-generator."""
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

    def meridian(self, int which_cusp):
        """
        Returns a word representing a conjugate of the current meridian for
        the given cusp.  Guaranteed to commute with the longitude for the same
        cusp.
        """
        return self.c_word_as_string(
            fg_get_meridian(self.c_group_presentation, which_cusp))

    def longitude(self, int which_cusp):
        """
        Returns a word representing a conjugate of the current
        longitude for the given cusp.  Guaranteed to commute with the
        meridian for the same cusp.  Note: for Klein bottle cusps,
        longitude must be defined carefully.
        """
        return self.c_word_as_string(
            fg_get_longitude(self.c_group_presentation, which_cusp))

    def peripheral_curves(self):
        """
        Returns a list of meridian-longitude pairs for all cusps.
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

    def _magma_init_(self):
        return self.magma_string()

class FundamentalGroup(CFundamentalGroup):
    """
    A FundamentalGroup represents a presentation of the fundamental
    group of a SnapPea Triangulation.  Group elements are described as
    words in the generators a,b,..., where the inverse of a is denoted
    A.  Words are represented by python strings (and the concatenation
    operator is named"+", according to Python conventions).

    Instantiate as FundamentalGroup(T), where T is a Triangulation.

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

cdef class CHolonomyGroup(CFundamentalGroup):
    def _matrices(self, word):
        """
        Returns (M,O) where M = SL2C(word) and O = O31(word).
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
            return sl2, o31
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

class HolonomyGroup(CHolonomyGroup):
    """
    A HolonomyGroup is a FundamentalGroup with added structure
    consisting of a holonomy representation into O(3,1), and an
    arbitrarily chosen lift of the holonomy representation to SL(2,C).
    The holonomy is determined by the shapes of the tetrahedra, so a
    HolonomyGroup is associated to a Manifold, while a Triangulation
    only has a FundamentalGroup.  Methods are provided to evaluate the
    representations on a group element.
    
    Instantiate as HolonomyGroup(M), where M is a Manifold.

    Methods (in addition to those inherited from FundamentalGroup):
        O31(word)        --> evaluates the holonomy of the word
        SL2C(word)       --> evaluates the chosen lift of the holonomy
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
            raise RuntimeError, 'Dirichet construction failed.'
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
        keys 'vertices', 'distance, 'closest', hue.  The distance
        from the origin the value for 'distance', and the value for
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
            raise RuntimeError, 'Please install PyOpenGL and numpy to use this feature.'

class DirichletDomain(CDirichletDomain):
    """
    A DirichletDomain object represents a Dirichlet Domain of 
    a hyperbolic manifold, typically centered at a point which
    is a local maximum of injectivity radius.  It will have ideal
    vertices if the manifold is not closed.

    Instantiate as DirichletDomain(M) where M is a Manifold to
    obtain a Dirichlet Domain centered at a point which maximizes
    injectivity radius.

    Other options can be provided to customize the computation:
    DirichletDomain(M,
                    vertex_epsilon=10.0**-8,
                    displacement = [0.0, 0.0, 0.0],
                    centroid_at_origin=True,
                    maximize_injectivity_radius=True)

    Methods:

    """
    pass

# Cusp Neighborhoods

cdef class CCuspNeighborhood:
    cdef CuspNeighborhoods *c_cusp_neighborhood
    cdef c_Triangulation *c_triangulation

    def __new__(self, Manifold manifold):
        copy_triangulation(manifold.c_triangulation, &self.c_triangulation)
        self.c_cusp_neighborhood = initialize_cusp_neighborhoods(
            self.c_triangulation)
        if self.c_cusp_neighborhood == NULL:
            raise RuntimeError, 'Cusp Neighborhood construction failed.'
        self.manifold_name = manifold.name()

    def __dealloc__(self):
        if self.c_triangulation != NULL:
            free_triangulation(self.c_triangulation)
        if self.c_cusp_neighborhood != NULL:
            free_cusp_neighborhoods(self.c_cusp_neighborhood)

    def __repr__(self):
        N = self.num_components()
        return 'Cusp Neighborhood with %d component%s'%(
            N, N != 1 and 's' or '')

    def check_index(self, cusp_index):
        N = int(cusp_index)
        if 0 <= N <= self.num_components():
            return N
        else:
            raise ValueError, 'Invalid cusp index %s. There are %d cusps'%(
                cusp_index, self.num_components()) 

    def num_components(self):
        """
        Return the number of cusps.
        """
        return get_num_cusp_neighborhoods(self.c_cusp_neighborhood)

    def topology(self, cusp_index):
        """
        Return the topological type of the specified cusp.
        """
        N = self.check_index(cusp_index)
        topology = get_cusp_neighborhood_topology(self.c_cusp_neighborhood,N)
        return CuspTopology[topology]

    def get_displacement(self, cusp_index):
        """
        Return the displacement of the specified cusp.
        """
        N = self.check_index(cusp_index)
        return get_cusp_neighborhood_displacement(self.c_cusp_neighborhood,N)

    def set_displacement(self, cusp_index, new_displacement):
        """
        Set the displacement of the specified cusp.
        """
        N = self.check_index(cusp_index)
        set_cusp_neighborhood_displacement(self.c_cusp_neighborhood,
                                           N,
                                           new_displacement)

    def get_tie(self, cusp_index):
        """
        Return True if the specified cusp is a member of the tied group. 
        """
        N = self.check_index(cusp_index)
        return get_cusp_neighborhood_tie(self.c_cusp_neighborhood,N)

    def set_tie(self, cusp_index, new_tie):
        """
        Mark the specified cusp as a member of the tied group. 
        """
        N = self.check_index(cusp_index)
        set_cusp_neighborhood_tie(self.c_cusp_neighborhood, N, new_tie)

#    extern double get_cusp_neighborhood_cusp_volume(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern double get_cusp_neighborhood_manifold_volume(CuspNeighborhoods *cusp_neighborhoods)
#    extern c_Triangulation *get_cusp_neighborhood_manifold(CuspNeighborhoods *cusp_neighborhoods)
#    extern double get_cusp_neighborhood_reach(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern double get_cusp_neighborhood_max_reach(CuspNeighborhoods *cusp_neighborhoods)
#    extern double get_cusp_neighborhood_stopping_displacement(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern int get_cusp_neighborhood_stopper_cusp_index(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern void get_cusp_neighborhood_translations(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Complex *meridian, Complex *longitude)
#    extern CuspNbhdHoroballList *get_cusp_neighborhood_horoballs(CuspNeighborhoods *cusp_neighborhoods, int cusp_index, Boolean full_list, double cutoff_height)
#    extern void free_cusp_neighborhood_horoball_list(CuspNbhdHoroballList *horoball_list)
#    extern CuspNbhdSegmentList *get_cusp_neighborhood_triangulation(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern CuspNbhdSegmentList *get_cusp_neighborhood_Ford_domain(CuspNeighborhoods *cusp_neighborhoods, int cusp_index)
#    extern void free_cusp_neighborhood_segment_list(CuspNbhdSegmentList *segment_list)


class CuspNeighborhood(CCuspNeighborhood):
    """
    A CuspNeighborhood object represents an equivariant collection of disjoint
    horoballs that project to cusp neighborhoods.
 
    Methods:
      num_components()   How many cusps are there?
    """
    pass

# get_triangulation

split_filling_info = re.compile("(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))+)")
is_census_manifold = re.compile("([msvxy])([0-9]+)$")
is_torus_bundle = re.compile("b([+-no])([+-])([lLrR]+)$")
is_knot_complement = re.compile("(?P<crossings>[0-9]+)_(?P<index>[0-9]+)$")
is_link_complement1 = re.compile("(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$")
is_link_complement2 = re.compile("(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$")
is_link_complement3 = re.compile("[lL]([0-9]+)")
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)')
is_braid_complement = re.compile("braid(\[[0-9, -]+\])")

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


# The below is copied by hand into the Triangulation and Manifold classes above

triangulation_help =  """
    A %s is specified by a string, according to the
    following conventions:

    1. Numbers in parens at the end specify Dehn fillings.  For example
    'm125(1,2)(4,5)' means do (1,2) filling on the first cusp and (4,5)
    filling on the second cusp of the census manifold m125.

    2. Strings of the form 'm123', 's123', 'v123', and so on refer to the
    SnapPea Cusped Census manifolds.

    3. Strings of the form '4_1', '04_1', '4_01', '5^2_6', '6_4^7',
    etc, refer to complements of links in Rolfsen's table.  Similarly
    for 'L20935', 'l104001', etc.

    4. Strings of the form 'b++LLR', 'b+-llR', 'bo-RRL', 'bn+LRLR'
    refer to the correponding punctured torus bundle.

    5. Strings of the form '11a17' or '12n345' refer to complements of
    knots in the Hoste-Thistlethwaite tables.

    6. Strings of the form 'braid[1,2,-3,4]' creates fibered manifold
    corresponding to the given braid.  In other words, think of a braid
    as givening an element of the mapping class group of the
    numStrands-punctured disc.  This function returns the
    corresponding mapping torus.  If you want the braid closure, you
    have to do (1,0) filling of the last cusp.

    If the string is not in any of the above forms it is assumed to be
    the name of a SnapPea manifold file.  The file will be loaded
    if found in the current directory or the path given by the user
    variable SNAPPEA_MANIFOLD_DIRECTORY.
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
        pathname =  os.path.join(link_directory, filename)
        if os.path.isfile(pathname):
            c_triangulation = read_triangulation(pathname)
        else:
            raise IOError, """
        The link complement %s was not found."""%real_name
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

    # Step 5. See if a (fibered) braid compelement is requested

    m = is_braid_complement.match(real_name)
    if m:
        word = eval(m.group(1))
        num_strands = max([abs(x) for x in word]) + 1
        c_triangulation = get_fibered_manifold_associated_to_braid(num_strands, word)
        set_cusps(c_triangulation, fillings)
        return c_triangulation

    # Step 6. If all else fails, try to load a manifold from a file.
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
        
    # Step 7. Give up.
    raise IOError, """
        The manifold file %s was not found.  Sorry.\n%s"""%(
        real_name, triangulation_help%'Triangulation or Manifold')
        
cdef int set_cusps(c_Triangulation* c_triangulation, fillings) except -1:
    if c_triangulation == NULL:
        return 0
    if len(fillings) > 0:
        num_cusps = get_num_cusps(c_triangulation) 
        if len(fillings) > num_cusps:
            raise ValueError, """
        The number of fillings specified exceeds the number of cusps."""
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

# These are the gzipped files holding the knot tables.
Alternating_table = gzip.open(os.path.join(table_directory, 'alternating.gz') )
Nonalternating_table = gzip.open(os.path.join(table_directory, 'nonalternating.gz') )

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

cdef c_Triangulation* get_HT_knot(crossings, alternation, index) except ? NULL:
    cdef int* DT_array
    cdef int i
    cdef c_Triangulation* c_triangulation
    DT = get_HT_knot_DT(crossings, alternation, index)
    DT_array = <int*>malloc(len(DT)*sizeof(int))
    for i from 0 <= i < len(DT):
        DT_array[i] = DT[i]
    c_triangulation = DT_int_to_triangulation(len(DT), DT_array)
    name = "%d" % crossings + alternation + "%d" % index
    set_triangulation_name(c_triangulation, name)
    free(DT_array)
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
            raise IndexError, """
        Index out of range."""
        c_triangulation = GetCuspedCensusManifold(
            manifold_path, num_tet, self.orientability, census_index)
        if c_triangulation == NULL:
            print num_tet, census_index
            raise RuntimeError, """
        SnapPea failed to read census manifold."""
        result = Manifold(spec='empty')
        result.set_c_triangulation(c_triangulation)
        return result

class OrientableCuspedCensus(CuspedCensus):
    """
    Iterator/Sequence for orientable manifolds in the SnapPea
    Cusped Census.
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

try:
    prompt = sys.ps1
    print """
    Hi.  I'm SnapPy.  
    SnapPy is based on the SnapPea kernel, written by Jeff Weeks.
    """

    if prompt.startswith('>>>'):
        print "Type doc() for help, or doc(X) for help on X."
except:
    pass
