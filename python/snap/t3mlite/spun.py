"""
Spun-normal surfaces in cusped manifolds.

Despite living in the t3mlite directory, it is essentially independent
of it (though depending very much on snappy), except for borrowing
some linear algebra code.
"""

import snappy, FXrays
if snappy._within_sage:
    from sage.all import gcd
    from sage.all import vector as Vector
    from sage.all import matrix as Matrix
else:
    from snappy.snap.t3mlite.linalg import Vector, Matrix, gcd

def weak_normalize_slope(slope):
    """For a tuple (a, b), scale it so that gcd(a,b)=1"""
    a, b = [int(s) for s in slope]
    if a == b == 0:
        return (0, 0)
    g = gcd(a,b)
    a, b = a//g, b//g
    return (a,b)
    
def normalize_slope(slope):
    """
    For a tuple (a, b), scale it so that gcd(a,b)=1 and it lies in the
    right half plane.

    >>> normalize_slope( (-10, 5) )
    (2, -1)
    >>> normalize_slope( (-3, 0) )
    (1, 0)

    The corner case of (0, b) is handled like this:
    
    >>> normalize_slope( (0, -10) )
    (0, 1)
    """
    a, b = weak_normalize_slope(slope)
    if a == b == 0:
        return (0, 0) 
    if a < 0:
        a, b = -a, -b
    elif a == 0 and b < 0:
        b = -b
    return (a,b)

def shift_matrix(n):
    """
    The edge shifts corresponding to each of three quad types, see Figure
    2.1 and 2.2 of `[DG] <http://arxiv.org/abs/1102.4588>`_
    """
    shifts = Matrix(3*n, 3*n)
    for i in range(0, 3*n, 3):
        for j in range(3):
            shifts[i+j, i+j] = 1
            shifts[i+j, i+((j+1) % 3)] = -1
    return shifts

def quad_vector_to_type_and_coeffs(quad_vector):
    """
    For an n-tetrahedra manifold, take a full quad vector
    of length 3n and store the quad type and weight for
    each tetrahedron.  
    """
    quad_types, coefficients = [], []
    quad_vector = list(quad_vector)
    for i in range(len(quad_vector)//3):
        one_tet = quad_vector[3*i:3*(i+1)]
        pos = [ (i, c) for i, c in enumerate(one_tet) if c > 0]
        assert len(pos) <= 1
        if pos:
            q, c = pos[0]
            quad_types.append(q)
            coefficients.append(c)
        else:
            quad_types.append(None)
            coefficients.append(0)

    return quad_types, Vector(coefficients)

class SpunSurface:
    """
    A spun normal surface in an ideal triangulation, as introduced by
    Thurston.

    For an quick sketch of this theory see `[DG]
    <http://arxiv.org/abs/1102.4588>`_ and for more details see
    `[Tillmann] <http://arxiv.org/abs/math/0406271>`_.

    The quad conventions are (Q02, Q03, Q01) corresponding to
    z -> 0, z' -> 0, and z'' -> 0 respectively, as per Figure 3.1 of
    `[DG] <http://arxiv.org/abs/1102.4588>`_.  The quad types
    are numbered 0, 1, 2; the "None" quad type means a
    tetrahedron contains no quads at all.  
    """ 
    def __init__(self, manifold, quad_vector=None, quad_types=None, index=None):
        self._manifold = manifold
        self._index=index
        if quad_types != None:
            coefficients = quad_vector
            quad_vector = []
            for c, q in zip(coefficients, quad_types):
                three = [0, 0, 0]
                three[q] = c
                quad_vector += three
        quad_vector = Vector(quad_vector)
        eqns = manifold._normal_surface_equations()
        assert eqns.is_solution(quad_vector)
        self._quad_vector = quad_vector
        self._quad_types, self._coefficients = quad_vector_to_type_and_coeffs(quad_vector)
        self._boundary_slopes = eqns.boundary_slope_of_solution(quad_vector)

    def quad_vector(self):
        return self._quad_vector

    def quad_types(self):
        return self._quad_types

    def coefficients(self):
        return self._coefficients

    def boundary_slopes(self):
        return self._boundary_slopes
    
    def is_compatible(self, other):
        for a, b in zip(self._quad_types, other._quad_types):
            if not (a == b or None in (a, b)):
                return False
        return True
                
    def __radd__(self, other):
        if other==0:
            return self
        
    def __add__(self, other):
        if other==0:
            return self
        if not self.is_compatible(other):
            raise ValueError('Normal surfaces are not compatible')
        return SpunSurface(self._manifold, self._quad_vector + other._quad_vector)

    def __repr__(self):
        return "<Surface %s: %s %s %s>" % (self._index, self._quad_types,
                                        list(self._coefficients), tuple(self._boundary_slopes))

class SpunNormalSurfaceEquations:
    def __init__(self, manifold):
        self.manifold = manifold
        n = manifold.num_tetrahedra()
        self.shift_matrix = shift_matrix(n)
        gluing_equations = list(manifold.gluing_equations())
        edge_equations = Matrix(gluing_equations[:n])
        self.quad_equations = edge_equations * self.shift_matrix
        self.cusp_equations = cusp_equations = [Vector(eqn) for eqn in gluing_equations[n:]]
        slope_matrix = []
        for i in range(manifold.num_cusps()):
            slope_matrix += [-cusp_equations[2*i+1], cusp_equations[2*i]]
        self.slope_matrix = Matrix(slope_matrix)

    def vertex_solutions(self, algorithm='FXrays'):
        if algorithm == 'FXrays':
            M = self.quad_equations
            return FXrays.find_Xrays(M.nrows(), M.ncols(), M.list(),
                                     0, print_progress=False)
        elif algorithm == 'regina':
            try:
                import regina
            except ImportError:
                raise ImportError('Regina module not available')
            M = self.manifold
            T = regina.NTriangulation(M._to_string())
            ans = []
            tets = range(M.num_tetrahedra())
            surfaces = regina.NNormalSurfaceList.enumerate(T, regina.NS_QUAD)
            for i in range(surfaces.getNumberOfSurfaces()):
                S = surfaces.getSurface(i)
                coeff_vector = [int(S.getQuadCoord(tet, quad).stringValue())
                                for tet in tets for quad in (1, 2, 0)]
                ans.append(coeff_vector)
            return ans
        else:
            raise ValueError("Algorithm should be one of {'FXrays', 'regina'}")


    def is_solution(self, quad_vector):
        return self.quad_equations * quad_vector == 0

    def boundary_slope_of_solution(self, quad_vector):
        return self.slope_matrix*self.shift_matrix*quad_vector



# The following methods get monkey patched into the manifold
# classes.  

def _normal_surface_equations(self):
    name = '_normal_surface_equations'
    if name not in self._cache:
        eqns = SpunNormalSurfaceEquations(self)
        self._cache[name] = SpunNormalSurfaceEquations(self)
    return self._cache[name] 

def normal_surfaces(self, algorithm='FXrays'):
    """
    All the vertex spun-normal surfaces in the current triangulation.

    >>> M = Manifold('m004')
    >>> M.normal_surfaces()    # doctest: +NORMALIZE_WHITESPACE
    [<Surface 0: [0, 0] [1, 2] (4, 1)>,
     <Surface 1: [0, 1] [1, 2] (4, -1)>,
     <Surface 2: [1, 2] [2, 1] (-4, -1)>,
     <Surface 3: [2, 2] [2, 1] (-4, 1)>]
    """
    if 'normal_surfaces' not in self._cache:
        eqns = self._normal_surface_equations()
        self._cache['normal_surfaces'] = [SpunSurface(self, qv, index=i)
                            for i, qv in enumerate(eqns.vertex_solutions(algorithm))]
    return self._cache['normal_surfaces']

def normal_boundary_slopes(self, subset='all', algorithm='FXrays'):
    """
    For a one-cusped manifold, returns all the nonempty boundary slopes of
    spun normal surfaces.  Provided the triangulation supports a
    genuine hyperbolic structure, then by `Thurston and Walsh
    <http://arxiv.org/abs/math/0503027>`_ any strict boundary slope
    (the boundary of an essential surface which is not a fiber or
    semifiber) must be listed here.

    >>> M = Manifold('K3_1')
    >>> M.normal_boundary_slopes()
    [(16, -1), (20, -1), (37, -2)]

    If the ``subset`` flag is set to ``'kabaya'``, then it only
    returns boundary slopes associated to vertex surfaces with a quad
    in every tetrahedron; by Theorem 1.1. of `[DG]
    <http://arxiv.org/abs/1102.4588>`_ these are all strict boundary
    slopes.

    >>> N = Manifold('m113')
    >>> N.normal_boundary_slopes()
    [(1, 1), (1, 2), (2, -1), (2, 3), (8, 11)]
    >>> N.normal_boundary_slopes('kabaya')
    [(8, 11)]

    If the ``subset`` flag is set to ``'brasile'`` then it returns
    only the boundary slopes that are associated to vertex surfaces
    giving isolated rays in the space of embedded normal surfaces.

    >>> N.normal_boundary_slopes('brasile')
    [(1, 2), (8, 11)]
    """
    if not self.is_orientable():
        raise ValueError('Manifold must be orientable')
    if self.num_cusps() != 1:
        raise ValueError('More than 1 cusp, so need to look at the surfaces directly.')

    surfaces = self.normal_surfaces(algorithm)
    if subset is 'kabaya':
        surfaces = [S for S in surfaces if min(S.coefficients()) > 0]
    elif subset is 'brasile':
        isolated_surfaces = []
        for S in surfaces:
            isolated = True
            for F in surfaces:
                if S != F and S.is_compatible(F):
                    isolated = False
                    break
            if isolated:
                isolated_surfaces.append(S)
        surfaces = isolated_surfaces
    else:
        if subset != 'all':
            raise ValueError("Subset must be one of 'all', 'kabaya', and 'brasile'")
    
    slopes = set([normalize_slope(S.boundary_slopes()) for S in surfaces])
    slopes.discard( (0, 0) )
    return sorted(slopes)
        
if __name__ == "__main__":
    import doctest
    names = {'Manifold':snappy.Manifold}
    doctest.testmod(extraglobs=names)
