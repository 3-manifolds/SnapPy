# $Id: mcomplex.py,v 1.14 2009/08/20 15:58:58 t3m Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from .simplex import *
from .tetrahedron import Tetrahedron
from .corner import Corner
from .arrow import Arrow
from .face import Face
from .edge import Edge
from .vertex import Vertex
from .surface import Surface, SpunSurface, ClosedSurface, ClosedSurfaceInCusped
from .perm4 import Perm4, inv
from . import files
from . import linalg
from . import homology
import sys
import random
import io

try:
    import snappy
except ImportError:
    snappy = None

VERBOSE = 0

# Globals needed for normal surfaces:

# The height shift dictionaries for the three quad types.
# Shift[ tet edge ] is a tuple of the shifts of the three quad
# types (Q03, Q13, Q23) along that edge.
#
# That is the first entry E01:(-1, 1, 0) means that Q03 shifts
# by -1 along E01, Q13 shifts by +1 along E01 and Q23 doesn't
# shift.

Shift = {E01:(-1,1,0), E02:(1,0,-1), E21:(0,-1,1),
         E32:(-1,1,0), E31:(1,0,-1), E03:(0,-1,1)}

# The solution vector templates for the four vertex types.

VertexVector = {V0:(1,0,0,0), V1:(0,1,0,0),
                V2:(0,0,1,0), V3:(0,0,0,1)}


def edge_and_arrow(edge_or_arrow):
    """
    Given and edge or an arrow, returns the corresponding compatible
    (edge, arrow) pair.
    """
    if isinstance(edge_or_arrow, Edge):
        edge = edge_or_arrow
        arrow = Arrow(edge.Corners[0].Subsimplex,
                  LeftFace[edge.Corners[0].Subsimplex],
                  edge.Corners[0].Tetrahedron)
    else:
        if not isinstance(edge_or_arrow, Arrow):
            raise ValueError('Input edge_or_arrow is neither')
        arrow = edge_or_arrow.copy()
        edge = arrow.axis()
    return edge, arrow


class Insanity(Exception):
    pass


class Mcomplex:
    """
    An Mcomplex is a union of tetrahedra with faces identified in
    pairs.  The edges (vertices) are equivalence classes under the
    induced equivalence relation on the set of edges (vertices) of the
    tetrahedra.

    >>> T = Mcomplex([Tetrahedron()])
    >>> len(T), len(T.Vertices)
    (1, 4)
    >>> T = Mcomplex('m004')
    >>> len(T)
    2
    >>> tet_data = [([0,1,0,1], [(2,1,0,3), (0,3,2,1), (2,1,0,3), (0,1,3,2)]),
    ...             ([1,1,0,0], [(1,0,2,3), (1,0,2,3), (0,1,3,2), (0,3,2,1)])]
    >>> S = Mcomplex(tet_data)
    >>> len(S)
    2
    """
    def __init__(self, tetrahedron_list=None):
        if tetrahedron_list is None:
            tetrahedron_list = []
        elif isinstance(tetrahedron_list, str) and snappy is None:
            tetrahedron_list = tets_from_data(files.read_SnapPea_file(file_name=tetrahedron_list))
        elif snappy:
            if isinstance(tetrahedron_list, str):
                tetrahedron_list = snappy.Triangulation(tetrahedron_list,
                                                        remove_finite_vertices=False)
            if hasattr(tetrahedron_list, '_get_tetrahedra_gluing_data'):
                tetrahedron_list = tets_from_data(
                    tetrahedron_list._get_tetrahedra_gluing_data())
        if isinstance(tetrahedron_list, (list, tuple)):
            if len(tetrahedron_list) > 0 and not isinstance(tetrahedron_list[0], Tetrahedron):
                tetrahedron_list = tets_from_data(tetrahedron_list)

        self.Tetrahedra = tetrahedron_list
        self.Edges = []
        self.Faces = []
        self.Vertices = []
        self.NormalSurfaces = []
        self.AlmostNormalSurfaces = []
        self.build()

    def copy(self, base_arrow=None):
        new_tets = []
        new_to_old = {}
        old_to_new = {}
        for tet in self.Tetrahedra:
            new_tet = Tetrahedron()
            old_to_new[tet] = new_tet
            new_to_old[new_tet] = tet
            new_tets.append(new_tet)
        for new_tet in new_tets:
            for face in TwoSubsimplices:
                new_tet.attach(face,
                               old_to_new[new_to_old[new_tet].Neighbor[face]],
                               new_to_old[new_tet].Gluing[face].tuple())
        if base_arrow is None:
            return self.__class__(new_tets)
        else:
            new_arrow = base_arrow.copy()
            new_arrow.Tetrahedron = old_to_new[base_arrow.Tetrahedron]
            return (self.__class__(new_tets), new_arrow)

    def build(self):
        for i in range(len(self.Tetrahedra)):
            self.Tetrahedra[i].Index = i
        self.build_face_classes()
        self.build_edge_classes()
        self.build_vertex_classes()
        self.build_one_skeleton()
        self.LinkGenera = [vertex.link_genus() for vertex in self.Vertices]

    def rebuild(self):
        for tet in self.Tetrahedra:
            tet.clear_Class()
        for face in self.Faces:
            face.erase()
        for edge in self.Edges:
            edge.erase()
        for vertex in self.Vertices:
            vertex.erase()
        self.Faces = []
        self.Edges = []
        self.Vertices = []
        self.build()

    def add_tet(self, tet):
        self.Tetrahedra.append(tet)

    def clear_tet(self,tet):
        """
        Remove the face, edge and vertex classes of a tetrahedron.
        This should destroy the faces, edges and vertices that meet
        the tetrahedron.  A call to build_face_classes,
        build_edge_classes or build_vertex_classes will then rebuild
        the neighborhood without having to rebuild the whole manifold.
        """
        for two_subsimplex in TwoSubsimplices:
            face = tet.Class[two_subsimplex]
            if face is not None:
                face.erase()
            try:
                self.Faces.remove(face)
            except ValueError:
                pass

        for one_subsimplex in OneSubsimplices:
            edge = tet.Class[one_subsimplex]
            if edge is not None:
                edge.erase()
            try:
                self.Edges.remove(edge)
            except ValueError:
                pass

        for zero_subsimplex in ZeroSubsimplices:
            vertex = tet.Class[zero_subsimplex]
            if vertex is not None:
                vertex.erase()
            try:
                self.Vertices.remove(vertex)
            except ValueError:
                pass

    def delete_tet(self, tet):
        """
        Clear a tetrahedron, then remove it from the Tetrahedron list.
        """
        self.clear_tet(tet)
        tet.erase()
        self.Tetrahedra.remove(tet)

    def new_arrow(self):
        """
        Add one new tetrahedron and return one of its arrows.
        """
        tet = Tetrahedron()
        self.add_tet(tet)
        return Arrow(E01,F3,tet)

    def new_arrows(self,n):
        return [self.new_arrow() for i in range(n)]

    def new_tet(self):
        tet = Tetrahedron()
        self.add_tet(tet)
        return tet

    def new_tets(self,n):
        return [self.new_tet() for i in range(n)]

    def _triangulation_data(self):
        ans = []
        # We don't assume that the indices of the Tetraheda are equal
        # to range(len(self))
        tet_to_index = {T:i for i, T in enumerate(self.Tetrahedra)}
        for T in self.Tetrahedra:
            neighbors, perms = [], []
            for v in TwoSubsimplices:
                if T.Neighbor[v] is None:
                    neighbor, perm = None, None
                else:
                    neighbor = tet_to_index[T.Neighbor[v]]
                    perm = T.Gluing[v].tuple()
                neighbors.append(neighbor)
                perms.append(perm)
            ans.append((neighbors, perms))
        return ans

    def __len__(self):
        """
        Return the number of tetrahedra
        """
        return len(self.Tetrahedra)

    def __getitem__(self, index):
        """
        M[i] refers to the ith Tetrahedron of the mcomplex M.
        """
        return self.Tetrahedra[index]

    def info(self, out=sys.stdout):
        """
        M.info() describes the Mcomplex.
        """
        try:
            out.write( "Mcomplex with %d Tetrahedra\n\n" % len(self) )
            for tet in self.Tetrahedra:
                tet.info(out)
            out.write("\nEdges:\n")
            for edge in self.Edges:
                edge.info(out)
        except IOError:
            pass

    def build_edge_classes(self):
        """
        Construct the edge classes and compute valences.
        """
        for tet in self.Tetrahedra:
            for one_subsimplex in OneSubsimplices:
                if ( tet.Class[one_subsimplex] is None ):
                    newEdge = Edge()
                    self.Edges.append(newEdge)
                    first_arrow = Arrow(one_subsimplex, RightFace[one_subsimplex], tet)
                    a = first_arrow.copy()
                    sanity_check = 0
                    boundary_hits = 0
                    while 1:
                        # Walk around the edge.
                        if sanity_check > 6*len(self.Tetrahedra):
                            raise Insanity('Bad gluing data: could not construct edge link.')
                        # Record the corners and edge classes as we go.
                        newEdge._add_corner(a)
                        a.Tetrahedron.Class[a.Edge] = newEdge
                        if a.next() is None:
                            # We hit the boundary!
                            # Go back to the beginning and walk to the right.
                            # If this is our second boundary hit, we are done.
                            if not boundary_hits == 0:
                                newEdge.RightBdryArrow = a.copy()
                                # Make sure the corner list goes from right to left.
                                newEdge.Corners.reverse()
                                break
                            else:
                                boundary_hits = 1
                                newEdge.LeftBdryArrow = a.copy()
                                newEdge.IntOrBdry = 'bdry'
                                a = first_arrow.copy()
                                a.reverse()
                                # Don't record the first corner twice.
                                del newEdge.Corners[0]
                                # Reverse the corner list since we are now going right.
                                newEdge.Corners.reverse()
                        else:
                            # Stop if we get back to where we started.
                            if a == first_arrow:
                                newEdge.IntOrBdry = 'int'
                                break
                        sanity_check = sanity_check + 1
        self.EdgeValences = [edge.valence() for edge in self.Edges]
        for i in range(len(self.Edges)):
            self.Edges[i].Index = i

    def build_vertex_classes(self):
        for tet in self.Tetrahedra:
            for zero_subsimplex in ZeroSubsimplices:
                if ( tet.Class[zero_subsimplex] is None ):
                    newVertex = Vertex()
                    self.Vertices.append(newVertex)
                    self.walk_vertex(newVertex,zero_subsimplex,tet)
        for i in range(len(self.Vertices)):
            self.Vertices[i].Index = i

    def walk_vertex(self,vertex,zero_subsimplex,tet):
        if (tet.Class[zero_subsimplex] is not None ):
            return
        else:
            tet.Class[zero_subsimplex] = vertex
            vertex.Corners.append(Corner(tet,zero_subsimplex))
            for two_subsimplex in TwoSubsimplices:
                if ( is_subset(zero_subsimplex,two_subsimplex)
                     and
                     tet.Gluing[two_subsimplex] is not None):
                    self.walk_vertex(vertex,
                                     tet.Gluing[two_subsimplex].image(zero_subsimplex),
                                     tet.Neighbor[two_subsimplex])

    def build_one_skeleton(self):
        """
        Construct the 1-skeleton, i.e. record which edges are
        connected to which vertices.  This assumes that Edges and Vertices
        have already been built.
        """
        for edge in self.Edges:
            tet = edge.Corners[0].Tetrahedron
            one_subsimplex = edge.Corners[0].Subsimplex
            tail = tet.Class[Tail[one_subsimplex]]
            head = tet.Class[Head[one_subsimplex]]
            edge.Vertices = [tail , head]
            tail.Edges.append(edge)
            head.Edges.append(edge)
            if edge.IntOrBdry == 'bdry':
                tail.IntOrBdry = 'bdry'
                head.IntOrBdry = 'bdry'
        for vertex in self.Vertices:
            if vertex.IntOrBdry == '':
                vertex.IntOrBdry = 'int'

    def build_face_classes(self):
        """
        Construct the faces.
        """
        for tet in self.Tetrahedra:
            for two_subsimplex in TwoSubsimplices:
                if ( tet.Class[two_subsimplex] is None ):
                    newFace = Face()
                    self.Faces.append(newFace)
                    newFace.Corners.append(Corner(tet,two_subsimplex))
                    tet.Class[two_subsimplex] = newFace
                    othertet = tet.Neighbor[two_subsimplex]
                    if othertet:
                        newFace.IntOrBdry = 'int'
                        othersubsimplex = tet.Gluing[two_subsimplex].image(two_subsimplex)
                        newFace.Corners.append(Corner(othertet, othersubsimplex))
                        othertet.Class[othersubsimplex] = newFace
                    else:
                        newFace.IntOrBdry = 'bdry'
        for i in range(len(self.Faces)):
            self.Faces[i].Index = i

    def orient(self):
        """
        The simplification moves below assume that the Mcomplex is oriented.
        Yes, oriented, not just orientable.  An Mcomplex has been oriented if
        all of the gluing permutations are odd.  The orient method walks through
        the manifold reorienting tetrahedra to try to get all of the gluing
        permutations to be odd.  Returns True on success, False if the manifold is
        not orientable.
        """
        for tet in self.Tetrahedra:
            tet.Checked = 0
        self.walk_and_orient(self[0], 1)
        self.rebuild()
        return self.is_oriented()

    def is_oriented(self):
        for tet in self.Tetrahedra:
            for two_subsimplex in TwoSubsimplices:
                if (not tet.Neighbor[two_subsimplex] is None
                        and tet.Gluing[two_subsimplex].sign() == 0):
                    return False
        return True

    def walk_and_orient(self, tet, sign):
        if tet.Checked == 1:
            return
        tet.Checked = 1
        if sign == 0:
            tet.reverse()
        for ssimp in TwoSubsimplices:
            if tet.Neighbor[ssimp] is not None:
                self.walk_and_orient(tet.Neighbor[ssimp], tet.Gluing[ssimp].sign())

    def build_matrix(self):
        """
        Convention is that the ordered quads are (Q03, Q13, Q23).
        """
        int_edges = [edge for edge in self.Edges if edge.IntOrBdry == 'int']
        self.QuadMatrix = linalg.Matrix(len(int_edges), 3*len(self))
        for edge in int_edges:
            for corner in edge.Corners:
                i = int_edges.index(edge)
                j = corner.Tetrahedron.Index
                for k in range(3):
                    self.QuadMatrix[i,3*j+k] += Shift[corner.Subsimplex][k]
        self.build_vertex_incidences()

    def build_vertex_incidences(self):
        for vertex in self.Vertices:
            vertex.IncidenceVector = linalg.Vector( 4*len(self) )
            for corner in vertex.Corners:
                j = corner.Tetrahedron.Index
                vertex.IncidenceVector[4*j:4*j+4] += VertexVector[corner.Subsimplex]

    def find_normal_surfaces(self, modp=0, print_progress=False,
                             algorithm='FXrays'):
        """
        Convention is that the ordered quads are (Q03, Q13, Q23).
        """
        self.NormalSurfaces = []
        self.build_matrix()
        if algorithm == 'FXrays':
            try:
                import FXrays
            except ImportError:
                raise ImportError("You need to install the FXrays module"
                                  "if you want to find normal surfaces.")
            coeff_list = FXrays.find_Xrays(self.QuadMatrix.nrows(),
                                           self.QuadMatrix.ncols(),
                                           self.QuadMatrix.entries(), modp,
                                           print_progress=print_progress)

        elif algorithm == 'regina':
            T = self.regina_triangulation()
            import regina
            coeff_list = []
            tets = range(len(self))
            surfaces = regina.NNormalSurfaceList.enumerate(T, regina.NS_QUAD)
            for i in range(surfaces.getNumberOfSurfaces()):
                S = surfaces.getSurface(i)
                coeff_vector = [int(S.getQuadCoord(tet, quad).stringValue())
                                for tet in tets for quad in (2, 1, 0)]
                coeff_list.append(coeff_vector)

        else:
            raise ValueError("Algorithm must be in {'FXrays', 'regina'}")

        for coeff_vector in coeff_list:
            if max(self.LinkGenera) == 0:
                self.NormalSurfaces.append(ClosedSurface(self, coeff_vector))
            elif self.LinkGenera.count(1) == len(self.LinkGenera):
                self.NormalSurfaces.append(SpunSurface(self, coeff_vector))
            else:
                self.NormalSurfaces.append(Surface(self, coeff_vector))

    def normal_surface_info(self, out=sys.stdout):
        try:
            for surface in self.NormalSurfaces:
                out.write("-------------------------------------\n\n")
                surface.info(self, out)
                out.write('\n')
        except IOError:
            pass

    def almost_normal_surface_info(self, out=sys.stdout):
        try:
            for surface in self.AlmostNormalSurfaces:
                out.write("-------------------------------------\n\n")
                surface.info(self, out)
                out.write('\n')
        except IOError:
            pass

    # Simplification Moves
    #
    # The simplification moves require that the list of edge classes
    # be up to date.  Edge classes are recomputed as part of each
    # move.  The vertex classes are not used, nor are they updated, by
    # these moves, with the exception of randomize.

    def _face_permits_two_to_three(self, a, b):
        S, T = a.Tetrahedron, b.Tetrahedron
        if S is None:
            return False, 'Tetrahedron not attached to face'
        if S == T:
            return False, 'Two tetrahedra are the same'
        return True, None

    def _two_to_three_move_hook(self, old_arrow, new_arrows):
        pass

    def two_to_three(self, face_or_arrow, tet=None,
                     return_arrow=False, must_succeed=False,
                     unsafe_mode=False):
        """
        Perform a 2-to-3 Pachner move on the face specified by
        (face_or_arrow, tet), replacing the two tetrahedra with three
        tetrahedra around an edge.

        Returns ``True`` or ``False`` depending on whether the
        requested move succeeded.  When ``must_succeed`` is ``True``,
        it instead raises an exception if the requested move is
        topologically impossible.

        When ``unsafe_mode`` is ``True`` it does not rebuild the edge
        classes; in any mode, it does not rebuild the vertex classes.
        """

        if isinstance(face_or_arrow, Arrow):
            assert tet is None
            arrow = face_or_arrow
            a = arrow.copy()
        else:
            arrow = None
            a = Arrow(PickAnEdge[face_or_arrow], face_or_arrow, tet)

        a = a.copy()
        b = a.glued()
        if not unsafe_mode:
            possible, reason = self._face_permits_two_to_three(a, b)
            if not possible:
                if must_succeed:
                    raise ValueError(reason)
                return False
        a_orig = a.copy()
        new = self.new_arrows(3)
        for i in range(3):
            new[i].glue(new[(i + 1) % 3])
        a.reverse()
        for c in new:
            c.opposite().glue(a.glued())
            c.reverse().glue(b.glued())
            a.rotate(-1)
            b.rotate(1)
        for c in new:
            c.reverse()
            c.opposite()
        self._two_to_three_move_hook(a_orig, new)
        self.delete_tet(a.Tetrahedron)
        self.delete_tet(b.Tetrahedron)
        if not unsafe_mode:
            self.build_edge_classes()
        if VERBOSE:
            print('2->3')
            print(self.EdgeValences)
        if return_arrow:
            return new[1].north_head().get_arrow()
        else:
            return True

    def _edge_permits_three_to_two(self, edge):
        if not edge.IntOrBdry == 'int':
            return False, 'Cannot do move on exterior edge'
        if edge.valence() != 3:
            return False, 'Edge has valence %d not 3' % edge.valence()
        if not edge.distinct():
            return False, 'Tets around edge are not distinct'
        return True, None

    def _three_to_two_move_hook(self, old_arrow, new_arrows):
        pass

    def three_to_two(self, edge_or_arrow, return_arrow=False,
                     must_succeed=False, unsafe_mode=False):
        """
        Replaces the star of an edge of valence 3 by two tetrahedra.

        Options and return value are the same as ``two_to_three``.
        """
        edge, a = edge_and_arrow(edge_or_arrow)
        if not unsafe_mode:
            possible, reason = self._edge_permits_three_to_two(edge)
            if not possible:
                if must_succeed:
                    raise ValueError(reason)
                return False

        a_orig = a.copy()
        b = self.new_arrow()
        c = self.new_arrow()
        b.glue(c)
        b_orig = b.copy()
        b.reverse()
        b_to_return = b.copy()
        for i in range(3):
            b.glue(a.opposite().glued())
            c.glue(a.reverse().glued())
            b.rotate(-1)
            c.rotate(1)
            a.reverse().opposite().next()

        self._three_to_two_move_hook(a_orig, (b_orig, b, c))
        if unsafe_mode:
            tet0 = a_orig.Tetrahedron
            tet1 = a_orig.next().Tetrahedron
            tet2 = a_orig.next().Tetrahedron
            self.delete_tet(tet0)
            self.delete_tet(tet1)
            self.delete_tet(tet2)
        else:
            for corner in edge.Corners:
                self.delete_tet(corner.Tetrahedron)
        if not unsafe_mode:
            self.build_edge_classes()
        if VERBOSE:
            print('3->2')
            print(self.EdgeValences)
        if return_arrow:
            return b_to_return
        return True

    def _arrow_permits_two_to_zero(self, arrow):
        edge = arrow.axis()
        if not edge.IntOrBdry == 'int':
            return False, 'Cannot do move on exterior edge'
        if edge.valence() != 2:
            return False, 'Edge has valence %d not 2' % edge.valence()
        if not edge.distinct():
            return False, 'Tets around edge are not distinct'
        if arrow.equator() == arrow.glued().equator():
            return False, 'Edges opposite the valence 2 edge are the same'
        # You'd think we should exclude the following, but bizarrely
        # everything is fine when this happens, which is quite
        # frequently in some settings.
        #
        # T0 = edge.Corners[0].Tetrahedron
        # T1 = edge.Corners[1].Tetrahedron
        # if (T0 in T0.Neighbor.values()) or (T1 in T1.Neighbor.values()):
        #    return False, 'One tet is glued to itself'
        return True, None

    def _two_to_zero_hook(self, old_arrow):
        pass

    def two_to_zero(self, edge_or_arrow, must_succeed=False, unsafe_mode=False):
        """
        Flatten the star of an edge of valence 2 to eliminate two
        tetrahedra.

        Options and return value are the same as ``two_to_three``.
        """
        edge, a = edge_and_arrow(edge_or_arrow)
        b = a.glued()

        possible, reason = self._arrow_permits_two_to_zero(a)
        if not possible:
            if must_succeed:
                raise ValueError(reason)
            return False

        self._two_to_zero_hook(a)
        a.opposite().glued().reverse().glue(b.opposite().glued())
        a.reverse().glued().reverse().glue(b.reverse().glued())

        for corner in edge.Corners:
            self.delete_tet(corner.Tetrahedron)
        if not unsafe_mode:
            self.build_edge_classes()
            if VERBOSE:
                print('2->0')
                print(self.EdgeValences)
        return True

    def zero_to_two(self, arrow1, gap):
        """
        Blow up two adjacent faces into a pair of tetrahedra.  The
        faces are specified by passing an arrow specifying the first
        face and an integer n.  The second face is obtained by
        reversing the arrow and applying next() n times.  Thus there
        are n faces between the two that are involved in the blow up.
        Returns ``True`` on success, ``False`` if the move cannot be
        performed.
        """
        arrow2 = arrow1.copy().reverse()
        count = 0
        while count < gap:
            if arrow2.next() is None:
                return False
            count = count + 1
        # Do we *also* need the old test, which was
        # b.Tetrahedron == arrow1.Tetrahedron ?
        if arrow1.face_class() == arrow2.face_class():
            return 0
        a = arrow1.glued()
        b = arrow2.glued()
        c = self.new_arrows(2)
        c[0].glue(c[1])
        c[1].glue(c[0])
        c[0].opposite().glue(a)
        c[0].reverse().glue(b)
        c[1].opposite().glue(arrow1.reverse())
        c[1].reverse().glue(arrow2.reverse())
        self.clear_tet(arrow1.Tetrahedron)
        self.clear_tet(arrow2.Tetrahedron)
        self.build_edge_classes()
        if VERBOSE:
            print('0->2')
            print(self.EdgeValences)
        return True

    def _edge_permits_four_to_four(self, edge):
        if not edge.IntOrBdry == 'int':
            return False, 'Cannot do move on exterior edge'
        if edge.valence() != 4:
            return False, 'Edge has valence %d not 4' % edge.valence()
        if not edge.distinct():
            return False, 'Tets around edge are not distinct'
        return True, None

    def _four_to_four_move_hook(self, old_arrow, new_arrows):
        pass

    def four_to_four(self, edge_or_arrow, must_succeed=False, unsafe_mode=False):
        """
        Replace an edge of valence 4 by another diagonal of the
        octahedron formed by the star of the edge.  There are two
        choices for this diagonal.  If you care which one is used then
        pass an arrow representing the edge of valence four.  The head
        of the arrow will be an endpoint of the new diagonal.  If you
        don't care, just pass an edge.  The choice of diagonal will
        then be made randomly.

        Options and return value are the same as ``two_to_three``.
        """
        edge, a = edge_and_arrow(edge_or_arrow)
        a_orig = a.copy()

        possible, reason = self._edge_permits_four_to_four(edge)
        if not possible:
            if must_succeed:
                raise ValueError(reason)
            return False

        c = self.new_arrows(4)
        c_orig = [x.copy() for x in c]
        for i in range(4):
            c[i].glue(c[(i + 1) % 4])
        b = a.glued().reverse()
        c[0].opposite().glue(a.rotate(1).glued())
        c[1].opposite().glue(b.rotate(-1).glued())
        c[2].opposite().glue(b.rotate(-1).glued())
        c[3].opposite().glue(a.rotate(1).glued())
        a.rotate(1).reverse().next()
        b.rotate(-1).reverse().next()
        c[0].reverse().glue(a.rotate(-1).glued())
        c[1].reverse().glue(b.rotate(1).glued())
        c[2].reverse().glue(b.rotate(1).glued())
        c[3].reverse().glue(a.rotate(-1).glued())

        self._four_to_four_move_hook(a_orig, c_orig)
        for corner in edge.Corners:
            self.delete_tet(corner.Tetrahedron)

        if not unsafe_mode:
            self.build_edge_classes()
            if VERBOSE:
                print('4->4')
                print(self.EdgeValences)

        return True

    def attack_valence_one(self):
        """
        Modify the triangulation near a valence 1 edge, creating a
        valence 2 edge that can likely be eliminated, reducing the
        number of tetrahedra by one.
        """
        if len(self) == 1:
            return False
        for e in self.Edges:
            if e.valence() == 1:
                corner = e.Corners[0]
                tet = corner.Tetrahedron
                sub = corner.Subsimplex
                other_faces = [face for face in TwoSubsimplices
                               if not is_subset(sub, face)]
                assert len(other_faces) == 2
                face = other_faces[0]
                self.two_to_three(face, tet, must_succeed=True)
                return True
        return False

    def eliminate_valence_two(self):
        """
        Perform a single ``two_to_zero`` move on a valence 2 edge, if
        any such is possible.
        """
        did_simplify = False
        progress = True
        while progress:
            progress = False
            for edge in self.Edges:
                if edge.valence() == 2:
                    if self.two_to_zero(edge):
                        progress, did_simplify = True, True
                        break
        return did_simplify

    def eliminate_valence_three(self):
        """
        Perform a single ``three_to_two`` move on a valence 3 edge, if
        any such is possible.
        """
        did_simplify = False
        progress = True
        while progress:
            progress = False
            for edge in self.Edges:
                if edge.valence() == 3:
                    if self.three_to_two(edge):
                        progress, did_simplify = True, True
                        break
        return did_simplify

    def easy_simplify(self):
        """
        Perform moves eliminating edges of valence 1, 2, and 3,
        monotonically reducing the number of tetrahedra until no
        further such moves are possible.  Returns whether or not the
        number of tetrahedra was reduced.

        >>> M = Mcomplex('zLALvwvMwLzzAQPQQkbcbeijmoomvwuvust'
        ...              'wwytxtyxyahkswpmakguadppmrssxbkoxsi')
        >>> M.easy_simplify()
        True
        >>> len(M)
        1
        >>> M.rebuild(); M.isosig()
        'bkaagj'
        """

        init_tet = len(self)
        progress = True
        while progress:
            curr_tet = len(self)
            while self.attack_valence_one():
                pass
            while self.eliminate_valence_two() | self.eliminate_valence_three():
                pass
            progress = len(self) < curr_tet

        return len(self) < init_tet

    def jiggle(self):
        """
        Do a random ``four_to_four`` move if one is possible.
        """
        fours = [edge for edge in self.Edges
                 if edge.valence() == 4 and edge.IntOrBdry == 'int']
        if len(fours) == 0:
            return False
        return self.four_to_four(random.choice(fours))

    JIGGLE_LIMIT = 6

    def simplify(self, jiggle_limit=None):
        """
        Try to simplify the triangulation using only moves that do not
        increase the total number of tetrahedra, using a combination
        of ``jiggle`` and ``easy_simplify``.

        When ``jiggle_limit`` is ``None``, it defaults to
        ``self.JIGGLE_LIMIT`` which is typically 6.
        """

        if jiggle_limit is None:
            jiggle_limit = self.JIGGLE_LIMIT

        init_tet = len(self)
        for j in range(jiggle_limit):
            self.easy_simplify()
            if not self.jiggle():
                break
        self.eliminate_valence_two()
        return len(self) < init_tet

    def blowup(self, n):
        """
        Do ``n`` randomly chosen ``two_to_three`` moves.
        """
        for i in range(n):
            rand_tet = self[ random.randint(0, len(self) - 1) ]
            rand_face = TwoSubsimplices[random.randint(0,3)]
            self.two_to_three(rand_face, rand_tet)
            self.eliminate_valence_two()
        return len(self)

    def blowup2(self, n):
        """
        Create ``n`` edges of valence 2 in random places, removing valence
        3 edges whenever they appear.
        """

        for i in range(n):
            rand_edge = self.Edges[ random.randint(0, len(self.Edges) - 1) ]
            j = random.randint(0, len(rand_edge.Corners) - 1)
            k = random.randint(0, len(rand_edge.Corners) - 1 - j)
            one_subsimplex = rand_edge.Corners[j].Subsimplex
            two_subsimplex = LeftFace[one_subsimplex]
            a = Arrow(one_subsimplex, two_subsimplex,
                      rand_edge.Corners[j].Tetrahedron)
            self.zero_to_two(a, k)
            self.eliminate_valence_three()
        return len(self)

    BLOW_UP_MULTIPLE = 6

    def randomize(self, blow_up_multiple=None):
        """
        Do ``blow_up_multiple`` times the current number of tetrahedra
        random ``two_to_three`` moves, and then ``simplify``.

        If ``blow_up_multiple`` is ``None``, it defaults to
        ``self.BLOW_UP_MULTIPLE`` which is typically 6.

        Unlike the other simplification methods, this one rebuilds the
        vertices.
        """
        if blow_up_multiple is None:
            blow_up_multiple = self.BLOW_UP_MULTIPLE
        self.blowup(blow_up_multiple * len(self))
        self.simplify()
        self.rebuild()
        return len(self)

    def bdry_neighbor(self, arrow):
        """
        Find a boundary face adjoining a given boundary face.
        Given an Arrow representing a boundary face, return the Arrow
        representing the boundary face that shares the Arrow's Edge.
        """
        if arrow.next() is not None:
            raise Insanity("That boundary face is not on the boundary!")
        edge = arrow.Tetrahedron.Class[arrow.Edge]
        if edge.LeftBdryArrow == arrow:
            return edge.RightBdryArrow
        else:
            return edge.LeftBdryArrow

    def add_fan(self, edge, n):
        """
        Adds a fan of ``n`` tetrahedra onto a boundary edge and rebuilds.
        """
        if not edge.IntOrBdry == 'bdry':
            return 0
        a = edge.LeftBdryArrow
        b = edge.RightBdryArrow.reverse()
        if n == 0:
            a.glue(b)
            return 1
        new = self.new_arrows(n)
        a.glue( new[0] )
        for j in range(len(new) - 1):
            new[j].glue( new[j + 1] )
        new[-1].glue( b )
        self.rebuild()
        return 1

    def split_star(self,edge):
        """
        Subdivides the star of an edge e.  If the edge has an embedded
        star then this operation first subdivides the edge, producing
        one new vertex and two new edges.  Next each tetrahedron which
        meets the edge is divided into two tetrahedra along a face
        which is the join of the new vertex to the edge opposite to e.
        The edge e must not be self-adjacent in any 2-simplex for this
        operation to be possible.  However, it is allowed for a
        tetrahedron to have two opposite edges identified to e.  In
        this case the tetrahedron is split into four tetrahedra,
        forming the join of two segments of length 2.  In order to
        deal with this situation we work our way around the edge
        making the identifications as we go.  The first time that we
        encounter a corner of a certain tetrahedron it gets split into
        two.  Those two are glued into place and may be encountered
        later in the process, at which time each of them get split in
        two.

        Returns an arrow associated to the "top half" of the original edge
        and the "first" tetrahedron adjacent to that edge, or 0 if the edge
        is self-adjacent.
        """

        if edge.selfadjacent():
            return 0
        # Collect the garbage as we go -- some of the new tets may
        # turn into garbage later on.
        garbage = []
        # Remember where we started.
        first_arrow = edge.get_arrow().next()
        first_bottom,first_top = self.new_arrows(2)
        a = first_arrow.copy()
        bottom = first_bottom.copy()
        top = first_top.copy()
        # Work around the edge.
        while 1:
            garbage.append(a.Tetrahedron)
            # Glue the two new tetrahedra together.
            bottom.glue(top)
            # Attach the top face of our new pair.
            a.opposite()
            above = a.glued()
            if above.is_null():
                # This may mean that our first tetrahedron had opposite edges attached
                # to our edge.  We are splitting it for the second time, which will
                # cause first_top and first_bottom to get dumped onto the garbage heap.
                # We have to create a new first_top or first_bottom here, or we won't
                # be able to close everything up at the end.  (On the other hand, we
                # may just have found a boundary face.)
                check = a.copy().opposite().reverse()
                new_first = top.copy().opposite().reverse()
                if check == first_top:
                    first_top = new_first
                elif check == first_bottom:
                    first_bottom = new_first
            else:
                top.glue(above)
            # Attach the bottom face of the new pair.
            bottom.reverse()
            a.reverse()
            below = a.glued()
            if below.is_null():
                # See comment above.
                check = a.copy().opposite().reverse()
                new_first = bottom.copy().opposite().reverse()
                if check == first_bottom:
                    first_bottom = new_first
                elif check == first_top:
                    first_top = new_first
            else:
                bottom.glue(below)
            bottom.reverse()
            # Now move on around the edge.
            a.reverse()
            a.opposite()
            a.next()
            if a == first_arrow:
                break
            next_bottom, next_top = self.new_arrows(2)
            top.opposite()
            bottom.opposite()
            top.glue(next_top)
            bottom.glue(next_bottom)
            top = next_top.opposite()
            bottom = next_bottom.opposite()
        # OK. We are back to the beginning.  Close it up.
        top.opposite()
        bottom.opposite()
        top.glue(first_top.opposite())
        bottom.glue(first_bottom.opposite())
        # Clean up the garbage.
        for tet in garbage:
            self.delete_tet(tet)
        self.rebuild()
        return first_top

    def smash_star(self, edge):
        """
        If an edge joins distinct vertices and has an embedded open
        star then the following method will smash each 3-simplex in
        the star down to a 2-simplex, and smash the edge to a vertex,
        reducing the number of vertices by 1.  Returns ``True`` on
        success, ``False`` on failure.
        """
        if not edge.distinct() or edge.Vertices[0] == edge.Vertices[1]:
            return False
        start = edge.get_arrow()
        a = start.copy()
        garbage = []
        while 1:
            garbage.append(a.Tetrahedron)
            top = a.opposite().glued()
            bottom = a.reverse().glued().reverse()
            bottom.glue(top)
            a.reverse().opposite().next()
            if a == start:
                break
        for tet in garbage:
            self.delete_tet(tet)
        self.rebuild()
        return True

    def smash_all_edges(self):
        """
        Collapse edges to reduce the number of vertices as much as
        possible. Returns whether the number of vertices has been
        reduced to one.
        """
        success = True
        while len(self.Vertices) > 1 and success:
            success = False
            edges = sorted(self.Edges, key=lambda E:E.valence(), reverse=True)
            edges = self.Edges
            for edge in edges:
                if self.smash_star(edge):
                    success = True
                    break

        return len(self.Vertices) == 1

    def replace_star(self, arrow, top_arrows, bottom_arrows):
        """
        This method takes an arrow and replaces its star with
        other_complex attaching that complex via top_arrows and
        bottom_arrows where: Let a be the arrow defining the same
        directed edge as arrow which is the ith such arrow counting
        around the star.  Then a.glued() is glued to top_arrow[i] and
        a.reverse().glued() is glued to bottom_arrow[i].

        NOTE:  If it fails, you need to delete any tets that you were
        trying to add.
        """

        edge = arrow.Tetrahedron.Class[arrow.Edge]
        a = arrow.copy().opposite()

        # check to make sure that the replacement will work

        if not edge.IntOrBdry == 'int':
            return None
        if not edge.distinct():
            return None
        valence = edge.valence()
        if len(top_arrows) != valence or len(bottom_arrows) != valence:
            return None

        # Attach other_complex to manifold replace star of arrow
        #
        # It's important that we do things incrementally as follows in
        # case two outside faces of the star are glued together.

        for i in range(valence):
            top_arrows[i].glue(a.glued())
            a.reverse()
            bottom_arrows[i].glue(a.glued())
            a.reverse()

            # Now advance a to represent the same edge but in the next
            # tetrahedra in the star.
            a.opposite()
            a.next()
            a.opposite()

        # Delete old star

        for corner in edge.Corners:
            self.delete_tet(corner.Tetrahedron)

        # Rebuild mcomplex
        self.build_edge_classes()
        self.orient()

        return True

    def suspension_of_polygon(self, num_sides_of_polygon):
        """
        This method adds the suspension of a triangulation of a
        polygon to self.Tetrahedra and returns::

          (top_arrows, bottom_arrows)

        Currently the choice of triangulation of the polygon is one
        that is the cone over an edge.  Probably this should be
        generalized.  top_arrows and bottom arrows are for gluing in
        this complex via the method ``replace_star``.
        """
        top_tets = self.new_tets(num_sides_of_polygon - 2)
        bottom_tets = self.new_tets(num_sides_of_polygon - 2)
        n = len(top_tets)

        # glue top and bottom together
        for i in range(n):
            top_tets[i].attach( F3, bottom_tets[i], (0, 2, 1, 3) )

        # glue each tet to its neighbor

        for i in range(n-1):
            top_tets[i].attach(F0, top_tets[i+1], (1, 0, 2, 3) )
            bottom_tets[i].attach(F0, bottom_tets[i+1], (2, 1 , 0, 3) )

        # make arrows

        top_arrows = [ Arrow( comp(E13), F1, top_tets[0]) ]
        bottom_arrows = [ Arrow( comp(E23), F2, bottom_tets[0]) ]
        for i in range(n):
            top_arrows.append(Arrow( comp(E23), F2, top_tets[i]))
            bottom_arrows.append(Arrow( comp(E13), F1, bottom_tets[i]))

        top_arrows.append(Arrow(comp(E03), F0, top_tets[i]))
        bottom_arrows.append(Arrow(comp(E03), F0, bottom_tets[i]))

        return (top_arrows, bottom_arrows)

    def save(self, filename, format="snappy"):
        """
        Nontypical example showing saving to a string buffer:

        >>> import io
        >>> buffer = io.StringIO()
        >>> T = Mcomplex('v3551')
        >>> T.save(buffer, 'snappy')
        >>> T.save(buffer, 'geo')
        >>> T.save(buffer, 'spine')
        >>> len(buffer.getvalue())
        1936
        """
        if not hasattr(filename, 'write'):
            file = open(filename, 'w')
            close = True
        else:
            file = filename
            close = False
        if format == "snappy":
            files.write_SnapPea_file(self, file)
        elif format == "geo":
            files.write_geo_file(self, file)
        elif format == "spine":
            files.write_spine_file(self, file)
        if close:
            file.close()

    def _snappea_file_contents(self):
        data = io.StringIO()
        data.name = 'from_t3m'
        files.write_SnapPea_file(self, data)
        return data.getvalue()

    def snappy_triangulation(self, remove_finite_vertices=True):
        """
        >>> Mcomplex('4_1').snappy_manifold().homology()
        Z

        WARNING: Code implicitly assumes all vertex links are orientable.
        """
        # We don't assume that the indices of the Tetraheda are equal
        # to range(len(self))
        tet_to_index = {T:i for i, T in enumerate(self.Tetrahedra)}

        # Initially set all to -1, which corresponds to a finite vertex
        to_cusp_index = {vertex:-1 for vertex in self.Vertices}
        torus_cusps = 0
        for vertex in self.Vertices:
            g = vertex.link_genus()
            if g > 1:
                raise ValueError('Link of vertex has genus more than 1.')
            if g == 1:
                to_cusp_index[vertex] = torus_cusps
                torus_cusps += 1

        tet_data, cusp_indices, peripheral_curves = [], [], []

        for tet in self.Tetrahedra:
            neighbors, perms = [], []
            for face in TwoSubsimplices:
                if tet.Neighbor[face] is None:
                    raise ValueError('SnapPy triangulations cannot have boundary')

                neighbor = tet_to_index[tet.Neighbor[face]]
                perm = tet.Gluing[face].tuple()
                neighbors.append(neighbor)
                perms.append(perm)
            tet_data.append((neighbors, perms))

            cusp_indices.append([to_cusp_index[tet.Class[vert]]
                                 for vert in ZeroSubsimplices])
            if hasattr(tet, 'PeripheralCurves'):
                for curve in tet.PeripheralCurves:
                    for sheet in curve:
                        one_curve_data = []
                        for v in ZeroSubsimplices:
                            for f in TwoSubsimplices:
                                one_curve_data.append(sheet[v][f])
                        peripheral_curves.append(one_curve_data)
            else:
                for i in range(4):
                    peripheral_curves.append(16*[0])

        M = snappy.Triangulation('empty')
        M._from_tetrahedra_gluing_data(tetrahedra_data=tet_data,
                                       num_or_cusps=torus_cusps,
                                       num_nonor_cusps=0,
                                       cusp_indices=cusp_indices,
                                       peripheral_curves=peripheral_curves,
                                       remove_finite_vertices=remove_finite_vertices)
        return M

    def snappy_manifold(self):
        return self.snappy_triangulation().with_hyperbolic_structure()

    def isosig(self):
        contents = self._snappea_file_contents()
        T = snappy.Triangulation(contents, remove_finite_vertices=False)
        return T.triangulation_isosig(decorated=False)

    def regina_triangulation(self):
        """
        >>> M = Mcomplex('K14n1234')
        >>> try:
        ...     T = M.regina_triangulation()
        ...     assert M.isosig() == T.isoSig()
        ... except ImportError:
        ...     pass
        """
        try:
            import regina
        except ImportError:
            raise ImportError('Regina module not available')

        T = regina.Triangulation3()
        regina_tets = {tet:T.newTetrahedron() for tet in self}
        self.rebuild()
        for face in self.Faces:
            if face.IntOrBdry == 'int':
                corner = face.Corners[0]
                tet0 = corner.Tetrahedron
                face0 = corner.Subsimplex
                tet1 = tet0.Neighbor[face0]
                perm = tet0.Gluing[face0]

                r_tet0 = regina_tets[tet0]
                r_tet1 = regina_tets[tet1]
                r_face = FaceIndex[face0]
                r_perm = regina.Perm4(*perm.tuple())
                r_tet0.join(r_face, r_tet1, r_perm)

        return T

    def boundary_maps(self):
        """
        The boundary maps in the homology chain complex of the
        underlying cell-complex of a Mcomplex.

        >>> M = Mcomplex('o9_12345')
        >>> len(M.boundary_maps()) == 3
        True
        """
        return homology.boundary_maps(self)

    def isomorphisms_to(self, other, orientation_preserving=False, at_most_one=False):
        """
        Return the list of isomorphisms between the MComplexes M and N.
        If `at_most_one` is `True`, only returns the first one found (but
        still as a list).

        >>> tri_data = [([0,1,0,1], [(2,1,0,3), (0,3,2,1), (2,1,0,3), (0,1,3,2)]),
        ...             ([1,1,0,0], [(1,0,2,3), (1,0,2,3), (0,1,3,2), (0,3,2,1)])]
        >>> M = Mcomplex(tri_data)
        >>> N = Mcomplex(M.isosig())
        >>> isos = M.isomorphisms_to(N); len(isos)
        4
        >>> isos[0]
        {0: [tet0, (0, 2, 1, 3)], 1: [tet1, (0, 2, 1, 3)]}
        >>> len(M.isomorphisms_to(N, orientation_preserving=True))
        2
        >>> M.two_to_three(Arrow(E01, F3, M[0])); M.rebuild()
        True
        >>> len(M), len(N)
        (3, 2)
        >>> M.isomorphisms_to(N)
        []
        >>> F = Mcomplex('m004')
        >>> N.isomorphisms_to(F)
        []
        >>> N = Mcomplex(M.isosig())
        >>> M.isomorphisms_to(N, at_most_one=True)[0]
        {0: [tet1, (0, 2, 3, 1)], 1: [tet2, (0, 2, 3, 1)], 2: [tet0, (0, 3, 1, 2)]}
        >>> M = Mcomplex(tri_data)
        >>> M.two_to_three(Arrow(E01, F3, M[0])); M.two_to_three(Arrow(E01, F3, M[1]))
        True
        True
        >>> M.rebuild()
        >>> len(M) == 4
        True
        >>> N = Mcomplex(M.isosig())
        >>> M.isomorphisms_to(N, at_most_one=True)[0]  # doctest: +NORMALIZE_WHITESPACE
        {0: [tet0, (1, 3, 0, 2)], 1: [tet1, (3, 0, 1, 2)],
         2: [tet3, (2, 0, 3, 1)], 3: [tet2, (3, 1, 2, 0)]}
        """
        M, N = self, other
        if not isinstance(N, Mcomplex):
            raise ValueError('The other triangulation must be an Mcomplex')

        if len(M) != len(N):
            return []
        t_M0 = M[0]

        if orientation_preserving:
            if not (M.is_oriented() and N.is_oriented()):
                raise ValueError('Asked for orientation preserving isomorphisms '
                                 'of unoriented triangulations')
            permutations = list(Perm4.A4())  # even perms only
        else:
            permutations = list(Perm4.S4())

        isomorphisms = []
        # We will try and build an isomorphism from M to N that sends t_M0 to t_N0
        for t_N0 in N:
            # for each way t_M can be identified with t_N
            for perm in permutations:
                # initially the map is not defined
                iso = {k:None for k in range(len(M))}
                # set up first map t_M -> t_N
                # temporary way of encoding gluing.
                iso[0] = [t_N0, perm]
                tet_queue = [t_M0]
                while tet_queue != []:
                    t_M = tet_queue.pop()
                    t_N = iso[t_M.Index][0]
                    perm = iso[t_M.Index][1]

                    # Now, for each face F of t_0, package tet that meets
                    # t_0 along F in list neighbors
                    neighbors_M = [t_M.Neighbor[face] for face in TwoSubsimplices]
                    # Need info in N too.
                    neighbors_N = [t_N.Neighbor[perm.image(face)] for face in TwoSubsimplices]

                    # record gluings for each face
                    gluings_M = [t_M.Gluing[face] for face in TwoSubsimplices]
                    gluings_N = [t_N.Gluing[perm.image(face)] for face in TwoSubsimplices]
                    # check compatibility
                    maps = [gluings_N[k]*perm*inv(gluings_M[k]) for k in [0,1,2,3]]

                    # now we try and update iso and hope there are no
                    # incompatibilities
                    for i in range(len(neighbors_M)):
                        t = neighbors_M[i]
                        s = neighbors_N[i]
                        map = maps[i]

                        if iso[t.Index] is not None:
                            if iso[t.Index][0] != s or iso[t.Index][1].tuple() != map.tuple():
                                # not an iso!
                                iso = {k:None for k in range(len(M))} # reset iso
                                tet_queue = [] # clear queue
                                break
                        else:
                            # iso[t] hasn't been set, so we set it and
                            # move forward with it on the queue
                            iso[t.Index] = [s,map]
                            tet_queue = tet_queue + [t]

                # did we succeed, or do we need to reset everything
                if None not in list(iso.values()):  # we succeed!
                    isomorphisms.append(iso.copy())
                    if at_most_one:
                        return isomorphisms
                # otherwise, we failed and the loop goes on

        return isomorphisms


def tets_from_data(fake_tets):
    """
    Takes a list where the ith element represents the gluing data
    for the ith tetraherda::

      ( [Neighbors], [Glueings] )

    and creates the corresponding glued Tetraherda.
    """
    fake_tets = fake_tets
    num_tets = len(fake_tets)
    tets = [Tetrahedron() for i in range(num_tets)]
    for i in range(num_tets):
        neighbors, perms = fake_tets[i]
        for k in range(4):
            tets[i].attach(TwoSubsimplices[k], tets[neighbors[k]], perms[k])
    return tets


def read_geo_file(filename):
    return Mcomplex(tets_from_data(files.read_geo_file(filename)))


def read_SnapPea_file(filename):
    return Mcomplex(tets_from_data(files.read_SnapPea_file(filename)))
