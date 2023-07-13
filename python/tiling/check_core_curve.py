from . import exceptions
from .line import R13LineWithMatrix
from ..snap.t3mlite import simplex, Tetrahedron
from ..hyperboloid.distances import distance_r13_lines
from ..hyperboloid.line import R13Line

from typing import Optional

def check_away_from_core_curve_iter(iterator, epsilon, obj_name = None):
    for tile in iterator:
        check_away_from_core_curve(
            tile.lifted_geometric_object.r13_line,
            tile.tet,
            simplex.T,
            epsilon,
            obj_name)

        yield tile

def check_away_from_core_curve(line : R13Line,
                               tet : Tetrahedron,
                               subsimplex : int,
                               epsilon,
                               obj_name = None):
    """
    If the geodesic is intersecting a core curve, the tracing would
    fail in that it would never reach the intersection point and thus
    either hit the iteration limit or breaking down because of
    rounding-errors.

    This function is catching this case to give a meaningful exception
    faster. It does so by computing the distance between the lift of
    the geodesic we are tracing and the lifts of the core curve
    corresponding to the vertices of the tetrahedra adjacent to the
    given face.
    """

    for v in simplex.ZeroSubsimplices:
        if not simplex.is_subset(v, subsimplex):
            continue
        core_curve : Optional[R13LineWithMatrix] = tet.core_curves.get(v, None)
        if core_curve is None:
            continue
        d = distance_r13_lines(
            core_curve.r13_line,
            line)
        if not d > epsilon:
            raise exceptions.ObjectCloseToCoreCurve(
                obj_name, tet.Class[v].Index)
