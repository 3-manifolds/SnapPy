from .line import R13LineWithMatrix

from ...snap.t3mlite import simplex, Tetrahedron
from ...hyperboloid.distances import distance_r13_lines
from ...hyperboloid.line import R13Line

from typing import Optional

class ObjectCloseToCoreCurve(RuntimeError):
    def __init__(self, obj_name, cusp_index, distance):
        self.obj_name = obj_name
        self.cusp_index = cusp_index
        self.distance = distance
        s = self.obj_name if self.obj_name else "Given geometric object"
        super().__init__(
            "%s is very close to the core curve "
            "of cusp %d and might intersect it. Distance: %r." % (
                s, cusp_index, distance))

def check_away_from_core_curve_iter(iterator, epsilon, obj_name = None):
    for tile in iterator:
        check_away_from_core_curve(
            tile.inverse_lifted_geometric_object,
            tile.lifted_tetrahedron.tet,
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
            raise ObjectCloseToCoreCurve(
                obj_name, tet.Class[v].Index, d)
