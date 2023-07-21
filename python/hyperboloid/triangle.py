from .line import R13Line
from . import (time_r13_normalise,
               space_r13_normalise,
               r13_dot)

from typing import Sequence

__all__ = ['R13IdealTriangle']

class R13IdealTriangle:
    def __init__(self,
                 plane, # one space-like normal vector
                 bounding_planes, # three space-like normal vectors
                 edges : Sequence[R13Line] # Same order as bounding_planes
                 ):
        self.plane = plane
        self.bounding_planes = bounding_planes
        self.edges = edges

def triangle_bounding_plane(v_opp, v0, v1):
    m = time_r13_normalise(
        -(r13_dot(v1, v_opp) * v0 + r13_dot(v0, v_opp) * v1))

    return _make_r13_unit_tangent_vector(m - v_opp, m)

def _make_r13_unit_tangent_vector(direction, point):
    s = r13_dot(direction, point)
    return space_r13_normalise(direction + s * point)

