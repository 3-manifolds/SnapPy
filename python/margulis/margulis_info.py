from ..SnapPy import Info

from typing import Optional

class MargulisInfo(Info):
    pass

class MargulisTubeInfo(MargulisInfo):
    def __repr__(self):
        return "Tube about geodesic %s%s of radius %r" % (
            self.word, _format_core_curve(self.core_curve), self.radius)

class MargulisCuspNeighborhoodInfo(MargulisInfo):
    def __repr__(self):
        return "Cusp Neighborhood for cusp %d of area %r" % (
            self.cusp_index, self.cusp_area)

def _format_core_curve(core_curve : Optional[int]):
    if core_curve is None:
        return ''
    return ' (Core curve of cusp %d)' % core_curve
