from ..geometric_structure.geodesic.fixed_points import r13_fixed_line_of_psl2c_matrix
from ..geometric_structure.geodesic.geodesic_start_point_info import sample_line, GeodesicStartPointInfo
from ..geometric_structure.geodesic.line import R13LineWithMatrix
from ..hyperboloid.line import R13Line
from ..hyperboloid import so13_to_psl2c
from ..upper_halfspace import complex_length_of_psl2c_matrix
from ..math_basics import lower # type: ignore
from ..snap.t3mlite import Mcomplex

from typing import List, Optional

class GeodesicInfoBase:
    """
    Basic information about a geodesic, consisting of word and matrix.

    Used intermediately during the computation of the length spectrum.
    Ordered by (lower bound of) real length.

    After de-duplication, it is converted to the user-facing
    LengthSpectrumGeodesicInfo.
    """
    def __init__(self,
                 word : List[int],
                 o13_matrix):
        self.word = word
        self.o13_matrix = o13_matrix

        self.psl2c_matrix = so13_to_psl2c(self.o13_matrix)
        self.length = complex_length_of_psl2c_matrix(self.psl2c_matrix)
        self._key = lower(self.length.real())

    def __lt__(self, other):
        """
        Ordering <.
        """
        return self._key < other._key

class CoreCurveGeodesicInfo(GeodesicInfoBase):
    """
    Information for a known core curve. That is, we go through each
    filled cusp and compute this information before starting to compute
    the length spectrum.

    It contains the index of the corresponding (filled) cusp.
    """
    def __init__(self,
                 word : List[int],
                 o13_matrix,
                 core_curve : int):
        super().__init__(word, o13_matrix)
        self.core_curve = core_curve

class GeodesicKeyInfo(GeodesicInfoBase):
    """
    Information for a geodesic which might potentially a multiple of another
    geodesic or even a core curve.

    Via geodesic_start_point_info, we can determine whether this geodesic
    corresponds to a multiple of a core curve.

    Given two geodesics that are not core curves, we can also use this
    information here to determine whether one is a conjugate of a multiple of
    the other.

    Also see get_geodesic_key_info_dict and get_geodesic_key_info_set.
    """
    def __init__(self,
                 mcomplex : Mcomplex,
                 word : List[int],
                 o13_matrix):
        super().__init__(word, o13_matrix)
        self.mcomplex = mcomplex

        self._r13_line_with_matrix : Optional[R13LineWithMatrix] = None
        self._info : Optional[GeodesicStartPointInfo] = None

    def r13_line_with_matrix(self) -> R13LineWithMatrix:
        """
        The actual line in H^3 with the matrix corresponding to the
        geodesic.
        """
        if self._r13_line_with_matrix is None:
            self._r13_line_with_matrix = (
                r13_fixed_line_of_psl2c_matrix(self.psl2c_matrix))

        return self._r13_line_with_matrix

    def r13_line(self) -> R13Line:
        """
        The actually line in H^3 corresponding to the geodesic.
        """
        return self.r13_line_with_matrix().r13_line

    def geodesic_start_point_info(self) -> GeodesicStartPointInfo:
        """
        Information to start developing about the geodesic.
        """
        if self._info is None:
            start_point = sample_line(self.r13_line())

            self._info = GeodesicStartPointInfo(
                mcomplex=self.mcomplex,
                word=self.word,
                trace=self.psl2c_matrix.trace(),
                unnormalised_start_point = start_point,
                unnormalised_end_point=self.o13_matrix * start_point,
                line=self.r13_line_with_matrix())
            self._info.find_tet_or_core_curve()

        return self._info
