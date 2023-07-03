from ..hyperboloid import r13_dot
from ..tiling.line import R13Line, R13LineWithMatrix

def balance_end_points_of_line(line_with_matrix : R13LineWithMatrix,
                               point) -> R13LineWithMatrix:
    return R13LineWithMatrix(
        R13Line(
            [ endpoint / -r13_dot(point, endpoint)
              for endpoint in line_with_matrix.r13_line.points]),
        line_with_matrix.o13_matrix)
