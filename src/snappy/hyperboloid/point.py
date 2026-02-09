__all__ = [ 'R13Point' ]

class R13Point:
    def __init__(self, point):
        self.point = point

    def transformed(self,
                    m): # O13-matrix
        return R13Point(m * self.point)
