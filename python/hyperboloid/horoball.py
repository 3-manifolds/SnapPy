__all__ = ['R13Horoball']

class R13Horoball:
    """
    Horoball defined by { x : r13_dot(x, l) > -1 } where l is a
    light-like vector.
    """

    def __init__(self,
                 defining_vec): # Light-like vector
        self.defining_vec = defining_vec

    def transformed(self,
                    m): # O13-matrix
        """
        Returns image of the horoball under given O13-matrix m.
        """

        return R13Horoball(m * self.defining_vec)
