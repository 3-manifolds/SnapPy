#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from operator import inv


def _make_opp_dict():
    def swap(t):
        return (t[1], t[0])
    dic = {(0, 1): (2, 3), (2, 0): (1, 3), (1, 2): (0, 3)}
    for k in list(dic):
        dic[dic[k]] = k
        dic[swap(dic[k])] = swap(k)
        dic[swap(k)] = swap(dic[k])
    return dic


# opposite[(i,j)] = (k,l), where i,j,k,l are distinct and the permutation
# 0->i,1->j,2->k,3->l is even.

opposite = _make_opp_dict()


class Perm4Basic:
    """
    Class Perm4Basic: A permutation of {0,1,2,3}.

    A permutation can be initialized with a length 4 dictionary or a
    4-tuple or a length 2 dictionary; in the latter case the sign of
    the permutation can be specified by setting "sign=0" (even) or
    "sign=1" (odd).  The default sign is odd, since odd permutations
    describe orientation-preserving gluings.

    This is the original version of Perm4, which was tablized for
    speed reasons.
    """

    def __init__(self, init, sign=1):
        self.dict = {}
        if isinstance(init, Perm4Basic) or len(init) == 4:
            for i in range(4):
                self.dict[i] = init[i]
        else:
            self.dict = init
            v = list(init.items())
            x = opposite[(v[0][0],v[1][0])]
            y = opposite[(v[0][1],v[1][1])]
            self.dict[x[0]] = y[sign]
            self.dict[x[1]] = y[1-sign]

    def image(self, bitmap):
        """
        A subset of {0,1,2,3} can be represented by a bitmap.  This
        computes the bitmap of the image subset.

        >>> Perm4Basic([2, 3, 1, 0]).image(10)
        9
        """
        image = 0
        for i in range(4):
            if bitmap & (1 << i):
                image = image | (1 << self.dict[i])
        return image

    def __repr__(self):
        return str(self.tuple())

    def __call__(self, a_tuple):
        """
        P((i, ... ,j)) returns the image tuple (P(i), ... , P(j))
        """
        image = []
        for i in a_tuple:
            image.append(self.dict[i])
        return tuple(image)

    def __getitem__(self, index):
        """
        P[i] returns the image of i, P(i)

        >>> Perm4Basic([2, 3, 1, 0])[3]
        0
        """
        return self.dict[index]

    def __mul__(self, other):
        """
        P*Q is the composition P*Q(i) = P(Q(i))

        >>> P = Perm4Basic((2, 3, 1, 0))
        >>> Q = Perm4Basic((1, 0, 2, 3))
        >>> P * Q
        (3, 2, 1, 0)
        >>> Q * P
        (2, 3, 0, 1)
        """
        composition = {}
        for i in range(4):
            composition[i] = self.dict[other.dict[i]]
        return Perm4Basic(composition)

    def __invert__(self):
        """
        inv(P) is the inverse permutation

        >>> inv(Perm4Basic([2, 1, 3, 0]))
        (3, 1, 0, 2)
        """
        inverse = {}
        for i in range(4):
            inverse[self.dict[i]] = i
        return Perm4Basic(inverse)

    def sign(self):
        """
        sign(P) is the sign: 0 for even, 1 for odd

        >>> Perm4Basic([0, 1, 3, 2]).sign()
        1
        >>> Perm4Basic([1, 0, 3, 2]).sign()
        0
        """
        sign = 0
        for (i,j) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
            sign = sign ^ (self.dict[i] < self.dict[j])
        return sign

    def tuple(self):
        """
        P.tuple() returns a tuple representing the permutation P

        >>> Perm4Basic([1, 2, 0, 3]).tuple()
        (1, 2, 0, 3)
        """
        return tuple(self.dict[i] for i in range(4))


S4_tuples = [(0,1,2,3), (0,1,3,2), (0,2,1,3), (0,2,3,1), (0,3,1,2), (0,3,2,1),
             (1,0,2,3), (1,0,3,2), (1,2,0,3), (1,2,3,0), (1,3,0,2), (1,3,2,0),
             (2,0,1,3), (2,0,3,1), (2,1,0,3), (2,1,3,0), (2,3,0,1), (2,3,1,0),
             (3,0,1,2), (3,0,2,1), (3,1,0,2), (3,1,2,0), (3,2,0,1), (3,2,1,0)]

A4_tuples = [(0,1,2,3), (0,2,3,1), (0,3,1,2), (1,0,3,2), (1,2,0,3), (1,3,2,0),
             (2,0,1,3), (2,1,3,0), (2,3,0,1), (3,0,2,1), (3,1,0,2), (3,2,1,0)]

KleinFour_tuples = [(0,1,2,3),  # Id
                    (1,0,3,2),  # (01)(23)
                    (2,3,0,1),  # (02)(13)
                    (3,2,1,0) ] # (03)(12)

perm_tuple_to_index = {t:i for i, t in enumerate(S4_tuples)}
perm_basic_by_index = [Perm4Basic(t) for t in S4_tuples]


def perm_basic_to_index(perm):
    return perm_tuple_to_index[perm.tuple()]


perm_signs_by_index = {i: perm.sign()
                       for i, perm in enumerate(perm_basic_by_index)}

bitmap_images = {(i, bitmap): perm.image(bitmap)
                 for bitmap in range(16)
                 for i, perm in enumerate(perm_basic_by_index)}

index_of_inverse_by_index = {i: perm_basic_to_index(inv(perm))
                             for i, perm in enumerate(perm_basic_by_index)}

index_mult_table_by_index = {(i, j): perm_basic_to_index(P * Q)
                       for i, P in enumerate(perm_basic_by_index)
                       for j, Q in enumerate(perm_basic_by_index)}


class Perm4():
    """
    Class Perm4: A permutation of {0,1,2,3}.

    A permutation can be initialized with a length 4 dictionary or a
    4-tuple or a length 2 dictionary; in the latter case the sign of
    the permutation can be specified by setting "sign=0" (even) or
    "sign=1" (odd).  The default sign is odd, since odd permutations
    describe orientation-preserving gluings.
    """

    def __init__(self, init, sign=1):
        if isinstance(init, int):
            self._index = init
            self._tuple = S4_tuples[init]
        elif isinstance(init, Perm4):
            self._index = init._index
            self._tuple = init._tuple
        elif len(init) == 4:
            self._tuple = tuple(init[i] for i in range(4))
            self._index = perm_tuple_to_index[self._tuple]
        else:
            self._tuple = Perm4Basic(init, sign).tuple()
            self._index = perm_tuple_to_index[self._tuple]

    def image(self, bitmap):
        """
        A subset of {0,1,2,3} can be represented by a bitmap.  This
        computes the bitmap of the image subset.

        >>> Perm4([2, 3, 1, 0]).image(10)
        9
        """
        return bitmap_images[self._index, bitmap]

    def __repr__(self):
        return str(self._tuple)

    def __call__(self, a_tuple):
        """
        P((i, ... ,j)) returns the image tuple (P(i), ... , P(j))

        >>> Perm4([2, 3, 1, 0])(range(3))
        (2, 3, 1)
        """
        image = []
        for i in a_tuple:
            image.append(self._tuple[i])
        return tuple(image)

    def __getitem__(self, index):
        """
        P[i] returns the image of i, P(i)

        >>> Perm4([2, 3, 1, 0])[3]
        0
        """
        return self._tuple[index]

    def __mul__(self, other):
        """
        P*Q is the composition P*Q(i) = P(Q(i))

        >>> P = Perm4((2, 3, 1, 0))
        >>> Q = Perm4((1, 0, 2, 3))
        >>> P * Q
        (3, 2, 1, 0)
        >>> Q * P
        (2, 3, 0, 1)
        """
        return mult_table_by_index[self._index, other._index]

    def __invert__(self):
        """
        inv(P) is the inverse permutation

        >>> inv(Perm4([2, 1, 3, 0]))
        (3, 1, 0, 2)
        """
        return inverse_by_index[self._index]

    def sign(self):
        """
        sign(P) is the sign: 0 for even, 1 for odd

        >>> Perm4([0, 1, 3, 2]).sign()
        1
        >>> Perm4([1, 0, 3, 2]).sign()
        0
        """
        return perm_signs_by_index[self._index]

    def tuple(self):
        """
        P.tuple() returns a tuple representing the permutation P

        >>> Perm4([1, 2, 0, 3]).tuple()
        (1, 2, 0, 3)
        """
        return self._tuple

    @staticmethod
    def S4():
        """"
        All permutations in S4

        >>> len(list(Perm4.S4()))
        24
        """
        for p in S4_tuples:
            yield Perm4(p)

    @staticmethod
    def A4():
        """
        All even permutations in A4

        >>> len(list(Perm4.A4()))
        12
        """
        for p in A4_tuples:
            yield Perm4(p)

    @staticmethod
    def KleinFour():
        """
        Z/2 x Z/2 as a subgroup of A4.

        >>> len(list(Perm4.KleinFour()))
        4
        """
        for p in KleinFour_tuples:
            yield Perm4(p)


inverse_by_index = {k:Perm4(v) for k, v in index_of_inverse_by_index.items()}
mult_table_by_index = {k:Perm4(v) for k, v in index_mult_table_by_index.items()}

__all__ = ["Perm4", "inv"]

if __name__ == '__main__':
    import doctest
    doctest.testmod()
