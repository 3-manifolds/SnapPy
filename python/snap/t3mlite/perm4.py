#$Id: perm4.py,v 1.3 2002/09/20 03:52:16 culler Exp $
#   t3m - software for studying triangulated 3-manifolds
#   Copyright (C) 2002 Marc Culler, Nathan Dunfield and others
#
#   This program is distributed under the terms of the 
#   GNU General Public License, version 2 or later, as published by
#   the Free Software Foundation.  See the file GPL.txt for details.

from operator import inv

def _make_opp_dict():
      def swap(t): return (t[1], t[0])
      dict = {(0,1) : (2,3), (2, 0):(1, 3), (1,2):(0,3)}
      for k in list(dict.keys()):
        dict[dict[k]] = k
        dict[swap(dict[k])] = swap(k)
        dict[swap(k)] = swap(dict[k])
      return dict

# Class Perm4: A permutation of {0,1,2,3}.
# A permutation can be initialized with a length 4 dictionary or a
# 4-tuple or a length 2 dictionary; in the latter case the sign of the
# permutation can be specified by setting "sign=0" (even) or "sign=1"
# (odd).  The default sign is odd, since odd permutations describe
# orientation-preserving gluings.

class Perm4:

  def __init__(self, init, sign=1):
    self.dict = {}
    if len(init) == 4:
      for i in range(4):
        self.dict[i] = init[i]
    else:
      self.dict = init
      v = init.items()
      x = self.opposite[(v[0][0],v[1][0])]
      y = self.opposite[(v[0][1],v[1][1])]
      self.dict[x[0]] = y[sign]
      self.dict[x[1]] = y[1-sign]  

  # opposite[(i,j)] = (k,l), where i,j,k,l are distinct and the permutation
  # 0->i,1->j,2->k,3->l is even.
  opposite = _make_opp_dict()

  # A subset of {0,1,2,3} can be represented by a bitmap.  This computes
  # the bitmap of the image subset.
  def image(self, bitmap):
    image = 0
    mask = 1
    for i in range(4):
      if bitmap & (1 << i):
        image = image | (1 << self.dict[i])
    return  image

  def __repr__(self):
     return(str(self.tuple()))

  # P((i, ... ,j)) returns the image tuple (P(i), ... , P(j))
  def __call__(self, a_tuple):
    image = []
    for i in a_tuple:
      image.append(self.dict[i])
    return tuple(image)

  # P[i] returns the image of i, P(i)
  def __getitem__(self, index):
    return self.dict[index]

  # P*Q is the composition P*Q(i) = P(Q(i))
  def __mul__(self, other):
    composition = {}
    for i in range(4):
      composition[i] = self.dict[other.dict[i]]
    return Perm4(composition)

  # inv(P) is the inverse permutation 
  def __invert__(self):
    inverse = {}
    for i in range(4):
      inverse[self.dict[i]] = i
    return Perm4(inverse)

  # sign(P) is the sign: 0 for even, 1 for odd
  def sign(self):
    sign = 0
    for (i,j) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
      sign = sign ^ (self.dict[i] < self.dict[j])
    return sign

  # P.tuple() returns a tuple representing the permutation P
  def tuple(self):
    the_tuple = ()
    for i in range(4):
     the_tuple = the_tuple + (self.dict[i],)
    return the_tuple

  _rawS4 = [ (0,1,2,3),
             (0,1,3,2),
             (0,2,1,3),
             (0,2,3,1),
             (0,3,1,2),
             (0,3,2,1),
             (1,0,2,3),
             (1,0,3,2),
             (1,2,0,3),
             (1,2,3,0),
             (1,3,0,2),
             (1,3,2,0),
             (2,0,1,3),
             (2,0,3,1),
             (2,1,0,3),
             (2,1,3,0),
             (2,3,0,1),
             (2,3,1,0),
             (3,0,1,2),
             (3,0,2,1),
             (3,1,0,2),
             (3,1,2,0),
             (3,2,0,1),
             (3,2,1,0) ]


  @staticmethod
  def S4():
    """" 
    All permutations in S4
    """
    
    for p in Perm4._rawS4:
      yield Perm4(p)
