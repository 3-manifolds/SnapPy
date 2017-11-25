""" 
Computing the defining equations of the SL(2, C) character variety
in terms of the trace coordinates.

Code contributed by Jean-Philippe Burelle
"""

import snappy
from snappy.pari import pari
import string
from itertools import combinations, combinations_with_replacement, product

def cycle_sort(l):
    """
    Utility function which takes a list l and returns the minimum
    for the alphabetical order among all cyclic permutations of the list.
    """
    s = l
    for i in range(0,len(l)):
        temp = l[i:] + l[0:i]
        if temp<s:
            s = temp
    return s

def pari_poly_variable(variable_name):
    """
    Ensures that PARI has the requested polynomial variable defined.
    If "variable_name" is already defined in PARI as something else,
    an exception is raised.

    >>> val = 3*pari_poly_variable('silly')**2; val
    3*silly^2
    >>> ten = pari('silly = 10')
    >>> pari_poly_variable('silly')
    Traceback (most recent call last):
    ...
    RuntimeError: In PARI, "silly" is already defined
    """
    p = pari(variable_name)
    success = p.type() == 't_POL' and p.variables() == [p]
    if not success:
        raise RuntimeError('In PARI, "%s" is already defined' % variable_name)
    return p


class TracePolynomialRing(object):
    """
    >>> R = TracePolynomialRing('abc')
    >>> R.var_names
    ['Ta', 'Tb', 'Tc', 'Tab', 'Tac', 'Tbc', 'Tabc']
    >>> R('Ta*Tb')
    Tb*Ta
    """
    def __init__(self, gens):
        self._set_var_names(gens)
        self.vars = [pari_poly_variable(v) for v in self.var_names]

    def _set_var_names(self, gens):
        if (len(set(gens)) != len(gens) or
            not set(gens).issubset(string.ascii_lowercase)):
            raise ValueError('Generators are unsuitable')
        poly_vars = list(gens) + list(combinations(gens, 2))
        poly_vars += list(combinations(gens, 3))
        self.var_names = ['T' + ''.join(v) for v in poly_vars]
        
    def __call__(self, poly):
        return pari(poly)

class Word(object):
    """
    The Word class is used to make objects which represent words in a
    free group. The words are represented by a string of letters, with
    capital letters standing for inverses. ex: 'abAB'. The string is
    contained in the 'letters' attribute of the class.
    """
    
    def __init__(self,letters):
        """Creates a Word from a string, automatically reduces"""
        if isinstance(letters, Word):
            s = letters.letters
        else:
            s = letters
            # As long as the word is not reduced, delete all substrings xX or Xx
            while not self.is_reduced(s):
                for i,j in zip(s,s[1:]):
                    if i != j and (i == j.upper() or i == j.lower()):
                        s = s.replace(i + j,"")
        self.letters = s

    def __repr__(self):
        return self.letters

    def __mul__(self,other):
        return Word(self.letters+other.letters)

    def is_reduced(self,s):
        """Returns true if and only if the string s represents a reduced word"""
        for i,j in zip(s,s[1:]):
            if i != j and (i == j.upper() or i == j.lower()):
                break
        else:
            return True
        return False

    def SL2_trace(self):
        """
        Returns the simplified SL(2) trace of the Word represented by this
        object. The format of the output is a Pari polynomial in the
        variables Tw where w is a word of length 3 or less

        Examples:
        >>> Word("a").SL2_trace()
        Ta
        >>> Word("A").SL2_trace()
        Ta
        >>> Word("aa").SL2_trace()
        Ta^2 - 2
        >>> Word("abAB").SL2_trace()
        Ta^2 - Tab*Tb*Ta + (Tb^2 + (Tab^2 - 2))
        >>> Word("abca").SL2_trace()
        Tabc*Ta - Tbc
        """
        if self.letters == "":
            return pari("2") # The SL(2,C) trace of the identity word is 2.

        # Cyclically permute the letters until they are minimal for the
        # lexicographic ordering.
        s = cycle_sort(self.letters)

        # Reduction of traces with a double letter
        for i,j in zip(s,s[1:]):
            if i == j:
                [w1,c,w2] = s.partition(i+i)
                return (Word(i).SL2_trace()*Word(i+w2+w1).SL2_trace()
                - Word(w2+w1).SL2_trace())

        # Reduction of traces with inverses
        for i in s:
            if i.isupper():
                [w1,c,w2] = s.partition(i)
                return (Word(i.lower()).SL2_trace()*Word(w2+w1).SL2_trace()
                - Word(w1+i.lower()+w2).SL2_trace())

        # Reductions of traces of length larger than 4
        if len(s) >= 4:
            [x,y,z,w] = [s[0],s[1],s[2],s[3:]]
            return (pari("1/2")*(tr(Word(x))*tr(Word(y))*tr(Word(z))*tr(Word(w))
            + tr(Word(x))*tr(Word(y+z+w)) + tr(Word(y))*tr(Word(x+z+w))
            + tr(Word(z))*tr(Word(x+y+w)) + tr(Word(w))*tr(Word(x+y+z))
            - tr(Word(x+z))*tr(Word(y+w)) + tr(Word(x+w))*tr(Word(y+z))
            + tr(Word(x+y))*tr(Word(z+w))
            - tr(Word(x))*tr(Word(y))*tr(Word(z+w))
            - tr(Word(x))*tr(Word(w))*tr(Word(y+z))
            - tr(Word(y))*tr(Word(z))*tr(Word(x+w))
            - tr(Word(z))*tr(Word(w))*tr(Word(x+y))))

        # If the word is length 3 but not lexicographically sorted, we use a
        # trace identity to express it in terms of lexicographically
        # ordered words
        if len(s) == 3 and s != ''.join(sorted(s)):
            [x,y,z] = s
            return (-Word(x+z+y).SL2_trace()
            + Word(x).SL2_trace()*Word(y+z).SL2_trace()
            + Word(y).SL2_trace()*Word(x+z).SL2_trace()
            + Word(z).SL2_trace()*Word(x+y).SL2_trace()
            - Word(x).SL2_trace()*Word(y).SL2_trace()*Word(z).SL2_trace())

        # Output the trace if is one of the generators (length 3 or less,
        # alphabetical order)
        if len(s) <= 3 and s.islower():
            return pari("T"+s)

def tr(w):
    """Shortcut for the SL2_trace method of a word object"""
    return w.SL2_trace()


class Presentation(object):
    """
    Class representing a presentation of a finitely presented group.
    gens is a list of Word objects representing the generators
    rels is a list of Word objects representing the relations.
    """
    def __init__(self,G,R):
        """Creates a Presentation from G:generators and R:relations"""
        self.gens = G
        self.rels = R

    def __repr__(self):
        r = ''
        for w in self.rels:
            r += ('\n' + str(w))
        return "Generators\n" + str(self.gens) + "\n"+"Relations" + r


############################
# Relations in the character variety

def mult_traceless(a1,a2,a3=None):
    """
    Takes 2 or 3 words and returns the trace of their product after
    making them traceless via M-> M-1/2tr(M)*I.
    """
    if a3 is None:
        return tr(a1*a2)-pari("1/2")*tr(a1)*tr(a2)
    else:
        return (pari("5/8")*tr(a1)*tr(a2)*tr(a3) - pari("1/2")*tr(a3)*tr(a1*a2)
                - pari("1/2")*tr(a2)*tr(a1*a3) - pari("1/2")*tr(a1)*tr(a2*a3)
                + tr(a1*a2*a3))

def s3(a1,a2,a3):
    """
    Accessory function to sum (with sign) "mult_traceless" over all
    permutations of three arguments. Used in defining rel1.
    """
    return (mult_traceless(a1,a2,a3) + mult_traceless(a2,a3,a1)
            + mult_traceless(a3,a1,a2) - mult_traceless(a1,a3,a2)
            - mult_traceless(a3,a2,a1) - mult_traceless(a2,a1,a3))

def det(M):
    """Determinant of a 3x3 matrix"""
    return (M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0]
            + M[1][0]*M[2][1]*M[0][2] - M[0][2]*M[1][1]*M[2][0]
            - M[0][1]*M[1][0]*M[2][2] - M[0][0]*M[1][2]*M[2][1])

def rel1(i):
    """Generates type 1 relations for generators (words) i1,i2,i3 j1,j2,j3"""
    [[i1,i2,i3],[j1,j2,j3]] = i
    return (s3(i1,i2,i3)*s3(j1,j2,j3) + 18*det(
        [[mult_traceless(i1,j1),mult_traceless(i1,j2),mult_traceless(i1,j3)],
         [mult_traceless(i2,j1),mult_traceless(i2,j2),mult_traceless(i2,j3)],
         [mult_traceless(i3,j1),mult_traceless(i3,j2),mult_traceless(i3,j3)]]))

def rel2(j):
    """Generates type 2 relations for generators (words) i,p0,p1,p2,p3"""
    [i,[p0,p1,p2,p3]] = j
    return (mult_traceless(i,p0)*s3(p1,p2,p3)
            - mult_traceless(i,p1)*s3(p0,p2,p3)
            + mult_traceless(i,p2)*s3(p0,p1,p3) - mult_traceless(i,p3)*s3(p0,p1,p2))

def rels_from_rel(R, G):
    """
    Returns the relations in the character variety coming from a relation
    in the group presentation. The input is:

    R - a word object, the relation in the group
    G - a list of words, the set of generators of the group
    """
    relations = [tr(R*g)-tr(g) for g in G]
    relations = relations + [tr(R)-tr(Word(""))]
    return [r for r in relations if r != 0]


def character_variety(gens, rels=None):
    """
    Takes a list of generators and relators, either as Words or as
    plain strings, and returns a Presentation object containing
    generators and relations for the character variety of the group
    generated by gens and with relations rels. You can also give a
    SnapPy fundamental group as the sole argument.

    Examples:

    >>> character_variety([Word("a"),Word("b")],[Word("aba")])
    Generators
    [Ta, Tb, Tab]
    Relations
    Tab*Ta^2 + (-Tb - 1)*Ta - Tab
    (-Tb + (Tab^2 - 2))
    Tab*Ta + (-Tb - 2)
    >>> character_variety(["a","b"],["abAB"])
    Generators
    [Ta, Tb, Tab]
    Relations
    Ta^3 - Tab*Tb*Ta^2 + (Tb^2 + (Tab^2 - 4))*Ta
    Ta^2 - Tab*Tb*Ta + (Tb^2 + (Tab^2 - 4))
    >>> character_variety("abc",[])
    Generators
    [Ta, Tb, Tc, Tab, Tac, Tbc, Tabc]
    Relations
    36*Ta^2 + ((36*Tabc*Tc - 36*Tab)*Tb + (-36*Tac*Tc - 36*Tabc*Tbc))*Ta + (36*Tb^2 + (-36*Tbc*Tc - 36*Tabc*Tac)*Tb + (36*Tc^2 - 36*Tabc*Tab*Tc + (36*Tab^2 + 36*Tbc*Tac*Tab + (36*Tac^2 + (36*Tbc^2 + (36*Tabc^2 - 144))))))
    >>> len(character_variety("abcd",[]).rels)
    14

    >>> character_variety("xy",["xy"])
    Generators
    [Tx, Ty, Txy]
    Relations
    (Txy - 1)*Tx - Ty
    -Tx + (Txy - 1)*Ty
    Txy - 2

    >>> H = snappy.Manifold('dLQacccbjkg')        # Hopf link exterior.
    >>> character_variety(H.fundamental_group())  # Answer copied from above.
    Generators
    [Ta, Tb, Tab]
    Relations
    Ta^3 - Tab*Tb*Ta^2 + (Tb^2 + (Tab^2 - 4))*Ta
    Ta^2 - Tab*Tb*Ta + (Tb^2 + (Tab^2 - 4))

    >>> G = snappy.Manifold('L6a5').fundamental_group()
    >>> V = character_variety(G)
    >>> len(V.gens), len(V.rels)
    (7, 9)
    """
    if rels is None:  # SnapPy group
        G = gens
        gens, rels = G.generators(), G.relators()
    gens = [Word(gen) for gen in gens]
    rels = [Word(R) for R in rels]
    ring = TracePolynomialRing([g.letters for g in gens])
    
    #Type 1
    triples = list(combinations(gens,3))
    pairsoftriples = list(combinations_with_replacement(triples,2))

    t1 = [rel1(i) for i in pairsoftriples]

    #Type 2
    fours = list(combinations(gens,4))
    indices = product(gens,fours)

    t2 = [rel2(j) for j in indices]

    #Relations from relations
    r = []
    for R in rels:
        r += rels_from_rel(R,gens)

    return Presentation(ring.vars, t1+t2+r)

def character_variety_ideal(gens, rels=None):
    """
    sage: M = snappy.Manifold('m004')
    sage: I = character_variety_ideal(M.fundamental_group())
    sage: I.dimension()
    1
    sage: len(I.radical().primary_decomposition())
    2
    """
    pres = character_variety(gens, rels)
    from sage.all import PolynomialRing, QQ
    R = PolynomialRing(QQ, [repr(v) for v in pres.gens])
    I = R.ideal([R(p) for p in pres.rels])
    return I

def total_answer_length(I):
    return sum([len(list(p)) for p in I.gens()])

if __name__ == "__main__":
   from snappy.sage_helper import _within_sage, doctest_modules
   import sys
   current_module = sys.modules[__name__]
   doctest_modules([current_module])

   
