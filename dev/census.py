from __future__ import print_function

"""
Finding a hyperbolic manifold in a list or census, and for generating
such censuses.

If one wants to find a single manifold in a list of n manifolds, then
there's no way to avoid O(n) behavior. However, when finding many
manifolds in the same fixed list, one can reduce the marginal cost of
a finding a manifold to essentially O(1) by computing a good hash for
each of the manifolds in the list.

The SnapPea isometry checker as several issues that we need to work
around.  When "is_isometric" reports true, one is guaranteed that a
combinatorial homeomorphism of the manifolds exists.  However, it
*can* report false negatives (by miscomputing the canonical
triangulation).  Also, the isometry checker can simply fail, throwing
a RuntimeError.  This happens most often for non-isometric closed
manifolds.  Thus one needs to keep the number of actual isometry
comparisons to a minimum. Our approach for both issues is to use
non-geometric hashes based on covers.
"""

import snappy, re, os, sys


#-----------------------------------------
#
# Misc helper functions
#
#-----------------------------------------

def unique(L):
    return sorted(list(set(L)))

def repetition_count(L):
    return [(x, L.count(x)) for x in unique(L)]

#-----------------------------------------
#
# Hash functions:  manifolds -> strings 
#
#-----------------------------------------

def to_str_at_prec(x, d):
    return ('%.' + repr(d) + 'f') % x

def basic_hash(manifold, digits=6):
    return to_str_at_prec(manifold.volume(), digits) + " " + repr(manifold.homology())

def chern_hash(manifold, digits=6):
    try:
        return  to_str_at_prec(abs(manifold.chern_simons()), digits)
    except ValueError:  #non-hyperbolic
        return 'na'

def basic_plus(manifold, digits=6):
    return basic_hash(manifold, digits) + " " + chern_hash(manifold, digits)

def cover_type(manifold):
    return re.findall("~reg~|~irr~|~cyc~", manifold.name())[-1][1:-1]
                     
def cover_hash(manifold, degree):
    return repr(sorted([(cover_type(C), C.homology()) for C in manifold.covers(degree)]))

def cover_hash_fn(degree):
    return lambda M : cover_hash(M, degree)

def cover_hash_fns(degrees):
    return [cover_hash_fn(d) for d in degrees]

def link_hash(manifold):
    return basic_hash(manifold) + " " + cover_hash(manifold, 2)

class HashFunctions():
    """
    A list of at least one hash function. The first function is the "primary"
    one and is what's evaluated when the HashFunctions is called.
    """
    def __init__(self, hash_fns):
        self.hash_fns = hash_fns

    def __call__(self, manifold):
        return self.hash_fns[0](manifold)

    def all_hashes(self, manifold):
        return [h(manifold) for h in self.hash_fns]

    def combined_hash(self, manifold):
        return " &and& ".join( self.all_hashes(manifold) )

    def can_distinguish(self, M, N):
        for h in self.hash_fns:
            if h(M) != h(N):
                return True

        return False

basic_hashes = HashFunctions([basic_hash])
standard_hashes = HashFunctions([basic_hash, cover_hash_fn(2), cover_hash_fn(3)])
big_cover_hashes = HashFunctions([basic_hash, chern_hash] + cover_hash_fns(range(2, 6)) )

#-----------------------------------------
#
# Simple O(n) code for finding a manifold in a list
#
#-----------------------------------------

def appears_hyperbolic(M):
    acceptable = ['all tetrahedra positively oriented',
                  'contains negatively oriented tetrahedra']
    return M.solution_type() in acceptable and M.volume() > 0.8

def are_isometric(M, N, hash_fns=None, certify=False):
    """
    Return values:

      True :  manifolds are provably isometric
      False : manifolds are (likely/definitely) not isometric
      None : test inclusive/failed

    When certify_distinct is False, a return value of False
    only means that the SnapPea kernel thinks they are
    not isometric; when certify is True, it this
    function only reports False when one of the hashes
    differs.

    A RuntimeError from the kernel is included in the
    "test inclusive" case.

    Manifolds are tested for isometry if the first hash
    function agrees; further hash functions are used
    when kernel fails or if reports False and
    certify is True.  
    """
    if not (appears_hyperbolic(M) and appears_hyperbolic(N)):
        raise ValueError("One of the manifolds is not hyperbolic")

    if hash_fns == None:
        hash_fns = basic_hashes

    if hash_fns(M) != hash_fns(N):
        return False

    try:
        ans = N.is_isometric_to(M)
        if ans:
            return True
    except RuntimeError:
        ans = None

    if certify and ans == False:
        ans = None

    if ans == None:
        if hash_fns.can_distinguish(M, N):
            return False
        
    return ans

def find_hyperbolic_manifold_in_list(M, manifolds, hash_fns=None, certify=False, ignore_errors=False):
    if len(manifolds) == 0:
        return None

    isom_error = False
    for N in manifolds:
        if isinstance(N, str):
            N = snappy.Manifold(N)
        ans = are_isometric(N, M, hash_fns, certify)
        if ans:
            return N

        if ans is None:
            isom_error = True
                
    if isom_error and not ignore_errors:
        raise RuntimeError("SnapPea failed to determine whether the manifolds are isometric.")
    
    return None

#-----------------------------------------------
#
# Storing manifolds by hashes for quick comparison
#
#-----------------------------------------------

class ObjectsByHashes:
    """
    For objects "o" store them by their hashes "h".
    There can be multiple objects with the same hash.
    """
    def __init__(self, objects_with_hashes=None):
        self.hashes = {}
        if objects_with_hashes:
            for o, h in objects_with_hashes:
                self.add(o,h)

    def __getitem__(self, h):
        if h in self.hashes:
            return self.hashes[h]
        return tuple()
        
    def add(self, o, h):
        self.hashes[h] = self[h] + (o,)

    def num_hashes(self):
        return len(self.hashes.keys())

    def num_objects(self):
        return sum([len(x) for x in self.hashes.values()])

    def object_group_sizes(self):
        sizes = [len(v) for v in self.hashes.values()]
        return repetition_count(sizes)

    def __iter__(self):
        return self.hashes.iteritems()

    def __repr__(self):
        return "<ObjHashes %d hashes %d items>" % (self.num_hashes(), self.num_objects())


class ManifoldRecognizer():
    """
    Class for finding manifolds in census or other
    collection of manifolds.

    Manifolds can be add dynamically and searching
    the initial list of manifolds is lazy.  
    """
    def __init__(self, init_manifolds=None, hash_fns=None):
        self.obj_hashes = ObjectsByHashes()
        self.init_manifolds = init_manifolds if init_manifolds else []
        self.init_manifold_iter = iter(init_manifolds)
        self.init_manifolds_left = len(init_manifolds)
        if hash_fns == None:
            hash_fns = basic_hashes
        self.hash = hash_fns

    def find(self, manifold):
        h = self.hash(manifold)
        manifolds = self.obj_hashes[h]
        init_all_added = self.init_manifolds_left == 0
        M = find_hyperbolic_manifold_in_list(manifold, manifolds, self.hash,
                                             certify=init_all_added, ignore_errors=not init_all_added) 
        if M == None and not init_all_added:
            self.add_init_manifolds_until_hash(h)
            return self.find(manifold)
        return M

    def add(self, manifold):
        """
        Adds the manifold to the collection and returns its hash.
        """
        h = self.hash(manifold)
        self.obj_hashes.add(manifold, h)
        return h

    def add_init_manifold(self):
        try:
            M = self.init_manifold_iter.next()
            self.init_manifolds_left += -1
            return self.add(M)
        except StopIteration:
            return None

    def add_init_manifolds_until_hash(self, target_hash=None):
        while self.init_manifolds_left > 0:
            if self.add_init_manifold() == target_hash:
                return

    def basic_hash_collisions(self):
        return sum( [list(mflds) for h, mflds in self.obj_hashes if len(mflds) > 1], [])

    def full_hash_collisions(self, hash_fns = None):
        if hash_fns == None:
            hash_fns = self.hash
        h = HashFunctions([hash_fns.combined_hash])
        obj_h = ObjectsByHashes( [ (M, h(M)) for M in self.basic_hash_collisions()] )
        return obj_h

    def test_hashes(self, hash_fns=None):
        print(self)
        coarse = self.obj_hashes.object_group_sizes()
        c, ms = coarse[0]
        already_resolved = ms if c == 1 else 0
        print(coarse)
        obj_h = self.full_hash_collisions(hash_fns)
        print(obj_h.object_group_sizes(), already_resolved + obj_h.num_hashes())
            
    def __repr__(self):
        h, hm, l = self.obj_hashes.num_hashes(), self.obj_hashes.num_objects(), self.init_manifolds_left
        return "<MfldRecog: %d hashes of %d mflds with %d mflds unhashed>" % (h, hm, l)
        


class OrientableCuspedCensusRecognizer(ManifoldRecognizer):
    """
    17,661 manifolds.

    Basic hash is volume which gives 17,037 classes.

    Other hashes add CS and homology of covers of degree 2 and 3 and gives 17,660 classes.
    """
    def __init__(self):
        ManifoldRecognizer.__init__(self, snappy.OrientableCuspedCensus(), standard_hashes)

class OrientableClosedCensusRecognizer(ManifoldRecognizer):
    """
    11,031 manifolds.

    Basic hash is volume and gives 10,885 classes.

    Other hashes add CS and homology of covers of degree 2 and 3 and gives 11,022 classes.
    """
    def __init__(self):
        ManifoldRecognizer.__init__(self, snappy.OrientableClosedCensus(), standard_hashes)


def add_recognizer_to_census(recog, census):
    census._regonizer = recog()
    def find(self, manifold):
        return self._regonizer.find(manifold)
    census.find = find

add_recognizer_to_census(OrientableClosedCensusRecognizer, snappy.OrientableClosedCensus)
add_recognizer_to_census(OrientableCuspedCensusRecognizer, snappy.OrientableCuspedCensus)


def test_census_recognizer(census):
    import random
    mflds = [M for M in census()]
    random.shuffle(mflds)
    CR = census._regonizer
    for M in mflds:
        print( M, CR.find(M), CR )

    CR.test_hashes()

#test_census_recognizer(snappy.OrientableCuspedCensus)
