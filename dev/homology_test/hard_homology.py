import os, sys, re, gzip, snappy

class TestManifold:
    def __init__(self, volume, name, degree, hom):
        self.volume, self.name, self.degree, self.hom = volume, name, degree, hom

    def cover(self):
        return self.base().cover(self.hom)

    def base(self):
        return snappy.Manifold(self.name)
    
    def __repr__(self):
        return "<%s %d %f>" % (self.name, self.degree, self.volume)

def manifolds():
    return [ TestManifold(*[eval(a) for a in line[:-1].split("\t")])
             for line in gzip.open("hard_homology.data.gz").readlines()]

def main():
    count = 0
    for M in manifolds():
        # Some of the examples don't work, probably because of presentation issues.
        try:
            C = M.cover()
        except:
            continue
        count += 1
        if count > 675:
            break
        hom = C.big_homology()
        #assert C.small_homology() == hom
        #hom = C.small_homology()
        print count, C.num_tetrahedra(), hom
         

if __name__ == "__main__":
    main()
    
