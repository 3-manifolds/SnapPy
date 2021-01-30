from hashlib import md5
import array
import re

# Codecs we use.
try:
    unichr
    def encode_torsion(divisors):
        return ''.join(unichr(x) for x in divisors).encode('utf8')
except: # Python3
    def encode_torsion(divisors):
        return ''.join(chr(x) for x in divisors).encode('utf8')


def decode_torsion(utf8):
    return [ord(x) for x in utf8.decode('utf8')]


def encode_matrices(matrices):
    """
    Convert a list of 2x2 integer matrices into a sequence of bytes.
    """
    # The tricky thing here is converting signed integers to bytes.
    return bytes(array.array('b', sum(sum(matrices,[]),[])).tostring())
    # NOTE: tostring is deprecated in python3, but for now
    # it does the same thing as tobytes.


def decode_matrices(byteseq):
    """
    Convert a sequence of 4n bytes into a list of n 2x2 integer matrices.
    """
    m = array.array('b')
    m.fromstring(byteseq)
    return [ [ list(m[n:n+2]), list(m[n+2:n+4]) ]
             for n in range(0, len(m), 4) ]

# Some hash functions for manifolds:

def old_basic_hash(mfld, digits=6):
    return '%%%df'%digits%mfld.volume() + " " + repr(mfld.homology())


def basic_hash(mfld, digits=6):
    if mfld.solution_type() != 'contains degenerate tetrahedra':
        volume = '%%%df'%digits%mfld.volume()
    else:
        volume = 'degenerate'
    return volume + " " + repr(mfld.homology())


def cover_type(mfld):
    return re.findall("~reg~|~irr~|~cyc~", mfld.name())[-1][1:-1]


def cover_hash(mfld, degrees):
    return [ repr(sorted(
        [(cover_type(C), C.homology()) for C in mfld.covers(degree)]
        )) for degree in degrees ]


def old_combined_hash(mfld):
    hash = str(" &and& ".join( [old_basic_hash(mfld)] + cover_hash(mfld, (2,3)) ))
    return hash.encode('utf8')


def combined_hash(mfld):
    hash = str(" &and& ".join( [basic_hash(mfld)] +
                               cover_hash(mfld, (2,3)) ))
    return hash.encode('utf8')


# This one is the hash used in the first version of the database.
def old_db_hash(mfld):
    return md5(old_combined_hash(mfld)).hexdigest()


# This one is used now.
def db_hash(mfld):
    return md5(combined_hash(mfld)).hexdigest()
