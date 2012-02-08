from snappy import *
import os
import sqlite3
import binascii
from hashlib import md5
import array
from census import *

schemas = [
"""
CREATE TABLE cusped_census (
 id integer primary key,
 name text,
 volume real,
 chernsimons real,
 hash blob,
 triangulation blob)
"""
]

cusped_insert_query = """insert into cusped_census
(name, volume, chernsimons, hash, triangulation)
values ('%s', %s, %s, X'%s', X'%s')"""

USE_COBS = 1 << 7
USE_LONG_STRING = 1 << 6
epsilon = 0.000001

def flatten_matrices(matrices):
    """
    Convert a list of 2x2 integer matrices into a sequence of bytes.
    """
    # The tricky point here is converting signed integers to bytes.
    # tostring is deprecated in python3, but it's the same as tobytes.
    return bytes(array.array('b', sum(sum(matrices,[]),[])).tostring())
    
def create_manifold_db(connection):
    """
    Create the empty tables for our manifold database.
    """
    for schema in schemas:
        connection.execute(schema)
        connection.commit()

def ambiguity_exists(M):
    """
    Does this manifold seem non-hyperbolic, or does it seem to have
    a square or hexagonal cusp?
    """
    if M.solution_type() != 'all tetrahedra positively oriented':
        return True
    for shape in [c.shape for c in M.cusp_info()]: 
        if abs(shape**3+1) < epsilon or abs(shape**2+1) < epsilon:
            return True
    return False

def get_header(mfld, is_link=False, use_string=False):
    """
    Build the header byte for the triangulation string.
    The high order bit indicates that basis change is needed.
    The next bit indicates that the manifold should be built from a
    string containing a triangulation file, rather than a terse
    triangulation string.
    The rest of the byte holds the number of cusps.
    """
    header = mfld.num_cusps()
    if is_link or ambiguity_exists(mfld):
        header |= USE_COBS
    if use_string:
        header |= USE_LONG_STRING
    return header&USE_COBS, bytes(bytearray([header]))
    
def insert_cusped_manifold(connection, mfld, is_link=False, use_string=False):
    name = mfld.name()
    volume = mfld.volume()
    cs = mfld.chern_simons()
    use_cobs, triangulation = get_header(mfld, is_link, use_string)
    if use_cobs:
        cobs = mfld.set_peripheral_curves('combinatorial')
        mfld.set_peripheral_curves(cobs) # undo the basis change
        triangulation += flatten_matrices(cobs)
    if use_string:
        triangulation += mfld._to_string()
    else:
        triangulation += mfld._to_bytes()
    triangulation = binascii.hexlify(triangulation)
    hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
    query = cusped_insert_query%(name, volume, cs, hash, triangulation)
    try:
        connection.execute(query)
    except:
        print query

def make_cusped(connection):
    for M in OrientableCuspedCensus():
        insert_cusped_manifold(connection, M)
    connection.commit()

def make_closed(connection):
    # Doesn't work yet!
    for M in OrientableClosedCensus():
        insert_closed_manifold(connection, M)
    connection.commit()

if __name__ == '__main__':
    if os.path.exists('manfolds.sqlite'):
        print '%s already exists!'%'manifolds.sqlite'
    connection = sqlite3.connect('manifolds.sqlite')
    create_manifold_db(connection)
    make_cusped(connection)
