from snappy import *
import os
import sqlite3
import binascii
from hashlib import md5
import array
from census import *

cusped_schema ="""
CREATE TABLE %s (
 id integer primary key,
 name text,
 volume real,
 chernsimons real,
 hash blob,
 triangulation blob)
"""

cusped_insert_query = """insert into %s
(name, volume, chernsimons, hash, triangulation)
values ('%s', %s, %s, X'%s', X'%s')"""

USE_COBS = 1 << 7
USE_STRING = 1 << 6
epsilon = 0.000001

def flatten_matrices(matrices):
    """
    Convert a list of 2x2 integer matrices into a sequence of bytes.
    """
    # The tricky point here is converting signed integers to bytes.
    return bytes(array.array('b', sum(sum(matrices,[]),[])).tostring())
    # NOTE: tostring is deprecated in python3, but for now
    # it does the same thing as tobytes.
    
def create_manifold_db(connection):
    """
    Create the empty tables for our manifold database.
    """
    for table in ['orientable_cusped_census',
                  'link_exteriors',
                  'census_knots']:
        connection.execute(cusped_schema%table)
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
    if use_string:
        header |= USE_STRING
    elif is_link or ambiguity_exists(mfld):
        header |= USE_COBS
    return header|USE_COBS, bytes(bytearray([header]))
    
def insert_cusped_manifold(connection, table, mfld,
                           is_link=False,
                           use_string=False):
    name = mfld.name()
    volume = mfld.volume()
    try:
        cs = mfld.chern_simons()
    except ValueError:
        print 'Chern-Simons failed for %s'%name
        cs = 'NULL'
    use_cobs, triangulation = get_header(mfld, is_link, use_string)
    if use_cobs:
        cobs = mfld.set_peripheral_curves('combinatorial')
        mfld.set_peripheral_curves(cobs) # undo the basis change
        try:
            triangulation += flatten_matrices(cobs)
        except OverflowError:
            #fall back to the verbose string record
            header = mfld.num_cusps() | USE_STRING
            triangulation = bytes(bytearray([header]))
            use_string = True
    if use_string:
        triangulation += mfld._to_string()
    else:
        triangulation += mfld._to_bytes()
    triangulation = binascii.hexlify(triangulation)
    hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
    query = cusped_insert_query%(table, name, volume, cs, hash, triangulation)
    try:
        connection.execute(query)
    except:
        print query

def make_cusped(connection):
    table = 'orientable_cusped_census'
    for M in OrientableCuspedCensus():
        insert_cusped_manifold(connection, table, M)
    connection.commit()

def make_links(connection):
    table = 'link_exteriors'
    for n in range(1, 6):
        for M in LinkExteriors(n):
            insert_cusped_manifold(connection, table, M,
                                   is_link=True)
    connection.commit()

def make_census_knots(connection):
    table = 'census_knots'
    for M in CensusKnots():
        insert_cusped_manifold(connection, table, M,
                               is_link=True)
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
    make_links(connection)
    make_census_knots(connection)
