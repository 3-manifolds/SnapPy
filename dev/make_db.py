from snappy import *
from db_utils import encode_torsion, encode_matrices, db_hash
import os, sys, time
import sqlite3
import binascii
import re

cusped_schema ="""
CREATE TABLE %s (
 id integer primary key,
 name text,
 cusps int,
 betti int,
 torsion blob,
 volume real,
 chernsimons real,
 tets int, 
 hash blob,
 triangulation blob)
"""

cusped_insert_query = """insert into %s
(name, cusps, betti, torsion, volume, chernsimons, tets, hash, triangulation)
values ('%s', %s, %s, X'%s', %s, %s, %s, X'%s', X'%s')"""

closed_schema ="""
CREATE TABLE %s (
 id integer primary key,
 cusped text,
 m int,
 l int,
 betti int,
 torsion blob,
 volume real,
 chernsimons real,
 hash blob
)
"""

closed_insert_query = """insert into %s
(cusped, m, l, betti, torsion, volume, chernsimons, hash)
values ('%s', %d, %d, %d, X'%s', %s, %s, X'%s')"""

nono_cusped_schema ="""
CREATE TABLE %s (
 id integer primary key,
 name text,
 cusps int,
 betti int,
 torsion blob,
 volume real,
 tets int, 
 hash blob,
 triangulation blob)
"""

nono_cusped_insert_query = """insert into %s
(name, cusps, betti, torsion, volume, tets, hash, triangulation)
values ('%s', %s, %s, X'%s', %s, %s, X'%s', X'%s')"""

nono_closed_schema ="""
CREATE TABLE %s (
 id integer primary key,
 cusped text,
 m int,
 l int,
 betti int,
 torsion blob,
 volume real,
 hash blob
)
"""

nono_closed_insert_query = """insert into %s
(cusped, m, l, betti, torsion, volume, hash)
values ('%s', %d, %d, %d, X'%s', %s, X'%s')"""

USE_COBS = 1 << 7
USE_STRING = 1 << 6
epsilon = 0.000001

closed_re = re.compile('(.*)\((.*),(.*)\)')

def create_manifold_tables(connection):
    """
    Create the empty tables for our manifold database.
    """
    for table in ['orientable_cusped_census',
                  'link_exteriors',
                  'census_knots']:
        connection.execute(cusped_schema%table)
        connection.commit()
    for table in ['orientable_closed_census']:
        connection.execute(closed_schema%table)
        connection.commit()
    for table in ['nonorientable_cusped_census']:
        connection.execute(nono_cusped_schema%table)
        connection.commit()
    for table in ['nonorientable_closed_census']:
        connection.execute(nono_closed_schema%table)
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
    triangulation string. (In this case, change of basis is not
    needed.) The rest of the byte holds the number of cusps.
    """
    header = mfld.num_cusps()
    if use_string:
        header |= USE_STRING
    elif is_link or ambiguity_exists(mfld):
        header |= USE_COBS
    return header|USE_COBS, bytes(bytearray([header]))

def insert_closed_manifold(connection, table, mfld):
    """
    Insert a closed manifold into the specified table.
    """
    cusped, m, l = closed_re.match(repr(mfld)).groups()
    homology = mfld.homology()
    betti = homology.betti_number()
    divisors = [x for x in homology.elementary_divisors() if x > 0]
    torsion = binascii.hexlify(encode_torsion(divisors))
    volume = mfld.volume()
    if mfld.is_orientable():
        try:
            chernsimons = mfld.chern_simons()
        except:
            chernsimons = 'NULL'
    hash = db_hash(mfld)
    if mfld.is_orientable():
        query = closed_insert_query%(
            table, cusped, int(m), int(l), int(betti),
            torsion, volume, chernsimons, hash)
    else:
        query = nono_closed_insert_query%(
            table, cusped, int(m), int(l), int(betti),
            torsion, volume, hash)
    connection.execute(query)
    
def insert_cusped_manifold(connection, table, mfld,
                           is_link=False,
                           use_string=False):
    """
    Insert a cusped manifold into the specified table.
    """
    name = mfld.name()
    cusps = mfld.num_cusps()
    homology = mfld.homology()
    betti = homology.betti_number()
    divisors = [x for x in homology.elementary_divisors() if x > 0]
    torsion = binascii.hexlify(encode_torsion(divisors))
    volume = mfld.volume()
    if mfld.is_orientable():
        try:
            cs = mfld.chern_simons()
        except ValueError:
            print 'Chern-Simons failed for %s'%name
            cs = 'NULL'
    tets = mfld.num_tetrahedra()
    use_cobs, triangulation = get_header(mfld, is_link, use_string)
    if use_cobs:
        cobs = mfld.set_peripheral_curves('combinatorial')
        mfld.set_peripheral_curves(cobs) # undo the basis change
        try:
            triangulation += encode_matrices(cobs)
        except OverflowError:
            #fall back to the verbose string record
            header = mfld.num_cusps() | USE_STRING
            triangulation = bytes(bytearray([header]))
            use_string = True
    if use_string:
        triangulation += mfld.without_hyperbolic_structure()._to_string()
    else:
        triangulation += mfld._to_bytes()
    triangulation = binascii.hexlify(triangulation)
    hash = db_hash(mfld)
    if mfld.is_orientable():
        query = cusped_insert_query%(
            table, name, cusps, betti, torsion,
            volume, cs, tets, hash, triangulation)
    else:
        query = nono_cusped_insert_query%(
            table, name, cusps, betti, torsion,
            volume, tets, hash, triangulation)
    connection.execute(query)

def make_cusped(connection):
    table = 'orientable_cusped_census'
    for M in OrientableCuspedCensus():
        M.set_name(M.name().split('(')[0])
        insert_cusped_manifold(connection, table, M)
    connection.commit()
    # This index makes it fast to join this table on its name column.
    # Without the index, the join is very slow.
    connection.execute(
        'create index o_cusped_by_name on orientable_cusped_census (name)')

def make_links(connection):
    table = 'link_exteriors'
    for n in range(1, 6):
        for M in LinkExteriors(n):
            M.set_name(M.name().split('(')[0])
            insert_cusped_manifold(connection, table, M,
                                   is_link=True)
    connection.commit()

def make_morwen_links(connection):
    table = 'morwen_links'
    for n in range(1, 8):
        m=1
        for M in MorwenLinks(n):
            M.set_name(M.name().split('(')[0])
            print '%s %s %s'%(n, m, M.name())
            m += 1
            insert_cusped_manifold(connection, table, M,
                                   is_link=True)
    connection.commit()

def make_census_knots(connection):
    table = 'census_knots'
    for M in CensusKnots():
        M.set_name(M.name().split('(')[0])
        insert_cusped_manifold(connection, table, M,
                               is_link=True)
    connection.commit()

def make_closed(connection):
    table = 'orientable_closed_census'
    for M in OrientableClosedCensus():
        insert_closed_manifold(connection, table, M)
    connection.commit()

def make_nono_cusped(connection):
    table = 'nonorientable_cusped_census'
    for M in NonorientableCuspedCensus():
        insert_cusped_manifold(connection, table, M)
    connection.commit()

def make_nono_closed(connection):
    table = 'nonorientable_closed_census'
    for M in NonorientableClosedCensus():
        insert_closed_manifold(connection, table, M)
    connection.commit()

if __name__ == '__main__':
    dbfile = 'manifolds.sqlite'
    if os.path.exists(dbfile):
        os.remove(dbfile)
    connection = sqlite3.connect(dbfile)
    create_manifold_tables(connection)
    make_closed(connection)
    make_cusped(connection)
    make_links(connection)
    make_census_knots(connection)
    make_nono_cusped(connection)
    make_nono_closed(connection)
    # There are two reasons for using views.  One is that views
    # are read-only, so we have less chance of deleting our data.
    # The second is that they allow joins to be treated as if they
    # were tables, which we need for the closed census.
    connection.execute("""create view orientable_cusped_view as
    select * from orientable_cusped_census""")
    connection.execute("""create view link_exteriors_view as
    select * from link_exteriors""")
    connection.execute("""create view census_knots_view as
    select * from census_knots""")
    connection.execute("""create view nonorientable_cusped_view as
    select * from nonorientable_cusped_census""")
    connection.execute("""create view orientable_closed_view as
    select a.id, b.name, a.m, a.l, a.betti, a.torsion, a.volume,
    a.chernsimons, a.hash, b.triangulation
    from orientable_closed_census a
    left join orientable_cusped_census b
    on a.cusped=b.name""")
    connection.execute("""create view nonorientable_closed_view as
    select a.id, b.name, a.m, a.l, a.betti, a.torsion, a.volume,
    a.hash, b.triangulation
    from nonorientable_closed_census a
    left join nonorientable_cusped_census b
    on a.cusped=b.name""")
    connection.close()
