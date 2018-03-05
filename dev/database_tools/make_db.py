from snappy import Manifold
from snappy.db_utilities import encode_torsion, encode_matrices, db_hash
import snappy.SnapPy
import os, sys, time, sqlite3, binascii, re, tarfile, gzip
from multiprocessing import Process, Lock, cpu_count
from DT_tools import *

snappy.SnapPy.matrix = snappy.SnapPy.SimpleMatrix
# Make the paths point in to the current directory
snappy.SnapPy.manifold_path = manifold_path = './'
snappy.SnapPy.closed_census_directory = os.path.join(manifold_path,
                                                     'ClosedCensusData')
snappy.SnapPy.link_directory = os.path.join(manifold_path, 'ChristyLinks')
snappy.SnapPy.link_archive = link_archive = os.path.join(manifold_path, 'ChristyLinks.tgz')
snappy.SnapPy.census_knot_archive = census_knot_archive = os.path.join(manifold_path, 'CensusKnots.tgz')
#snappy.SnapPy.Census_Morwen8 = tarfile.open(os.path.join(manifold_path, 'morwen8.tgz'), 'r:*')
snappy.SnapPy.Christy_links = tarfile.open(link_archive, 'r:*')
snappy.SnapPy.Census_Knots = tarfile.open(census_knot_archive, 'r:*')

from snappy.SnapPy import ObsOrientableCuspedCensus, ObsNonorientableCuspedCensus, ObsLinkExteriors, ObsCensusKnots, ObsOrientableClosedCensus, ObsNonorientableClosedCensus, ObsMorwenLinks
    
cusped_schema ="""
CREATE TABLE %s (
 id integer primary key,
 name text,
 cusps int,
 perm int,
 betti int,
 torsion blob,
 volume real,
 chernsimons real,
 tets int, 
 hash blob,
 triangulation blob)
"""

link_schema ="""
CREATE TABLE %s (
 id integer primary key,
 name text,
 cusps int,
 perm int,
 DT text,
 betti int,
 torsion blob,
 volume real,
 chernsimons real,
 tets int, 
 hash blob,
 triangulation blob)
"""

cusped_insert_query = """insert into %s
(name, cusps, perm, betti, torsion, volume, chernsimons, tets, hash, triangulation)
values ('%s', %s, %s, %s, X'%s', %s, %s, %s, X'%s', X'%s')"""

link_insert_query = """insert into %s
(name, cusps, perm, DT, betti, torsion, volume, chernsimons, tets, hash, triangulation)
values ('%s', %s, %s, '%s', %s, X'%s', %s, %s, %s, X'%s', X'%s')"""

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
 perm int,
 betti int,
 torsion blob,
 volume real,
 tets int, 
 hash blob,
 triangulation blob)
"""

nono_cusped_insert_query = """insert into %s
(name, cusps, perm, betti, torsion, volume, tets, hash, triangulation)
values ('%s', %s, %s, %s, X'%s', %s, %s, X'%s', X'%s')"""

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
    Create the empty tables for our basic manifold database.
    """
    for table in ['orientable_cusped_census','census_knots']:
        connection.execute(cusped_schema%table)
        connection.commit()
    for table in ['link_exteriors']:
        connection.execute(link_schema%table)
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

def make_views(connection):    
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
        # MC this looks rather pointless at the moment
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
    volume = float(mfld.volume())
    if mfld.is_orientable():
        try:
            chernsimons = float(mfld.chern_simons())
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

def tuples(isom):
    result = []
    for m in isom.cusp_maps():
        a, b, c, d = m[0,0], m[0,1], m[1,0], m[1,1]
        if a*d - b*c != 1:
            return False
        else:
            result.append( (a,b,c,d) )
    return result

def bytes_n_cobs(mfld):
    """
    Return a bytestring encoding of the manifold and a list of basis
    changes that convert the combinatorial basis of the decoded
    bytestring back to the original peripheral basis of the manifold,
    and a permutation to be applied to the cusp indices.
    """
    cobs = mfld.set_peripheral_curves('combinatorial', return_matrices=True)
    bytestring = mfld._to_bytes()
    encoded_perm = 0
    N = Manifold('empty')
    N._from_bytes(bytestring)
    N.set_peripheral_curves('combinatorial')
    mfld.set_peripheral_curves(cobs) # put it back the way it was
    isoms = mfld.isomorphisms_to(N)
    abcd = False
    while isoms:
        pick_one = isoms.pop()
        abcd = tuples(pick_one)
        if abcd:
            break
    if not abcd:
        print('No orientation preserving isometries????')
        return bytestring, cobs, 0
    perm = pick_one.cusp_images()
    for n in range(mfld.num_cusps()):
        a, b, c, d = abcd[n]
        cobs[perm[n]] = [[a,c],[b,d]]
        encoded_perm |= (n << (perm[n]<<2))
    return bytestring, cobs, encoded_perm

def insert_cusped_manifold(connection, table, mfld,
                           mfld_hash=db_hash,
                           is_link=False,
                           use_string=False,
                           DTcode=None):
    """
    Insert a cusped manifold into the specified table.
    """
    name = mfld.name()
    cusps = mfld.num_cusps()
    homology = mfld.homology()
    betti = homology.betti_number()
    divisors = [x for x in homology.elementary_divisors() if x > 0]
    torsion = binascii.hexlify(encode_torsion(divisors))
    volume = float(mfld.volume())
    if mfld.is_orientable():
        try:
            cs = float(mfld.chern_simons())
        except ValueError:
#            print 'Chern-Simons failed for %s'%name
            cs = 'NULL'
    tets = mfld.num_tetrahedra()
    use_cobs, triangulation = get_header(mfld, is_link, use_string)
    if use_cobs:
        bytestring, cobs, perm = bytes_n_cobs(mfld)
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
        triangulation += bytestring
    triangulation = binascii.hexlify(triangulation)
    try:
        hash_value = mfld_hash(mfld)
    except:
        print 'failed to hash %s (%s)'%(mfld, DTcode)
        hash_value = '00000000000000000000000000000000'
    if mfld.is_orientable():
        if DTcode is not None:
            query = link_insert_query%(
                table, name, cusps, perm, DTcode, betti,
                torsion, volume, cs, tets, hash_value, triangulation)
        else:
            query = cusped_insert_query%(
                table, name, cusps, perm, betti, torsion,
                volume, cs, tets, hash_value, triangulation)
    else:
        query = nono_cusped_insert_query%(
            table, name, cusps, perm, betti, torsion,
            volume, tets, hash_value, triangulation)
    connection.execute(query)

def copy_table_to_disk(connection, table, dbfile):
    print 'copying %s'%table
    connection.execute('attach "%s" as disk'%dbfile)
    connection.execute('begin transaction')
    connection.execute('insert into disk.%s select * from %s'%(table, table))
    connection.commit()
    connection.close()
    print 'finished %s'%table

def make_cusped(dbfile):
    connection = setup_db(":memory:")
    table = 'orientable_cusped_census'
    print 'making %s'%table
    for M in ObsOrientableCuspedCensus():
        M.set_name(M.name().split('(')[0])
        insert_cusped_manifold(connection, table, M, is_link=True)
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)

def make_cusped_nine(dbfile):
    import bz2, regina
    connection = setup_db(":memory:")
    table = 'orientable_cusped_census'
    print 'making 9 tetrahedra census'
    curr_len = 17661
    connection.execute("insert into %s (id, name) values (%d,'None')" % (table, curr_len))
    for line in bz2.BZ2File('cusped-9-or-bab.bz2'):
        parts = line.strip().split()
        name, vol, isosig = parts[0], parts[1], parts[6]
        M = snappy.Manifold(regina.NTriangulation(isosig).snapPea())
        M.set_name(name)
        M.set_peripheral_curves('shortest')
        assert abs(M.volume() - float(vol)) < 1e-10
        insert_cusped_manifold(connection, table, M, is_link=True)

    connection.execute("delete from %s where id = %d" % (table, curr_len))

    connection.commit()
    copy_table_to_disk(connection, table, dbfile)

def make_links(dbfile):
    dt_codes = dict(re.findall( '(\S+)\s+(\S+)$',
                                gzip.open('ChristyDT.gz').read(), re.MULTILINE))
    connection = setup_db(":memory:")
    table = 'link_exteriors'
    print 'making %s'%table
    for n in range(1, 6):
        for M in ObsLinkExteriors(n):
            M.set_name(M.name().split('(')[0])
            insert_cusped_manifold(connection, table, M, is_link=True, DTcode=dt_codes[M.name()])
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)
    
def make_census_knots(dbfile):
    connection = setup_db(":memory:")
    table = 'census_knots'
    print 'making %s'%table
    for M in ObsCensusKnots():
        M.set_name(M.name().split('(')[0])
        insert_cusped_manifold(connection, table, M, is_link=True)
    for line in gzip.open('knots_8_tet.gz'):
        knot_name, census_name, peripheral_curves = line.split('\t')
        M = snappy.Manifold(census_name)
        M.set_name(knot_name)
        M.set_peripheral_curves(eval(peripheral_curves))
        insert_cusped_manifold(connection, table, M, is_link=True)        
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)
    
def make_closed(dbfile):
    connection = setup_db(":memory:")
    table = 'orientable_closed_census'
    print 'making %s'%table
    for M in ObsOrientableClosedCensus():
        insert_closed_manifold(connection, table, M)
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)
    
def make_nono_cusped(dbfile):
    connection = setup_db(":memory:")
    table = 'nonorientable_cusped_census'
    for M in ObsNonorientableCuspedCensus():
        insert_cusped_manifold(connection, table, M)
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)

def make_nono_closed(dbfile):
    connection = setup_db(":memory:")
    table = 'nonorientable_closed_census'
    for M in ObsNonorientableClosedCensus():
        insert_closed_manifold(connection, table, M)
    connection.commit()
    copy_table_to_disk(connection, table, dbfile)

def make_indexes(dbfile, columns_to_index=['name']):
    """
    Add an index for each table on the name column. This is necessary for
    joins as well as looking up manifolds quickly.
    """
    connection = sqlite3.connect(dbfile)
    cur = connection.execute('select name from sqlite_master where type="table"')
    tables = [row[0] for row in cur.fetchall()]
    for table in tables:
        cols = [row[1] for row in cur.execute('PRAGMA table_info(%s)' % table).fetchall()]
        for col in columns_to_index:
            if col in cols:
                cur.execute('CREATE INDEX %s_%s_index ON %s (%s)' % (table, col, table, col))
    connection.close()
    
def setup_db(dbfile):
    if os.path.exists(dbfile):
        os.remove(dbfile)
    connection = sqlite3.connect(dbfile)
    create_manifold_tables(connection)
    make_views(connection)
    return connection

def make_basic_db():
    dbfile = 'manifolds.sqlite'
    setup_db(dbfile)
    workers = []
    workers.append(Process(target=make_nono_closed,  args=(dbfile,)))
    workers.append(Process(target=make_nono_cusped,  args=(dbfile,)))
    workers.append(Process(target=make_closed,       args=(dbfile,)))
    workers.append(Process(target=make_census_knots, args=(dbfile,)))
    workers.append(Process(target=make_links,        args=(dbfile,)))
    workers.append(Process(target=make_cusped,       args=(dbfile,)))
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join()
    make_indexes(dbfile)

def create_extended_tables(connection):
    """
    Create the empty tables for our big manifold database.
    """
    connection.execute(link_schema%'HT_links')
    connection.commit()
    
def make_extended_views(connection):
    connection.execute("""create view HT_links_view as
    select * from HT_links""")


def setup_extended_db(dbfile):
    if os.path.exists(dbfile):
        os.remove(dbfile)
    connection = sqlite3.connect(dbfile)
    create_extended_tables(connection)
    make_extended_views(connection)
    # Add the poor missing trefoil
    M = Manifold('3_1')
    M.set_name('K3a1')
    insert_cusped_manifold(connection, 'HT_links', M, is_link=True,
                           DTcode='cacbca')
    connection.commit()
    connection.close()

def triangulation_sort_key(M):
    sol_type = M.solution_type(enum=True)
    syms = -M.symmetry_group().order() if sol_type in {1, 2} else None
    return (sol_type, syms, M.num_tetrahedra())

def link_exterior_tri(DT_code, starts=4, random=4):
    tris = []
    for s in range(starts):
        M = Manifold('DT[%s]'%DT_code)
        for r in range(random+1):
            tris.append(M.copy())
            for i in range(r):
                M.randomize()
    return min(tris, key=triangulation_sort_key)
    
    
def make_HT_links(my_list, my_lock, next_lock, dbfile):
    # Used to make the next process wait before writing to disk.
    next_lock.acquire()
    table = 'HT_links'
    connection = sqlite3.connect(":memory:")
    connection.execute(link_schema%table)
    connection.commit()
    print os.getpid(), 'starting'
    connection.execute('begin transaction')
    for code, name in my_list:
        M = link_exterior_tri(code)
        M.set_name(name)
        insert_cusped_manifold(connection, table, M, is_link=True,
                               DTcode=code)
    connection.commit()
    print os.getpid(), 'waiting'

    # The first process doesn't need to wait, but the others
    # must acquire their lock before writing to disk.
    if my_lock is not None:
        my_lock.acquire()
        my_lock.release()
    print os.getpid(), 'copying'
    connection.execute("attach database '%s' as disk"%dbfile)
    connection.execute("begin transaction")
    connection.execute("""
    insert into disk.HT_links (
      name, cusps, perm, DT, betti, torsion, volume, chernsimons,
      tets, hash, triangulation
      )
    select
      name, cusps, perm, DT, betti, torsion, volume, chernsimons,
      tets, hash, triangulation
    from HT_links""")
    connection.commit()
    connection.close()
    # Now tell the next process to go ahead.
    next_lock.release()
    print os.getpid(), 'finished'

def make_extended_db():
    dbfile = 'more_manifolds.sqlite'
    procs = cpu_count()
    setup_extended_db(dbfile)
    links = all_links()
    totalsize = len(links)
    blocksize = 1 + totalsize/procs
    if procs == 4:
        chunks = [0, 58000, 100000, 143000, totalsize]
    elif procs == 8: # untested - please tune
        chunks = [0, 28000, 56000, 92000, 108000, 135000, 159000, 170000,
                  totalsize]
    else:
        chunks = [n*blocksize for n in range(procs+1)]
    blocksize = 1 + totalsize/procs
    locks = [None] + [Lock() for n in range(procs)]
    processes = [
        Process(target=make_HT_links,
                args=( links[chunks[n]:chunks[n+1]],
                       locks[n], locks[n+1], dbfile) )
        for n in range(procs) ]
    for process in processes:
        process.start()
    for process in processes:
        process.join()
    make_indexes(dbfile)
    
if __name__ == '__main__':
    make_basic_db()
    #make_extended_db()
    #make_census_knots('manifolds_old.sqlite')
    #make_cusped_nine('manifolds.sqlite')
    pass
