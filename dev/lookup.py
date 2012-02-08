from snappy import *
import sqlite3, bz2, re
from hashlib import md5
import array
from census import standard_hashes, appears_hyperbolic, find_hyperbolic_manifold_in_list

USE_COBS = 1 << 7
USE_STRING = 1 << 6
CUSP_MASK = 0x3f

def inflate_matrices(byteseq):
    """
    Convert a sequence of 4n bytes into a list of n 2x2 integer matrices.
    """
    m = array.array('b')
    m.fromstring(byteseq)
    return [ [ list(m[n:n+2]), list(m[n+2:n+4]) ]
             for n in range(0, len(m), 4) ]

class ManifoldDatabase:
    """
    Object for querying an sqlite3 database of manifolds.  Initialize
    with a database filename and a table name.  The table schema is
    required to include a text field called 'name' and a blob field
    called 'triangulation', which holds the result of M._to_bytes().
    """
    def __init__(self, dbfile='', table=''):
        self.connection = conn = sqlite3.connect(dbfile)
        cursor = conn.execute("pragma table_info('%s')"%table)
        rows = cursor.fetchall()
        self.schema = dict([(row[1],row[2].lower()) for row in rows])
        assert self.schema['name'] == 'text' and \
               self.schema['triangulation'] == 'blob', \
               'Not a valid Manifold table.'
        conn.row_factory = self._manifold_factory
        self.query = ('select name, triangulation from XXX where %s '
                      'order by %s limit %s').replace('XXX', table)

    def _manifold_factory(self, cursor, row):
        """
        Our queries will always return manifolds.
        """
        # Our rows contain only the name and triangulation fields.
        buf = row[1]
        header = ord(buf[0])
        use_cobs, use_string = header&USE_COBS, header&USE_STRING
        num_cusps = header&CUSP_MASK
        if use_cobs:
            cobs = inflate_matrices(buf[1:4*num_cusps + 1])
            triangulation = bytes(buf[4*num_cusps +1:])
        else:
            triangulation = bytes(buf[1:])
        M = Manifold('empty')
        if use_string:
            M._from_string(triangulation)
        else:
            M._from_bytes(triangulation)
        M.set_name(row[0])
        if use_cobs:
            M.set_peripheral_curves(cobs)
        return M

    def keys(self):
        return self.schema.keys()
    
    def find(self, where='0=1', order_by='id', limit=25):
        """
        Find up to limit manifolds in the census satisfying the
        where clause, ordered by the order_by clause.
        """
        cursor = self.connection.execute(self.query%(where, order_by, limit))
        return cursor.fetchall()
    
    def find_by_volume(self, vol, tolerance, limit=25):
        """
        Find up to limit manifolds whose volume is equal to vol to
        within the specified tolerance, ordered by volume.
        """
        where = 'volume > %g and volume < %g'%(vol-tolerance, vol+tolerance)
        order_by = 'volume'
        return self.find(where=where, order_by=order_by)

    def siblings(self, mfld):
        """
        Return all manifolds in the census which have the same hash.
        """
        hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
        return self.find(where="hash = X'%s'"%hash)

    def identify(self, mfld):
        return find_hyperbolic_manifold_in_list(mfld, self.siblings(mfld))

    def __getitem__(self, index):
        try:
            where = 'id=%d' % (index + 1) 
        except TypeError:
            where = 'name="' + index + '"'

        matches = self.find(where)
        if len(matches) != 1:
            raise IndexError
        return matches[0]
    
DB = ManifoldDatabase(dbfile='manifolds.sqlite', table='cusped_census')

def test_census_database():
    L = ManifoldDatabase(dbfile='manifolds.sqlite', table='cusped_census')
    for M in CensusKnots():
        print M, L.identify(M)

def SmallHTWKnots():
    for census in [ AlternatingKnotExteriors(), NonalternatingKnotExteriors()]:
        for M in census:
            if re.match('12', M.name()):
                break
            yield M

def test_link_database():
    # Broken, as the tables do not exist yet.
    #print len([M for M in SmallHTWKnots()]), len([M for M in LinkExteriors(1)])
    L = ManifoldVerboseDatabase(dbfile='links.sqlite', table='census')
    K = ManifoldVerboseDatabase(dbfile='new_knots.sqlite', table='census')
    for census, db in [ (SmallHTWKnots(), L), (LinkExteriors(1), K) ]:
        non_hyp, missing = [], []
        count = 0
        for M in census:
            if M.volume() > 0.5:
                N = db.identify(M)
                if N == None:
                    missing.append(M)
            else:
                non_hyp.append(M)

            count += 1

        print count, len(db.find('cusps=1', limit=1000)), missing, len(non_hyp), non_hyp


if __name__ == '__main__':
    test_census_database()
    #ans = test_link_database()
