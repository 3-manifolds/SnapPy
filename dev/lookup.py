from snappy import *
import sqlite3, bz2, re
from hashlib import md5
from census import standard_hashes, appears_hyperbolic, find_hyperbolic_manifold_in_list
            
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
        M = Manifold('empty')
        M._from_bytes(bytes(row[1]))
        M.set_name(row[0])
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

class ManifoldVerboseDatabase(ManifoldDatabase):
    def _manifold_factory(self, cursor, row):
        """
        Our queries will always return manifolds.
        """
        M = Manifold('empty')
        M._from_string(bz2.decompress(row[1]))
        M.set_name(row[0])
        return M   
            
    
#DB = ManifoldDatabase(dbfile='census.sqlite', table='census')
#DL = ManifoldVerboseDatabase(dbfile='links.sqlite', table='census')
#print DL.find('1=1')


def test_census_database():
    L = ManifoldDatabase(dbfile='census.sqlite', table='census')
    for M in CensusKnots():
        print M, L.identify(M)

def SmallHTWKnots():
    for census in [NonalternatingKnotExteriors(), AlternatingKnotExteriors()]:
        for M in census:
            if re.match('12', M.name()):
                break
            yield M

    
    
def test_link_database():
    #print len([M for M in SmallHTWKnots()]), len([M for M in LinkExteriors(1)])
    L = ManifoldVerboseDatabase(dbfile='links.sqlite', table='census')
    K = ManifoldVerboseDatabase(dbfile='new_knots.sqlite', table='census')
    for census, db in [ (SmallHTWKnots(), K), (LinkExteriors(1), L) ]:
        non_hyp, missing = [], []
        count = 0
        for M in census:
            if appears_hyperbolic(M):
                N = db.identify(M)
                if N == None:
                    missing.append(N)
            else:
                non_hyp.append(N)

            count += 1

        print count, missing, len(non_hyp), non_hyp

if __name__ == '__main__':
    # test_census_database()
    test_link_database()
