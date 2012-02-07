from snappy import *
import sqlite3, bz2
from hashlib import md5
from census import standard_hashes
            
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
        print bz2.decompress(bytes(row[1]))
        return Manifold('m004')
        # Our rows contain only the name and triangulation fields.
        #M = Manifold('empty')
        #M._from_string()
        #M.set_name(row[0])
        #return M   
            
            
    
#DB = ManifoldDatabase(dbfile='census.sqlite', table='census')
DL = ManifoldVerboseDatabase(dbfile='census.sqlite', table='census')
print DL.find('1=1')
