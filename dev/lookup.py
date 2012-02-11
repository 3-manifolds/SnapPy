from snappy import *
import sqlite3, bz2, re
from hashlib import md5
import array
from census import standard_hashes, appears_hyperbolic, find_hyperbolic_manifold_in_list

# This module uses a single sqlite3 database with multiple tables.
# The path to the database file is specified at the module level.
database_path = 'manifolds.sqlite'

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

class ManifoldTable:
    """
    Object for querying an sqlite3 table of cusped manifolds and
    iterating through its rows.  Initialize with the table name.  The
    table schema is required to include a text field called 'name' and
    a blob field called 'triangulation'.  The blob holds the result of
    M._to_bytes() or M._to_string(), optionally preceded by a change
    of basis matrix for the peripheral curves.  The structure of the
    blob is determined by its first byte.
    """
    def __init__(self, table='', betti=None, num_cusps=None):
        self.table = table
        self.connection = conn = sqlite3.connect(database_path)
        cursor = conn.execute("pragma table_info('%s')"%table)
        rows = cursor.fetchall()
        self.schema = dict([(row[1],row[2].lower()) for row in rows])
        assert self.schema['name'] == 'text' and \
               self.schema['triangulation'] == 'blob', \
               'Not a valid Manifold table.'
        cursor = conn.execute("select count(*) from %s"%self.table)
        self._length = cursor.fetchone()[0]
        conn.row_factory = self._manifold_factory
        filters = []
        if betti:
            filters.append('betti=%d'%betti)
        if num_cusps:
            filters.append('cusps=%d'%num_cusps)
        if filters:
            self.filter_clause = ' and %s '%' and '.join(filters)
        else:
            self.filter_clause = ''
        self.select = 'select name, triangulation from %s '%(table)

    def __len__(self):
        return self._length
        
    def __iter__(self):
        return self.connection.execute(
            'select name, triangulation from %s %s'%(self.table,
                                                     self.filter_clause)
            )

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step:
                raise IndexError('Slices with steps are not supported.')
            start, stop = index.start, index.stop
            if (isinstance(start, (float, type(None)) )
                and
                isinstance(stop, (float, type(None)) ) ):
                # Slice by volume.
                limits = []
                if start:
                    limits.append('volume >= %g'%start)
                if stop:
                    limits.append('volume < %g'%stop)
                stop_clause = 'volume < %g'%stop if stop else ''    
                query = (self.select + 'where %s '%' and '.join(limits))
                return self.connection.execute(query)
            elif (isinstance(start, (int, type(None)) )
                  and
                  isinstance(stop, (int, type(None)) ) ):
                if start and start < 0:
                    start = self._length + start
                if stop and stop < 0:
                    stop = self._length + stop
                if self.filter_clause == '':
                    # With no filter we can slice by the id field;
                    start = 0 if start is None else start 
                    limit_clause = 'limit %d'%(stop - start) if stop else ''
                    query = (self.select + 'where id >= %d  %s '%(
                                 start + 1,
                                 limit_clause))
                    return self.connection.execute(query)
                # otherwise we just trash the rows at the beginning. :^(
                else:
                    query = (self.select + 'where 0=0 %s limit %d'%(
                                 self.filter_clause,
                                 stop))
                    cursor = self.connection.execute(query)
                    if start:
                        cursor.row_factory = lambda x, y : None
                        cursor.fetchmany(start)
                        cursor.row_factory = self._manifold_factory
                    return cursor
            else:
                raise IndexError(
                    'Use two ints or two floats for start and stop.')
        elif isinstance(index, int):
            matches = self.find('id=%d'%(index + 1))
            if len(matches) != 1:
                raise IndexError('Manifold index is out of bounds')
        elif isinstance(index, str):
            matches = self.find("name='%s'"%index)
            if len(matches) != 1:
                raise KeyError('Did not find a manifold named %s.'%index)
        else:
            raise IndexError('%s is not a valid index type for manifolds.'%
                             type(index))
        return matches[0]
    
    def _manifold_factory(self, cursor, row):
        """
        Factory for "select name, triangulation" queries.
        Returns a Manifold.
        """
        buf = row[1]
        header = ord(buf[0])
        use_cobs, use_string = header&USE_COBS, header&USE_STRING
        num_cusps = header&CUSP_MASK
        M = Manifold('empty')
        if use_string:
            M._from_string(buf[1:])
        else:
            M._from_bytes(bytes(buf[4*num_cusps +1:]))
            if use_cobs:
                cobs = inflate_matrices(buf[1:4*num_cusps + 1])
                M.set_peripheral_curves('combinatorial')
                M.set_peripheral_curves(cobs)
        M.set_name(row[0])
        return M

    def filtered_by(self, betti=None, num_cusps=None):
        return ManifoldTable(self.table, betti=betti, num_cusps=num_cusps)
    
    def keys(self):
        """
        Return the list of column names for this manifold table.
        """
        return self.schema.keys()
    
    def find(self, where, order_by='id', limit=25):
        """
        Return a list of up to limit manifolds stored in this manifold
        table, satisfying the where clause, and ordered by the order_by
        clause.  If limit is None, all matching manifolds are returned.
        A where clause is required.
        """
        if limit is None:
            suffix = 'where %s order by %s'%(where, order_by)
        else:
            suffix = 'where %s order by %s limit %d'%(where, order_by, limit)
        cursor = self.connection.execute(self.select + suffix)
        return cursor.fetchall()

    def siblings(self, mfld):
        """
        Return all manifolds in the census which have the same hash value.
        """
        hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
        return self.find(where="hash = X'%s'"%hash)

    def identify(self, mfld):
        """
        Return the manifold in the able which is isometric to the
        argument, if it exists.
        """
        return find_hyperbolic_manifold_in_list(mfld,
                                                self.siblings(mfld))

# Instantiate our tables.
OrientableCuspedDB = ManifoldTable(table='orientable_cusped_census')
LinkExteriorDB = ManifoldTable(table='link_exteriors')
CensusKnotsDB = ManifoldTable(table='census_knots')
                                  
# Test routines.
def test_census_database():
    L = OrientableCuspedDB
    for M in CensusKnots():
        print M, L.identify(M)

def SmallHTWKnots():
    for census in [ AlternatingKnotExteriors(), NonalternatingKnotExteriors()]:
        for M in census:
            if re.match('12', M.name()):
                break
            yield M

def test_link_database():
    # Broken at the moment
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

def manifolds_match(M, N):
    isoms = M.is_isometric_to(N, True)
    for i in isoms:
        n = M.num_cusps()
        if i.cusp_images() == range(n):
            if [m for m in i.cusp_maps() if m.tolist() != [[1,0],[0,1]]] == []:
                return True
    return False

def test():
    pairs = [ (CensusKnots(), CensusKnotsDB),
              (OrientableCuspedCensus(), OrientableCuspedDB)]
    for census, db in pairs:
        for M in census:
            if M.name()[0] == 't':
                break
            N = db.identify(M)
            assert repr(M) == repr(N)
            G, H = M.fundamental_group(), N.fundamental_group()
            if G.relators() != H.relators() or G.peripheral_curves() != H.peripheral_curves():
                print M

if __name__ == '__main__':
    test()
