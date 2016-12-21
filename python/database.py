from __future__ import print_function
# NB: this module uses the Manifold class from the SnapPy extension
# module. After the Manifold class has been defined, the SnapPy
# extension module sets database.Manifold = Manifold .
from .db_utilities import decode_torsion, decode_matrices, db_hash
from spherogram.codecs import DTcodec
import sqlite3, re, os, random

try:
    unicode
    byte_to_int = ord
except NameError: # Python 3
    byte_to_int = int

try:
    import sage.all
    def is_int(slice):
        return isinstance(slice, (sage.all.Integer,int))
        
    def is_int_or_none(slice):
        return isinstance(slice, (sage.all.Integer,int, type(None)))

    def is_float_or_none(slice):
        return isinstance(slice, (float, sage.all.RealDoubleElement,
                                  sage.rings.real_mpfr.RealNumber, type(None)))

except ImportError:
    def is_int(slice):
        return isinstance(slice, int)
    
    def is_int_or_none(slice):
        return isinstance(slice, (int, type(None)))

    def is_float_or_none(slice):
        return isinstance(slice, (float, type(None)))

# This module uses sqlite3 databases with multiple tables.
# The path to the database file is specified at the module level.
from .manifolds import __path__ as manifolds_paths
manifolds_path = manifolds_paths[0]
database_path = os.path.join(manifolds_path, 'manifolds.sqlite')
# Temporary - should get this from preferences.
alt_database_path = os.path.join(manifolds_path, 'more_manifolds.sqlite')
platonic_database_path = os.path.join(manifolds_path, 'platonic_manifolds.sqlite')

USE_COBS = 1 << 7
USE_STRING = 1 << 6
CUSP_MASK = 0x3f

split_filling_info = re.compile('(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))*$)')

def mfld_hash(manifold):
    """
    We cache the hash to speed up searching for one manifold in
    multiple tables.
    """
    if 'db_hash' not in manifold._cache:
        manifold._cache['db_hash'] = db_hash(manifold)
    return manifold._cache['db_hash']

class ManifoldTable(object):
    """
    Iterator for cusped manifolds in an sqlite3 table of manifolds.

    Initialize with the table name.  The table schema is required to
    include a text field called 'name' and a blob field called
    'triangulation'.  The blob holds the result of M._to_bytes() or
    M._to_string(), optionally preceded by a change of basis matrix
    for the peripheral curves.  The structure of the blob is
    determined by its first byte.

    Both mapping from the manifold name, and lookup by index are
    supported.  Slicing can be done either by numerical index or by
    volume.

    The __contains__ method is supported, so M in T returns True if M
    is isometric to a manifold in the table T.  The method
    T.identify(M) will return the matching manifold from the table.
    """
    # basic select clause.  Can be overridden, e.g. to add additional columns
    _select = 'select name, triangulation, perm from %s '

    def __init__(self, table='', db_path=database_path,
                 mfld_hash=mfld_hash, **filter_args):
        self._table = table
        self.mfld_hash = mfld_hash
        self._connection = sqlite3.connect(db_path)
        self._connection.row_factory = self._manifold_factory
        # Sometimes we need a connection without the row factory
        self._connection2 = conn = sqlite3.connect(db_path)
        self._set_schema()
        self._check_schema()
        self._configure(**filter_args)
        self._get_length()
        self._get_max_volume()
        self._select = self._select%table

    def _set_schema(self):
        conn, table = self._connection2, self._table
        cursor = conn.execute("pragma table_info('%s')" % table)
        rows = cursor.fetchall()
        self.schema = dict([(row[1],row[2].lower()) for row in rows])

    def _check_schema(self): 
        assert (self.schema['name'] == 'text' and 
                self.schema['triangulation'] == 'blob'
                ), 'Not a valid Manifold table.'

    @property
    def filter(self):
        return self._filter

    def _get_length(self):
        where_clause = 'where ' + self._filter if self._filter else '' 
        length_query = 'select count(*) from %s %s' % (self._table,
                                                       where_clause)
        cursor = self._connection2.execute(length_query)
        self._length = cursor.fetchone()[0]

    def _get_max_volume(self):
        where_clause = 'where ' + self._filter if self._filter else '' 
        vol_query = 'select max(volume) from %s %s' % (self._table,
                                                       where_clause)
        cursor = self._connection2.execute(vol_query)
        self._max_volume = cursor.fetchone()[0]
        
    def _configure(self, **kwargs):
        """
        Set up the filter.
        """
        conditions = []

        if 'filter' in kwargs:
            conditions.append(kwargs['filter'])
        if 'betti' in kwargs:
            conditions.append('betti=%d ' % kwargs['betti'])
        if 'num_cusps' in kwargs:
            conditions.append('cusps=%d ' % kwargs['num_cusps'])
        if 'cusps' in kwargs:
            conditions.append('cusps=%d ' % kwargs['cusps'])
        if 'num_tets' in kwargs:
            conditions.append('tets=%d ' % kwargs['num_tets'])
        if 'tets' in kwargs:
            conditions.append('tets=%d ' % kwargs['tets'])
        self._filter = ' and '.join(conditions)
         
    def __repr__(self):
        class_name = self.__class__.__name__
        if self._filter == '':
            return '%s without filters'%class_name
        else:
            return '%s with filter: %s'%(class_name, self._filter)
        
    def __call__(self, **kwargs):
        return self.__class__(**kwargs)
    
    def __len__(self):
        return self._length
        
    def __iter__(self):
        query = self._select
        if self._filter:
            query += ' where %s '%self._filter
        return self._connection.execute(query)

    def __contains__(self, mfld):
        try:
            M = self.identify(mfld)
            # duck test
            return M.num_tetrahedra() > 0
        except:
            return False
        
    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step:
                raise IndexError('Slices with steps are not supported.')
            start, stop = index.start, index.stop
            if is_float_or_none(start) and is_float_or_none(stop):
                # Slice by volume.
                conditions = []
                if self._filter:
                    conditions.append(self._filter)
                if start:
                    conditions.append('volume >= %f' % start)
                if stop:
                    conditions.append('volume < %f' % stop)
                where_clause = ' and '.join(conditions)
                if where_clause:
                    where_clause = 'where ' + where_clause
                query = (self._select + where_clause)
                return self._connection.execute(query).fetchall()
            elif (is_int_or_none(start) and is_int_or_none(stop)):
                if start and start < 0:
                    start = int(self._length + start)
                if stop and stop < 0:
                    stop = int(self._length + stop)
                if self._filter == '':
                    # With no filter we can slice by the id field;
                    start = 0 if start is None else start 
                    limit_clause = ' limit %d '%(stop - start) if stop else ''
                    query = (self._select + 'where id >= %d  %s ' % (
                                 start + 1,
                                 limit_clause))
                    return self._connection.execute(query)
                # otherwise we just trash the rows at the beginning. :^(
                else:
                    limit_clause = ' limit %d'%stop if stop else ''
                    query = (self._select + 'where %s %s'%(
                                 self._filter, limit_clause))
                    cursor = self._connection.execute(query)
                    if start:
                        cursor.row_factory = lambda x, y : None
                        cursor.fetchmany(start)
                        cursor.row_factory = self._manifold_factory
                    return cursor
            else:
                raise IndexError(
                    'Use two ints or two floats for start and stop.')
        elif is_int(index):
            if self.filter == '':
                if index < 0:
                    matches = self.find('id=%d'%(index + self._length + 1))
                else:
                    matches = self.find('id=%d'%(index + 1))
            else:
                matches = list(self[index:index+1])
            if len(matches) != 1:
                raise IndexError('Manifold index is out of bounds')
        elif isinstance(index, str):
            matches = self.find("name='%s'"%index)
            if len(matches) != 1:
                raise KeyError('The manifold %s was not found.'%index)
        else:
            raise IndexError('%s is not a valid index type for manifolds.'%
                             type(index))
        return matches[0]
    
    def _manifold_factory(self, cursor, row, M=None):
        """
        Factory for "select name, triangulation" queries.
        Returns a Manifold.
        """
        if M is None:
            M = Manifold('empty')
        buf = bytes(row[1])
        header = byte_to_int(buf[0])
        use_cobs, use_string = header&USE_COBS, header&USE_STRING
        num_cusps = header&CUSP_MASK
        if use_string:
            M._from_string(buf[1:])
        else:
            M._from_bytes(buf[4*num_cusps + 1:])
            if use_cobs:
                cobs = decode_matrices(buf[1:4*num_cusps + 1])
                M.set_peripheral_curves('combinatorial')
                M.set_peripheral_curves(cobs)
        self._finalize(M, row)
        return M

    def _finalize(self, M, row):
        """
        Give the manifold a name and make last-minute adjustments to
        the manifold before it leaves the factory, e.g. reordering the
        cusps.  Override this method for custom manifold production.
        """
        M.set_name(row[0])
        num = M.num_cusps()
        encoded_perm = row[2]
        if encoded_perm:
            perm = [(encoded_perm >> (n<<2)) & 0xf for n in range(num)]
            M._reindex_cusps(perm)
        # This seems to be necessary to make the triangulation
        # structure consistent.
        #M.dehn_fill([(0,0)]*num)

    def _one_manifold(self, name, M):
        """
        Inflates the given empty Manifold with the table manifold
        with the specified name.
        """
        cursor = self._connection2.execute(self._select + "where name='" + name + "'")
        rows = cursor.fetchall()
        if len(rows) != 1:
            raise KeyError('The manifold %s was not found.'%name)
        return self._manifold_factory(None, rows[0], M)
                
    def keys(self):
        """
        Return the list of column names for this manifold table.
        """
        return self.schema.keys()
    
    def find(self, where, order_by='id', limit=None):
        """
        Return a list of up to limit manifolds stored in this table,
        satisfying the where clause, and ordered by the order_by
        clause.  If limit is None, all matching manifolds are
        returned.  The where clause is a required parameter.
        """
        where_clause = where
        if self._filter:
            where_clause += ' and ' + self._filter
        if limit is None:
            suffix = 'where %s order by %s'%(where_clause, order_by)
        else:
            suffix = 'where %s order by %s limit %d'%(
                where_clause, order_by, limit)
        cursor = self._connection.execute(self._select + suffix)
        return cursor.fetchall()

    def siblings(self, mfld):
        """
        Return all manifolds in the census which have the same hash value.
        """
        return self.find("hash = X'%s'"%self.mfld_hash(mfld))

    def identify(self, mfld, extends_to_link=False):
        """
        Look for a manifold in this table which is isometric to the
        argument.

        Return the matching manifold, if there is one which SnapPea
        declares to be isometric.

        Return False if no manifold in the table has the same hash.

        Return None in all other cases (for now).

        If the flag "extends_to_link" is True, requires that the isometry
        sends meridians to meridians.   If the input manifold is closed
        this will result in no matches being returned.  
        """
        if hasattr(mfld, 'volume'):
            if mfld.volume() > self._max_volume + 1:
                return False 

        if extends_to_link and not (True in mfld.cusp_info('complete?')):
            return False

        sibs = self.siblings(mfld)
        if len(sibs) == 0:
            return False # No hash values match
        
        mfld = mfld.copy()
        mflds = [mfld]
        for i in range(4):
            mfld = mfld.copy()
            mfld.randomize()
            mflds.append(mfld)
        
        # Check for isometry
        for mfld in mflds:
            for N in sibs:
                try:
                    if not extends_to_link:
                        if mfld.is_isometric_to(N):
                            return N
                    else:
                        isoms = mfld.is_isometric_to(N, True)
                        if True in [i.extends_to_link() for i in isoms]:
                            return N
                except RuntimeError:
                    pass

        mfld = mflds[0]
        # Check for identical triangulations.
        if (not False in mfld.cusp_info('is_complete')) and not extends_to_link:
            for n in range(100):
                for N in sibs:
                    if mfld == N:
                        return N
                mfld.randomize()
        
        return None

    def random(self):
        return self[random.randrange(len(self))]
        
class ClosedManifoldTable(ManifoldTable):

    _select = 'select name, triangulation, m, l from %s '

    def __call__(self, **kwargs):
        return ClosedManifoldTable(self._table, **kwargs)

    def _finalize(self, M, row):
        """
        Give the closed manifold a name and do the Dehn filling.
        """
        M.set_name(row[0])
        M.dehn_fill(row[2:4])

class OrientableCuspedTable(ManifoldTable):
    """
    Iterator for all orientable cusped hyperbolic manifolds that
    can be triangulated with at most 9 ideal tetrahedra.

    >>> for M in OrientableCuspedCensus[3:6]: print(M, M.volume())
    ... 
    m007(0,0) 2.56897060
    m009(0,0) 2.66674478
    m010(0,0) 2.66674478
    >>> for M in OrientableCuspedCensus[-9:-6]: print(M, M.volume())
    ...
    o9_44241(0,0) 8.96323909
    o9_44242(0,0) 8.96736842
    o9_44243(0,0) 8.96736842
    >>> for M in OrientableCuspedCensus[4.10:4.11]: print(M, M.volume())
    ... 
    m217(0,0) 4.10795310
    m218(0,0) 4.10942659
    >>> for M in OrientableCuspedCensus(num_cusps=2)[:3]:
    ...   print(M, M.volume(), M.num_cusps())
    ... 
    m125(0,0)(0,0) 3.66386238 2
    m129(0,0)(0,0) 3.66386238 2
    m202(0,0)(0,0) 4.05976643 2
    >>> M = Manifold('m129')
    >>> M in LinkExteriors
    True
    >>> LinkExteriors.identify(M)
    5^2_1(0,0)(0,0)
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='orientable_cusped_view',
                                     **kwargs) 

class NonorientableCuspedTable(ManifoldTable):
    """
    Iterator for all orientable cusped hyperbolic manifolds that
    can be triangulated with at most 5 ideal tetrahedra.

    >>> for M in NonorientableCuspedCensus(betti=2)[:3]:
    ...   print(M, M.homology())
    ... 
    m124(0,0)(0,0)(0,0) Z/2 + Z + Z
    m128(0,0)(0,0) Z + Z
    m131(0,0) Z + Z
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='nonorientable_cusped_view',
                                     **kwargs)

class LinkTable(ManifoldTable):
    """
    Link exteriors usually know a DT code describing the assocated link.
    """
    _select = 'select name, triangulation, perm, DT from %s '

    def _manifold_factory(self, cursor, row, M=None):
        """
        Factory for "select name, triangulation" queries.
        Returns a Manifold with a DT code.
        """
        if M is None:
            M = Manifold('empty')
        buf = bytes(row[1])
        header = byte_to_int(buf[0])
        use_cobs, use_string = header&USE_COBS, header&USE_STRING
        num_cusps = header&CUSP_MASK
        knot = DTcodec(row[3])
        M._set_DTcode(knot)
        if use_string:
            M._from_string(buf[1:])
        else:
            M._from_bytes(buf[4*num_cusps + 1:])
            if use_cobs:
                cobs = decode_matrices(buf[1:4*num_cusps + 1])
                M.set_peripheral_curves('combinatorial')
                M.set_peripheral_curves(cobs)
        self._finalize(M, row)
        return M

class DTcodeTable(ManifoldTable):
    """
    Simple LinkTable which only looks up DT codes.
    Intended for use by Spherogram.
    """
    _select = 'select DT from %s '

    def _manifold_factory(self, cursor, row, M=None):
        """
        Return a DTcodec object initialized from the saved DT code.
        """
        return DTcodec(row[0])

    def _check_schema(self): 
        assert ('DT' in self.schema), 'Not a valid Link table.'

    def __getitem__(self, name):
        """
        Return the DTcode for the link with the given name.
        Raise an IndexError if there is no link with that name.
        """
        return self.find(where="name = '%s'"%name)[0]

class RolfsenTable(LinkTable):
    """
    Iterator for all knots with at most 11 crossings and links with
    at most 10 crossings, using the Rolfsen notation.  The triangulations
    were computed by Joe Christy.

    >>> for K in LinkExteriors(num_cusps=3)[-3:]:
    ...   print(K, K.volume())
    ... 
    10^3_72(0,0)(0,0)(0,0) 14.35768903
    10^3_73(0,0)(0,0)(0,0) 15.86374431
    10^3_74(0,0)(0,0)(0,0) 15.55091438
    >>> M = Manifold('8_4')
    >>> OrientableCuspedCensus.identify(M)
    s862(0,0)

    By default, the 'identify' returns the first isometric manifold it finds;
    if the optional 'extends_to_link' flag is set, it insists that meridians
    are taken to meridians.

    >>> M = Manifold('7^2_8')
    >>> LinkExteriors.identify(M)
    5^2_1(0,0)(0,0)
    >>> LinkExteriors.identify(M, extends_to_link=True)
    7^2_8(0,0)(0,0)
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='link_exteriors_view',
                                     **kwargs)

    def __call__(self, *args, **kwargs):
        if args: # backwards compatibility for LinkExteriors
            if not isinstance(args[0], int) or len(args) > 1:
                raise TypeError('Invalid specification for num_cusps.')
            if not 'num_cusps' in kwargs:
                kwargs['num_cusps'] = args[0]
        return self.__class__(**kwargs)

    def _configure(self, **kwargs):
        """
        Process the ManifoldTable filter arguments and then add
        the ones which are specific to links.
        """
        ManifoldTable._configure(self, **kwargs)
        conditions = []
        
        flavor = kwargs.get('knots_vs_links', None)
        if flavor == 'knots':
            conditions.append('cusps=1')
        elif flavor == 'links':
            conditions.append('cusps>1')
        if 'crossings' in kwargs:
            N = int(kwargs['crossings'])
            conditions.append(
                "(name like '%d^%%' or name like '%d|_%%' escape '|')"%(N,N))
        if self._filter:
            if len(conditions) > 0:
                self._filter += (' and ' + ' and '.join(conditions))
        else:
            self._filter = ' and '.join(conditions)

class RolfsenDTcodeTable(DTcodeTable):
    """
    DT codes for the Rolfsen Links.
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self, table='link_exteriors_view', **kwargs)
            
class HTLinkTable(LinkTable):
    """
    Iterator for all knots and links up to 14 crossings as tabulated
    by Jim Hoste and Morwen Thistlethwaite.  In addition to the filter
    arguments supported by all ManifoldTables, this iterator provides
    alternating=<True/False>; knots_vs_links=<'knots'/'links'>; and
    crossings=N. These allow iterations only through alternating or
    non-alternating links with 1 or more than 1 component and a
    specified crossing number.

    >>> HTLinkExteriors.identify(LinkExteriors['8_20'])
    K8n1(0,0)
    >>> Mylist = HTLinkExteriors(alternating=False,knots_vs_links='links')[8.5:8.7]
    >>> len(Mylist)
    8
    >>> for L in Mylist:
    ...   print( L.name(), L.num_cusps(), L.volume() )
    ... 
    L11n138 2 8.66421454
    L12n1097 2 8.51918360
    L14n13364 2 8.69338342
    L14n13513 2 8.58439465
    L14n15042 2 8.66421454
    L14n24425 2 8.60676092
    L14n24777 2 8.53123093
    L14n26042 2 8.64333782
    >>> for L in Mylist:
    ...   print( L.name(), L.DT_code() )
    ... 
    L11n138 [(8, -10, -12), (6, -16, -18, -22, -20, -2, -4, -14)]
    L12n1097 [(10, 12, -14, -18), (22, 2, -20, 24, -6, -8, 4, 16)]
    L14n13364 [(8, -10, 12), (6, -18, 20, -22, -26, -24, 2, -4, -28, -16, -14)]
    L14n13513 [(8, -10, 12), (6, -20, 18, -26, -24, -4, 2, -28, -16, -14, -22)]
    L14n15042 [(8, -10, 14), (12, -16, 18, -22, 24, 2, 26, 28, 6, -4, 20)]
    L14n24425 [(10, -12, 14, -16), (-18, 26, -24, 22, -20, -28, -6, 4, -2, 8)]
    L14n24777 [(10, 12, -14, -18), (2, 28, -22, 24, -6, 26, -8, 4, 16, 20)]
    L14n26042 [(10, 12, 14, -20), (8, 2, 28, -22, -24, -26, -6, -16, -18, 4)]
    """

    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='HT_links_view',
                                     db_path=alt_database_path,
                                     **kwargs)

    def _configure(self, **kwargs):
        """
        Process the ManifoldTable filter arguments and then add
        the ones which are specific to links.
        """
        ManifoldTable._configure(self, **kwargs)
        conditions = []

        alt = kwargs.get('alternating', None)
        if alt == True:
            conditions.append("name like '%a%'")
        elif alt == False:
            conditions.append("name like '%n%'")
        flavor = kwargs.get('knots_vs_links', None)
        if flavor == 'knots':
            conditions.append('cusps=1')
        elif flavor == 'links':
            conditions.append('cusps>1')
        if 'crossings' in kwargs:
            N = int(kwargs['crossings'])
            conditions.append(
                "(name like '_%da%%' or name like '_%dn%%')"%(N,N))
        if self._filter:
            if len(conditions) > 0:
                self._filter += (' and ' + ' and '.join(conditions))
        else:
            self._filter = ' and '.join(conditions)

class HTLinkDTcodeTable(DTcodeTable):
    """
    DT codes for the Hoste-Thistlethwaite Links.
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='HT_links_view',
                                     db_path=alt_database_path,
                                     **kwargs)

class CensusKnotsTable(ManifoldTable):
    """
    Iterator for all of the knot exteriors in the SnapPea Census, as
    tabulated by Callahan, Dean, Weeks, Champanerkar, Kofman and
    Patterson.  These are the knot exteriors which can be triangulated
    by at most 7 ideal tetrahedra.
    
    >>> for M in CensusKnots[3.4:3.5]:
    ...   print(M, M.volume(), LinkExteriors.identify(M))
    ... 
    K4_3(0,0) 3.47424776 False
    K5_1(0,0) 3.41791484 False
    K5_2(0,0) 3.42720525 8_1(0,0)
    K5_3(0,0) 3.48666015 9_2(0,0)
    """
    def __init__(self, **kwargs):
       return ManifoldTable.__init__(self,
                                     table='census_knots_view',
                                     **kwargs) 

class OrientableClosedTable(ClosedManifoldTable):
    """
    Iterator for 11,031 closed hyperbolic manifolds from the census by
    Hodgson and Weeks.

    >>> len(OrientableClosedCensus)
    11031
    >>> len(OrientableClosedCensus(betti=2))
    1
    >>> for M in OrientableClosedCensus(betti=2):
    ...   print(M, M.homology())
    ... 
    v1539(5,1) Z + Z
    """
    def __init__(self, **kwargs):
       return ClosedManifoldTable.__init__(self,
                                           table='orientable_closed_view',
                                           **kwargs) 

class NonorientableClosedTable(ClosedManifoldTable):
    """
    Iterator for 17 nonorientable closed hyperbolic manifolds from the
    census by Hodgson and Weeks.
    
    >>> for M in NonorientableClosedCensus[:3]: print(M, M.volume())
    ... 
    m018(1,0) 2.02988321
    m177(1,0) 2.56897060
    m153(1,0) 2.66674478
    """
    def __init__(self, **kwargs):
       return ClosedManifoldTable.__init__(self,
                                           table='nonorientable_closed_view',
                                           **kwargs) 

class IsosigManifoldTable(ManifoldTable):
    """
    Iterator for cusped manifolds in an sqlite3 table of manifolds.

    Initialize with the table name.  The table schema is required to
    include a text field called 'name' and a text field called
    'triangulation'.  The text holds the result of
    M.triangulation_isosig(), M.triangulation_isosig(decorated = True), or
    M._to_string().

    Both mapping from the manifold name, and lookup by index are
    supported.  Slicing can be done either by numerical index or by
    volume.

    The __contains__ method is supported, so M in T returns True if M
    is isometric to a manifold in the table T.  The method
    T.identify(M) will return the matching manifold from the table.
    """

    # basic select clause.  Can be overridden, e.g. to additional columns
    _select = 'select name, triangulation from %s '

    def _check_schema(self):
        assert (self.schema['name'] == 'text' and 
                self.schema['triangulation'] == 'text'
                ), 'Not a valid Manifold table.'
        
    def _manifold_factory(self, cursor, row, M=None):
        """
        Factory for "select name, triangulation" queries.
        Returns a Manifold.
        """
        
        if M is None:
            M = Manifold('empty')
            
        # Get fillings, if any
        m = split_filling_info.match(row[1])
        isosig = m.group(1)
        M._from_isosig(isosig)

        fillings = eval( '[' + m.group(2).replace(')(', '),(')+ ']', {})

        if fillings:
            M.dehn_fill(fillings)

        self._finalize(M, row)
        return M

    def _finalize(self, M, row):
        """
        Give the manifold a name.  Override this method for custom manifold
        production.
        """
        M.set_name(row[0])

class IsosigPlatonicManifoldTable(IsosigManifoldTable):
    """
    Iterator for platonic hyperbolic manifolds.
    """

    def __init__(self, table = '', db_path = platonic_database_path,
                 **filter_args):

        IsosigManifoldTable.__init__(self, table = table, db_path = db_path,
                                     **filter_args)

    def _configure(self, **kwargs):
        IsosigManifoldTable._configure(self, **kwargs)
        conditions = []

        if 'solids' in kwargs:
            N = int(kwargs['solids'])
            conditions.append('solids = %d' % N)
            
        if self._filter:
            if len(conditions) > 0:
                self._filter += (' and ' + ' and '.join(conditions))
        else:
            self._filter = ' and '.join(conditions)

class TetrahedralOrientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the tetrahedral orientable cusped hyperbolic manifolds up to
    25 tetrahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic tetrahedra.

    >>> for M in TetrahedralOrientableCuspedCensus(solids = 5):
    ...     print(M, M.volume())
    otet05_00000(0,0) 5.07470803
    otet05_00001(0,0)(0,0) 5.07470803
    >>> TetrahedralOrientableCuspedCensus.identify(Manifold("m004"))
    otet02_00001(0,0)


    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            table = 'tetrahedral_orientable_cusped_census',
            **kwargs)

class TetrahedralNonorientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the tetrahedral non-orientable cusped hyperbolic manifolds up to
    21 tetrahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic tetrahedra.

    >>> len(TetrahedralNonorientableCuspedCensus)
    25194
    >>> TetrahedralNonorientableCuspedCensus[:1.3]
    [ntet01_00000(0,0)]

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'tetrahedral_nonorientable_cusped_census',
            **kwargs)

class OctahedralOrientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the octahedral orientable cusped hyperbolic manifolds up to
    7 octahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic octahedra.

        >>> OctahedralOrientableCuspedCensus.identify(Manifold("5^2_1"))
        ooct01_00001(0,0)(0,0)

    For octahedral manifolds that are also the complement of an `Augmented
    Knotted Trivalent Graph (AugKTG) <http://arxiv.org/abs/0805.0094>`_, the
    corresponding link is included::

        >>> M = OctahedralOrientableCuspedCensus['ooct04_00034']
        >>> M.link()
        <Link: 4 comp; 17 cross>

    The link can be viewed with ``M.plink()``. To only see complements of
    AugKTGs, supply ``isAugKTG = True``::

     >>> len(OctahedralOrientableCuspedCensus(isAugKTG = True))
     238
     >>> for M in OctahedralOrientableCuspedCensus(isAugKTG = True)[:5]:
     ...     print(M, M.link().DT_code(DT_alpha=True))
     ooct02_00001(0,0)(0,0)(0,0)(0,0) DT[mdbceceJamHBlCKgdfI]
     ooct02_00002(0,0)(0,0)(0,0) DT[lcgbcIkhLBJecGaFD]
     ooct02_00003(0,0)(0,0)(0,0) DT[icebbGIAfhcEdB]
     ooct02_00005(0,0)(0,0)(0,0) DT[hcdbbFHegbDAc]
     ooct04_00027(0,0)(0,0)(0,0)(0,0) DT[zdpecbBujVtiWzsLQpxYREadhOKCmFgN]


    """

    _select = 'select name, triangulation, DT from %s '

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'octahedral_orientable_cusped_census',
            **kwargs)
        
    def _configure(self, **kwargs):
        IsosigPlatonicManifoldTable._configure(self, **kwargs)
        conditions = []
            
        if 'isAugKTG' in kwargs:
            if kwargs['isAugKTG']:
                conditions.append('isAugKTG = 1')
            else:
                conditions.append('isAugKTG = 0')
                    
        if self._filter:
            if len(conditions) > 0:
                self._filter += (' and ' + ' and '.join(conditions))
        else:
            self._filter = ' and '.join(conditions)

    def _finalize(self, M, row):
        IsosigPlatonicManifoldTable._finalize(self, M, row)
        if row[2]:
            M._set_DTcode(DTcodec(row[2]))

class OctahedralNonorientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the octahedral non-orientable cusped hyperbolic manifolds up to
    5 octahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic octahedra.

    >>> for M in OctahedralNonorientableCuspedCensus(solids = 3, betti = 3,cusps = 4):
    ...     print(M, M.homology())
    noct03_00007(0,0)(0,0)(0,0)(0,0) Z/2 + Z + Z + Z
    noct03_00029(0,0)(0,0)(0,0)(0,0) Z/2 + Z + Z + Z
    noct03_00047(0,0)(0,0)(0,0)(0,0) Z/2 + Z + Z + Z
    noct03_00048(0,0)(0,0)(0,0)(0,0) Z/2 + Z + Z + Z

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'octahedral_nonorientable_cusped_census',
            **kwargs)

class CubicalOrientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the cubical orientable cusped hyperbolic manifolds up to
    7 cubes, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic octahedra.

    >>> M = TetrahedralOrientableCuspedCensus['otet05_00001']
    >>> CubicalOrientableCuspedCensus.identify(M)
    ocube01_00002(0,0)(0,0)

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'cubical_orientable_cusped_census',
            **kwargs)

class CubicalNonorientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the cubical non-orientable cusped hyperbolic manifolds up to
    5 cubes, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic octahedra.

    >>> for M in CubicalNonorientableCuspedCensus[-3:]:
    ...     print(M, M.volume())
    ncube05_30945(0,0) 25.37354016
    ncube05_30946(0,0)(0,0) 25.37354016
    ncube05_30947(0,0)(0,0) 25.37354016

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'cubical_nonorientable_cusped_census',
            **kwargs)

class DodecahedralOrientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the dodecahedral orientable cusped hyperbolic manifolds up to
    2 dodecahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic dodecahedra.

    Complement of one of the dodecahedral knots by Aitchison and Rubinstein::

      >>> M=DodecahedralOrientableCuspedCensus['odode02_00913']
      >>> M.dehn_fill((1,0))
      >>> M.fundamental_group()
      Generators:
      <BLANKLINE>
      Relators:
      <BLANKLINE>

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'dodecahedral_orientable_cusped_census',
            **kwargs)

class DodecahedralNonorientableCuspedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the dodecahedral non-orientable cusped hyperbolic manifolds up to
    2 dodecahedra, i.e., manifolds that admit a tessellation by regular ideal
    hyperbolic dodecahedra.

    >>> len(DodecahedralNonorientableCuspedCensus)
    4146

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'dodecahedral_nonorientable_cusped_census',
            **kwargs)

class IcosahedralNonorientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the icosahedral non-orientable closed hyperbolic manifolds up
    to 3 icosahedra, i.e., manifolds that admit a tessellation by regular finite
    hyperbolic icosahedra.

    >>> list(IcosahedralNonorientableClosedCensus)
    [nicocld02_00000(1,0)]

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'icosahedral_nonorientable_closed_census',
            **kwargs)

class IcosahedralOrientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the icosahedral orientable closed hyperbolic manifolds up
    to 4 icosahedra, i.e., manifolds that admit a tessellation by regula finite
    hyperbolic icosahedra.
    
    >>> IcosahedralOrientableClosedCensus[0].volume()
    4.68603427

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'icosahedral_orientable_closed_census',
            **kwargs)

class CubicalNonorientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the cubical non-orientable closed hyperbolic manifolds up
    to 10 cubes, i.e., manifolds that admit a tessellation by regular finite
    hyperbolic cubes.

    >>> len(CubicalNonorientableClosedCensus)
    93

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'cubical_nonorientable_closed_census',
            **kwargs)

class CubicalOrientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the cubical orientable closed hyperbolic manifolds up
    to 10 cubes, i.e., manifolds that admit a tessellation by regular finite
    hyperbolic cubes.

    >>> len(CubicalOrientableClosedCensus)
    69

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'cubical_orientable_closed_census',
            **kwargs)

class DodecahedralNonorientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the dodecahedral non-orientable closed hyperbolic manifolds up
    to 2 dodecahedra, i.e., manifolds that admit a tessellation by regular finite
    hyperbolic dodecahedra with a dihedral angle of 72 degrees.

    >>> DodecahedralNonorientableClosedCensus[0].volume()
    22.39812948

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'dodecahedral_nonorientable_closed_census',
            **kwargs)

class DodecahedralOrientableClosedTable(IsosigPlatonicManifoldTable):
    """
    Iterator for the dodecahedral orientable closed hyperbolic manifolds up
    to 3 dodecahedra, i.e., manifolds that admit a tessellation by regular finite
    hyperbolic dodecahedra with a dihedral angle of 72 degrees.

    The Seifert-Weber space::

      >>> M = DodecahedralOrientableClosedCensus(solids = 1)[-1]
      >>> M.homology()
      Z/5 + Z/5 + Z/5

    """

    def __init__(self, **kwargs):
        return IsosigPlatonicManifoldTable.__init__(
            self,
            'dodecahedral_orientable_closed_census',
            **kwargs)

# Instantiate our tables ...
OrientableCuspedCensus = OrientableCuspedTable()
NonorientableCuspedCensus = NonorientableCuspedTable()
OrientableClosedCensus = OrientableClosedTable()
NonorientableClosedCensus = NonorientableClosedTable()
LinkExteriors = RolfsenTable()
CensusKnots = CensusKnotsTable()
HTLinkExteriors = HTLinkTable()
RolfsenDTcodes = RolfsenDTcodeTable()
HTLinkDTcodes = HTLinkDTcodeTable()

TetrahedralOrientableCuspedCensus = TetrahedralOrientableCuspedTable()
TetrahedralNonorientableCuspedCensus = TetrahedralNonorientableCuspedTable()
OctahedralOrientableCuspedCensus = OctahedralOrientableCuspedTable()
OctahedralNonorientableCuspedCensus = OctahedralNonorientableCuspedTable()
CubicalOrientableCuspedCensus = CubicalOrientableCuspedTable()
CubicalNonorientableCuspedCensus = CubicalNonorientableCuspedTable()
DodecahedralOrientableCuspedCensus = DodecahedralOrientableCuspedTable()
DodecahedralNonorientableCuspedCensus = DodecahedralNonorientableCuspedTable()
IcosahedralNonorientableClosedCensus = IcosahedralNonorientableClosedTable()
IcosahedralOrientableClosedCensus = IcosahedralOrientableClosedTable()
CubicalNonorientableClosedCensus = CubicalNonorientableClosedTable()
CubicalOrientableClosedCensus = CubicalOrientableClosedTable()
DodecahedralNonorientableClosedCensus = DodecahedralNonorientableClosedTable()
DodecahedralOrientableClosedCensus = DodecahedralOrientableClosedTable()

# Identify a Manifold

__all_tables__ = (
    OrientableCuspedCensus,
    NonorientableCuspedCensus,
    OrientableClosedCensus,
    NonorientableClosedCensus,
    LinkExteriors,
    CensusKnots,
    HTLinkExteriors,
    TetrahedralOrientableCuspedCensus,
    TetrahedralNonorientableCuspedCensus,
    OctahedralOrientableCuspedCensus,
    OctahedralNonorientableCuspedCensus,
    CubicalOrientableCuspedCensus,
    CubicalNonorientableCuspedCensus,
    DodecahedralOrientableCuspedCensus,
    DodecahedralNonorientableCuspedCensus,
    IcosahedralNonorientableClosedCensus,
    IcosahedralOrientableClosedCensus,
    CubicalNonorientableClosedCensus,
    CubicalOrientableClosedCensus,
    DodecahedralNonorientableClosedCensus,
    DodecahedralOrientableClosedCensus
)

def identify(manifold):
    return dict(( T, (T.identify(manifold), T.identify(manifold, True)) )
                 for T in __all_tables__)

# Lookup a DT code

## Regex objects for recognizing link names (copied from the SnapPy extension) 
is_knot_complement = re.compile('([0-9]+_[0-9]+)$')
is_link_complement1_pat = '(?P<crossings>[0-9]+)[\^](?P<components>[0-9]+)[_](?P<index>[0-9]+)$'
is_link_complement2_pat = '(?P<crossings>[0-9]+)[_](?P<index>[0-9]+)[\^](?P<components>[0-9]+)$'
is_link_complement3_pat = '[lL](?P<components>[0-9]{1})(?P<crossings>[0-9]{2})(?P<index>[0-9]+)$'
is_link_complement1 = re.compile(is_link_complement1_pat)
is_link_complement2 = re.compile(is_link_complement2_pat)
is_link_complement3 = re.compile(is_link_complement3_pat)
rolfsen_link_regexs = [is_link_complement1, is_link_complement2, is_link_complement3]
is_HT_knot = re.compile('(?P<crossings>[0-9]+)(?P<alternation>[an])(?P<index>[0-9]+)$')
is_HT_link = re.compile('[KL][0-9]+[an]([0-9]+)$')

def lookup_DT(name):
        if is_HT_link.match(name):
            return HTLinkDTcodes[name]
        if is_knot_complement.match(name):
            return RolfsenDTcodes[name]
        m = is_HT_knot.match(name)
        if m:
            HT_name = 'K%d%s%d'%(int(m.group('crossings')),
                                 m.group('alternation'),
                                 int(m.group('index')))
            return HTLinkDTcodes[HT_name]
        for regex in rolfsen_link_regexs:
            m = regex.match(name)
            if m:
                if int(m.group('components')) > 1:
                    rolfsen_name = '%d^%d_%d' % (int(m.group('crossings')),
                                                 int(m.group('components')),
                                                 int(m.group('index')))
                else:
                    rolfsen_name = '%d_%d' % (int(m.group('crossings')),
                                              int(m.group('index')))
                return RolfsenDTcodes[rolfsen_name]
        raise IndexError('Unrecognized format for a link name.')


# Test routines.
def test_census_database():
    L = OrientableCuspedDB
    for M in CensusKnots():
        print(M, L.identify(M))

def test():
    import re
    pairs = [ (CensusKnots(), CensusKnotsDB),
              (OrientableCuspedCensus(), OrientableCuspedDB)]
    for census, db in pairs:
        for M in census:
            if M.name()[0] == 't':
                break
            N = db.identify(M)
            assert repr(M) == repr(N)
            G, H = M.fundamental_group(), N.fundamental_group()
            if (G.relators() != H.relators() or
                G.peripheral_curves() != H.peripheral_curves()):
                print(M)

if __name__ == '__main__':
    test()
