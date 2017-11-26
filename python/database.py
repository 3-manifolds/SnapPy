"""
This module uses the Manifold class from the SnapPy extension
module. After the Manifold class has been defined, the "__init__.py"
module for "snappy" itself sets::

  database.Manifold = Manifold

This file simply defines the ManifoldTable class.  The ManifoldTables
that come with actual data are defined in external packages, such as
the required "snappy_manifolds".  Such a package must provide a
"get_tables" function that accepts the below "ManifoldTable" as input
and returns a list of subclasses of "ManifoldTable".
"""

from __future__ import print_function
from .db_utilities import decode_torsion, decode_matrices, db_hash
from .sage_helper import _within_sage
from spherogram.codecs import DTcodec
import sys, sqlite3, re, os, random, importlib, collections

if _within_sage:
    import sage.all
    def is_int(slice):
        return isinstance(slice, (sage.all.Integer,int))
        
    def is_int_or_none(slice):
        return isinstance(slice, (sage.all.Integer,int, type(None)))

    def is_float_or_none(slice):
        return isinstance(slice, (float, sage.all.RealDoubleElement,
                                  sage.rings.real_mpfr.RealNumber, type(None)))
else:
    def is_int(slice):
        return isinstance(slice, int)
    
    def is_int_or_none(slice):
        return isinstance(slice, (int, type(None)))

    def is_float_or_none(slice):
        return isinstance(slice, (float, type(None)))

split_filling_info = re.compile('(.*?)((?:\([0-9 .+-]+,[0-9 .+-]+\))*$)')

def connect_to_db(db_path):
    """
    Open the given sqlite database, ideally in read-only mode.
    """
    if sys.version_info >= (3,4):
        uri = 'file:' + db_path + '?mode=ro'
        return sqlite3.connect(uri, uri=True)
    elif sys.platform.startswith('win'):
        try:
            import apsw
            return apsw.Connection(db_path, flags=apsw.SQLITE_OPEN_READONLY)
        except ImportError:
            return sqlite3.connect(db_path)
    else:
        return sqlite3.connect(db_path)

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
    # basic select clause.  Can be overridden, e.g. to add additional columns
    _select = 'select name, triangulation from %s '

    def __init__(self, table='', db_path=None,
                 mfld_hash=mfld_hash, **filter_args):
        self._table = table
        self.mfld_hash = mfld_hash
        self._connection = connect_to_db(db_path)
        self._cursor = self._connection.cursor()
        self._set_schema()
        self._check_schema()
        self._configure(**filter_args)
        self._get_length()
        self._get_max_volume()
        self._select = self._select%table

    def _set_schema(self):
        cursor, table = self._cursor, self._table
        rows = cursor.execute("pragma table_info('%s')" % table).fetchall()
        self.schema = dict([(row[1],row[2].lower()) for row in rows])

    def _check_schema(self):
        assert (self.schema['name'] == 'text' and 
                self.schema['triangulation'] == 'text'
                ), 'Not a valid Manifold table.'

    @property
    def filter(self):
        return self._filter

    def _get_length(self):
        where_clause = 'where ' + self._filter if self._filter else '' 
        length_query = 'select count(*) from %s %s' % (self._table,
                                                       where_clause)
        cursor = self._cursor.execute(length_query)
        self._length = cursor.fetchone()[0]

    def _get_max_volume(self):
        where_clause = 'where ' + self._filter if self._filter else '' 
        vol_query = 'select max(volume) from %s %s' % (self._table,
                                                       where_clause)
        cursor = self._cursor.execute(vol_query)
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
            query += ' where %s order by id'%self._filter
        for row in self._cursor.execute(query):
            yield self._manifold_factory(row)

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
                filter = ' and '.join(conditions)
                return self.__class__(filter=filter)
            elif (is_int_or_none(start) and is_int_or_none(stop)):
                if start is None:
                    start = 0
                elif start < 0:
                    start = int(self._length + start)
                if stop is None:
                    stop = self._length
                elif stop < 0:
                    stop = int(self._length + stop)
                conditions = []
                base_query = 'select id from %s ' % self._table
                if self._filter:
                    base_query += 'where %s ' % self._filter
                query = base_query + 'limit 1 offset %d' % start
                start_id = self._cursor.execute(query).fetchone()
                if start_id is not None:
                    conditions.append('id >= %d' % start_id[0])
                query = base_query + 'limit 1 offset %d' % stop
                stop_id = self._cursor.execute(query).fetchone()
                if stop_id is not None:
                    conditions.append('id < %d' % stop_id[0])
                if self._filter:
                    conditions.append(self._filter)
                return self.__class__(filter=' and '.join(conditions))
            else:
                raise IndexError(
                    'Use two ints or two floats for start and stop.')
        elif is_int(index):
            if index < 0:
                index = self._length + index
            matches = self.find(limit=1, offset=index)
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
    
    def _manifold_factory(self, row, M=None):
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
        Give the manifold a name.  Override this method for custom
        manifold production.
        """
        M.set_name(row[0])

    def _one_manifold(self, name, M):
        """
        Inflates the given empty Manifold with the table manifold
        with the specified name.
        """
        if hasattr(self, '_regex'):
            if self._regex.match(name) is None:
                raise KeyError('The manifold %s was not found.'%name)
        cursor = self._cursor.execute(self._select + "where name='" + name + "'")
        rows = cursor.fetchall()
        if len(rows) != 1:
            raise KeyError('The manifold %s was not found.'%name)
        return self._manifold_factory(rows[0], M)
                
    def keys(self):
        """
        Return the list of column names for this manifold table.
        """
        return self.schema.keys()
    
    def find(self, where=None, order_by='id', limit=None, offset=None):
        """
        Return a list of up to limit manifolds stored in this table,
        satisfying the where clause, and ordered by the order_by
        clause.  If limit is None, all matching manifolds are
        returned.  If the offset parameter is set, the first offset
        matches are skipped.
        """
        conditions = [cond for cond in [self._filter, where] if cond]
        suffix = ' where ' if conditions else ' '
        suffix += ' and '.join(conditions)
        suffix += ' order by %s' % order_by
        if limit is not None:
            suffix += ' limit %d' % limit
        if offset is not None:
            suffix += ' offset %d' % offset
        cursor = self._cursor.execute(self._select + suffix)
        return [self._manifold_factory(row) for row in cursor.fetchall()]

    def siblings(self, mfld):
        """
        Return all manifolds in the census which have the same hash value.
        """
        vol = mfld.volume()
        epsilon = vol/1e5
        v_lower, v_upper = vol - epsilon, vol + epsilon
        cusps = mfld.cusp_info('is_complete').count(True)
        H = mfld.homology()
        betti = H.betti_number()
        torsion = [c for c in H.elementary_divisors() if c!=0]
        initial_candidates = self.find(
          "volume between %f and %f and cusps=%d and betti=%d and torsion='%s'"
          % (v_lower, v_upper, cusps, betti, torsion))
        if len(initial_candidates) == 0:
            return []
        return self.find("hash = '%s'"%self.mfld_hash(mfld))

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
            if mfld.volume() > self._max_volume + 0.1:
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


# The below function is used to add ManifoldTables defined in external
# packages.

this_module = sys.modules[__name__]
__all_tables__ = []

def add_tables_from_package(package_name, must_succeed=True):
    """
    Given a string with the name of an importable Python package that
    implements a "get_tables" function, load all the tables provided
    and put the results where the rest of SnapPy can find them.
    Returns an ordered dictionary of pairs (table_name, table).
    """
    
    try:
        package = importlib.import_module(package_name)
    except ImportError:
        if not must_succeed:
            return dict()
        else:
            raise ImportError('ManifoldTable package %s not found'
                              % package_name)

    tables = collections.OrderedDict()

    for table in package.get_tables(ManifoldTable):
        tables[table.__class__.__name__] = table
        __all_tables__.append(table)

    # Add to namespace of this module:
    for name, table in tables.items():
        setattr(this_module, name, table)
    
    
    # We also store the tables here so that their doctests can be
    # checked.
    if not hasattr(this_module, '__test__'):
        this_module.__test__ = dict()
    for name, table in tables.items():
        this_module.__test__[name] = table.__class__

    return tables
