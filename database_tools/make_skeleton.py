import sqlite3

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

if __name__ == '__main__':
    connection = sqlite3.connect('manifolds.sqlite')
    create_manifold_tables(connection)
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
