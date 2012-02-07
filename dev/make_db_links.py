from snappy import *
import os, sqlite3, binascii, bz2
from hashlib import md5
from census import *
from lookup import SmallHTWKnots

snappy_schema = """
CREATE TABLE census (
 id integer primary key,
 name text,
 volume real,
 cusps int,
 hash blob, 
 triangulation blob)
"""

insert_query = """insert into census
(name, volume, cusps, hash, triangulation)
values ('%s', %s, %s, X'%s', ?)"""

def create_census(connection):
    connection.execute(snappy_schema)
    connection.commit()

def insert_manifold(connection, mfld):
    name = mfld.name()
    volume = mfld.volume()
    if volume < 1.0:
        volume = 0.0
    tri = mfld.without_hyperbolic_structure()
    triangulation = sqlite3.Binary(bytes(bz2.compress(tri._to_string())))
    hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
    cusps = mfld.num_cusps()
    query = insert_query%(name, volume, cusps, hash)
    connection.execute(query, (triangulation,))
    
def make_links():
    if os.path.exists('links.sqlite'):
        os.remove('links.sqlite')
    snappy_connection = sqlite3.connect('links.sqlite')
    create_census(snappy_connection)
    for n in range(1, 6):
        for M in LinkExteriors(n):
            insert_manifold(snappy_connection, M)
    snappy_connection.commit()

def make_new_links():
    if os.path.exists('new_knots.sqlite'):
        os.remove('new_knots.sqlite')
    snappy_connection = sqlite3.connect('new_knots.sqlite')
    create_census(snappy_connection)
    for M in SmallHTWKnots():
        insert_manifold(snappy_connection, M)
    snappy_connection.commit()

if __name__ == '__main__':
    make_links()
    make_new_links()
