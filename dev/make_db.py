from snappy import *
import os
import sqlite3
import binascii
from hashlib import md5
from census import *

snappy_schema = """
CREATE TABLE census (
 id integer primary key,
 name text,
 volume real,
 chernsimons real,
 hash blob,
 triangulation blob)
"""

insert_query = """insert into census
(name, volume, chernsimons, hash, triangulation)
values ('%s', %s, %s, X'%s', X'%s')"""

def create_census(connection):
    connection.execute(snappy_schema)
    connection.commit()
    
def insert_manifold(connection, mfld):
    name = mfld.name()
    volume = mfld.volume()
    cs = mfld.chern_simons()
    triangulation = binascii.hexlify(mfld._to_bytes())
    hash = md5(standard_hashes.combined_hash(mfld)).hexdigest()
    query = insert_query%(name, volume, cs, hash, triangulation)
    try:
        connection.execute(query)
    except:
        print query

def main():
    if os.path.exists('census.sqlite'):
        print '%s already exists!'%'census.sqlite'
    snappy_connection = sqlite3.connect('census.sqlite')
    create_census(snappy_connection)
    for M in OrientableCuspedCensus():
        insert_manifold(snappy_connection, M)
    snappy_connection.commit()
if __name__ == '__main__':
    main()
