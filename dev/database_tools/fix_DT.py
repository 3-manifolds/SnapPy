"""
2016/10/7: Stavros pointed out that there are inconsistencies in
the DT codes storred with LinkExteriors.  Specifically, sometimes the
DT code gives the mirror image of the link associated with the
triangulation.
"""

import sqlite3, re, os, sys, gzip
import pandas as pd


dt_codes = dict(re.findall( '(\S+)\s+(\S+)$',
                                gzip.open('ChristyDT.gz').read(),
                                re.MULTILINE))


con = sqlite3.connect('./manifolds.sqlite')

df = pd.read_sql_query('select name, DT, id from link_exteriors', con)
assert set(df.name) == set(dt_codes.keys())
df['fixed_DT'] = df.name.apply(lambda x:dt_codes[x])

if sum(df.DT != df.fixed_DT) == 0:
    print('Database already fixed!')
else:
    for i, row in df.iterrows():
        id, DT = row['id'], row['fixed_DT']
        con.execute('update link_exteriors set DT = "%s" where id=%d'
                    % (DT, id))
    con.commit()
    print('Fixed database!')


    
    
