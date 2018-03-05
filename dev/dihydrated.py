import snappy, string

def str_to_int(x):
    return string.ascii_lowercase.index(x)

def rehydrate(dehydration):
    dehydration = map(str_to_int, dehydration)
    num_tet = dehydration[0]
    i = 0
    which_tets = []
    while 4*i < num_tet + 3:
        which_tets.append(16*dehydration[2*i + 1]  + dehydration[2*i + 2])
        i += 1
    which_gluings = dehydration[-num_tet:]
    return bytearray([num_tet] + which_tets + which_gluings)
    
                         

b = rehydrate('cabbbbaeil')
#M = snappy.Manifold('empty')
#M._from_bytes(str(b))
