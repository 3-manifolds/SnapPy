import snappy
M1 = snappy.Manifold('m1.tri')
M2 = snappy.Manifold('m2.tri')
print M1.volume(), M2.volume()
print M2.is_isometric_to(M1)
