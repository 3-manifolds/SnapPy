import snappy
if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')

snappy.number.Number._accuracy_for_testing = 6

values = set()
while True:
    vol = repr(snappy.ManifoldHP('5_2').complex_volume())
    if vol in values:
        continue
    values.add(vol)
    print vol
    

