import snappy
if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')

snappy.number.Number._accuracy_for_testing = 6

# addl_code/complex_volume.c uses the dilog callback into pari
#
# This callback used to give the wrong result when the real part
# of the dilog was small and pari returned as string something like
# "5.234 E-3". pari puts a space before the "E" causing the
# quad double library to parse the result incorrectly.
# See gen2Complex in SnapPycore.pxi

c = 0

values = set()
while True:
    c += 1
    vol = repr(snappy.ManifoldHP('5_2').complex_volume())
    if vol in values:
        continue
    values.add(vol)
    print vol, "after", c, "iterations"
    

