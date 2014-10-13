import snappy
if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')

snappy.number.Number._accuracy_for_testing = 6

# complex_volume() of a cusped manifold ultimately
# calls into addl_code/complex_volume.c
# When using quad-double precision instead of normal precision, the
# method every once in a while returns a volume that is off by
# a multiple of pi^2/12.
#
# Possible explanation: the code in addl_code/complex_volume.c
# has a deterministic part and a non-deterministic part.
# The deterministic part uses the triangulation as is provided,
# but because the triangulation is not ordered, it might be off by
# a multiple of pi^2/6. The non-deterministic part is doing 1-4 moves
# 3-2 moves and randomly shoots the newly vertices to infinity.
# Because the deterministic part has higher precision, we use the
# result of the non-deterministic part to add a multiple of pi^2/12
# (should be pi^2/6) to the deterministic result.
#
# There are checks that the shooting the vertices to infinity are
# not producing degenerate tetrahedra. Apparently, these checks fail
# when using high-precision types.

c = 0

values = set()
while True:
    c += 1
    vol = repr(snappy.ManifoldHP('5_2').complex_volume())
    if vol in values:
        continue
    values.add(vol)
    print vol, "after", c, "iterations"
    

