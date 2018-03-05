"""
Run this in Sage and SnapPy and see that the last line of the result differs dramatically.

If we multiply two numbers with different precisions, we expect the result to be
of the lower precision. This is consistent between Sage and SnapPy.
However, if we cast the low precision to a type of higher precision and then multiply, we
expect the result to be of the type of the higher precision.
This works in Sage, but it doesn't work for SnapPy Numbers.

This affects the code computing cusp translations, look for "_SnapPyNumberHack".
"""

from snappy import Manifold
M = Manifold("m004")

# Let's take a real 53bit precision number
v = M.cusp_neighborhood().volume(0)
# And a complex 1000bit precision number
z = M.tetrahedra_shapes('rect', bits_prec = 1000)[0]
# High precision real field
RF = z.real().parent()

print("Expect low precision number: ", v)
print("Expect high precision number:", z)
print("Product, expect low precision number:", v * z)
print("Casting to high precision:", RF(v))
print("Product again, this time it should be high precision (!!!):", RF(v) * z)

