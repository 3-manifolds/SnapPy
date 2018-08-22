import snappy
from sage.all import RealField, ComplexField, I

M = snappy.ManifoldHP('m004')
shapes = M.tetrahedra_shapes('rect')
M.set_tetrahedra_shapes(filled_shapes=shapes, fillings=[(1, 0)])
RR = RealField(212)
CC = ComplexField(25)
pi = RR.pi()
M.set_target_holonomy(RR(0))
for i in range(201):
    M.set_target_holonomy(i*pi/50*I)
    if i % 25 == 0:
        shapes = M.tetrahedra_shapes()
        z = shapes[1]
        print(i, CC(z['rect']), CC(z['log']))

print(M.solution_type())
print([CC(sh) for sh in M.tetrahedra_shapes('rect')])

    
