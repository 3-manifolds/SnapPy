from snappy import Triangulation
from snappy.ptolemy.homology import *
import sys

def to_index(s):
    var, face_num, tet_num = s.split('_')
    if var != 's':
        raise Exception("Not s")
    return 4 * int(tet_num) + int(face_num)

trig = Triangulation(sys.argv[1], remove_finite_vertices = False)

# trig.reverse_orientation()

face_classes = trig._ptolemy_equations_identified_face_classes()

H = homology_basis_representatives_with_orders(
    trig._ptolemy_equations_boundary_map_2()[0],
    trig._ptolemy_equations_boundary_map_3()[0],
    0)

for h, o in H:
    
    weights = (4 * trig.num_tetrahedra()) * [ 0 ]

    for e, face_class in zip(h, face_classes):
        for i, face in enumerate(face_class[2:]):
            weights[to_index(face)] = (-1)**i * e

    print("Order: ", o)
    print(weights)

