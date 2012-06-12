import snappy
from . import shapes, polished_reps, find_field

def tetrahedra_shapes(manifold):
    def func(prec):
        return shapes.polished_tetrahedra_shapes(manifold, prec)
    return find_field.SetOfApproximateAlgebraicNumbers(func)
