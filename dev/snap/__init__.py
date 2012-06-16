import snappy
from .shapes import polished_tetrahedra_shapes
from .polished_reps import polished_holonomy
from .find_field import ListOfApproximateAlgebraicNumbers

#def tetrahedra_shapes(self, part=None, fixed_alignment=True, prec=None):
#    if part=='rect' and fixed_alignment==True and prec:
#        return snap.polished_tetrahedra_shapes(self, prec)
#    else:
#        return Manifold.tetrahedra_shapes(part, fixed_alignment, prec)

def tetrahedra_field_gens(manifold):
    def func(prec):
        return polished_tetrahedra_shapes(manifold, prec)
    return ListOfApproximateAlgebraicNumbers(func)

def trace_field_gens(manifold, fundamental_group_args = []):
    def func(prec):
        return polished_holonomy(manifold, prec, fundamental_group_args).trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

def invariant_trace_field_gens(manifold, fundamental_group_args = []):
    def func(prec):
        return polished_holonomy(manifold, prec, fundamental_group_args).invariant_trace_field_generators()
    return ListOfApproximateAlgebraicNumbers(func)

def holonomy_matrix_entries(manifold, fundamental_group_args = []):
    def func(prec):
        G = polished_holonomy(manifold, prec, fundamental_group_args)
        return sum( [G.SL2C(g).list() for g in G.generators()], [])
    return ListOfApproximateAlgebraicNumbers(func)
        
    
