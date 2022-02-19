from ..snap.t3mlite import Tetrahedron, Perm4, Mcomplex, simplex

from .barycentric_subdivider import perms, transpositions, perm_to_index

def crush_link(original_triangulation,
               subdivided_triangulation):

    num_tets = len(original_triangulation.Tetrahedra)

    tetrahedra = [
        Tetrahedron() if tet.Class[simplex.E01].link_component_index == 0 else None
        for tet in subdivided_triangulation ]

    for original_tet in original_triangulation:
        for i, perm in enumerate(perms):
            tet_index = 24 * original_tet.Index + i
            tet0 = tetrahedra[tet_index]
            if tet0:
                for face in range(3):
                    other_perm = perm * transpositions[face]
                    j = perm_to_index(other_perm)
                    other_tet_index = 24 * original_tet.Index + j
                    tet1 = tetrahedra[other_tet_index]
                    if face == 1 and not tet1:
                        other_perm = perm * Perm4((2,1,0,3))
                        j = perm_to_index(other_perm)
                        other_tet_index = 24 * original_tet.Index + j
                        tet1 = tetrahedra[other_tet_index]
                    if j > i:
                        tet0.attach(simplex.TwoSubsimplices[face], tet1, (0,1,2,3))

                face = perm.image(simplex.F3)
                other_tet = original_tet.Neighbor[face]
                other_perm = original_tet.Gluing[face] * perm
                j = perm_to_index(other_perm)
                other_tet_index = 24 * other_tet.Index + j
                if other_tet_index > tet_index:
                    tet1 = tetrahedra[other_tet_index]
                    tet0.attach(simplex.F3, tet1, (0,1,2,3))

                tet0.post_drill_infos = (
                    subdivided_triangulation.Tetrahedra[tet_index].post_drill_infos)
                tet0.PeripheralCurves = (
                    subdivided_triangulation.Tetrahedra[tet_index].PeripheralCurves)
                tet0.base_for_peripherals = (
                    subdivided_triangulation.Tetrahedra[tet_index].base_for_peripherals)

    return Mcomplex([ tet for tet in tetrahedra if tet])
