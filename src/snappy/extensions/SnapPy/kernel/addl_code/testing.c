#include "testing.h"

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

void testing_compute_vertex_link_orientabilities(
    Triangulation * manifold)
{
    for (Vertex * cusp = manifold->vertex_list_begin.next;
         cusp != &manifold->vertex_list_end;
         cusp = cusp->next)
        cusp->orientability = unknown_orientability;
    compute_vertex_link_orientabilities(manifold);
}

SNAPPEA_NAMESPACE_END_SCOPE
