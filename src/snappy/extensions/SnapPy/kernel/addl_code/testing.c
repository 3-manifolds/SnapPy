#include "testing.h"

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

void testing_compute_cusp_orientabilities(
    Triangulation * manifold)
{
    for (Cusp * cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)
        cusp->orientability = unknown_cusp_orientability;
    compute_cusp_orientabilities(manifold);
}

SNAPPEA_NAMESPACE_END_SCOPE
