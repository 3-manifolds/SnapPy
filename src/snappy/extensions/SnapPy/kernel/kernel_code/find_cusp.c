/*
 *  find_cusp.c
 *
 *  This function provides the following utility for use within the kernel.
 *
 *      Vertex *find_cusp(Triangulation *manifold, int cusp_index);
 *
 *  The UI refers to cusps by their indices, because it doesn't know about
 *  the Vertex data structure.  When a kernel function receives a cusp_index,
 *  it may call the function find_cusp() to convert the index to an actual
 *  pointer to the corresponding Vertex data structure.  If find_cusp() cannot
 *  find a Vertex with the given cusp_index, it calls uFatalError().
 */

#include "kernel.h"
SNAPPEA_NAMESPACE_BEGIN_SCOPE

Vertex *find_cusp(
    Triangulation   *manifold,
    int             cusp_index)
{
    Vertex  *cusp;

    for (cusp = manifold->vertex_list_begin.next;
         cusp != &manifold->vertex_list_end;
         cusp = cusp->next)

        if (cusp->cusp_index == cusp_index)

            return cusp;

    /*
     *  No Vertex with the given index was found.
     */
    uFatalError("find_cusp", "find_cusp");

    /*
     *  The C++ compiler would like a return value, even though
     *  we never return from the uFatalError() call.
     */
    return NULL;
}
SNAPPEA_NAMESPACE_END_SCOPE
