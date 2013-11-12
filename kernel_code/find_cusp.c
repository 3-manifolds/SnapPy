/*
 *  find_cusp.c
 *
 *  This function provides the following utility for use within the kernel.
 *
 *      Cusp *find_cusp(Triangulation *manifold, int cusp_index);
 *
 *  The UI refers to cusps by their indices, because it doesn't know about
 *  the Cusp data structure.  When a kernel function receives a cusp_index,
 *  it may call the function find_cusp() to convert the index to an actual
 *  pointer to the corresponding Cusp data structure.  If find_cusp() cannot
 *  find a Cusp with the given cusp_index, it calls uFatalError().
 */

#include "kernel.h"
#include "kernel_namespace.h"

Cusp *find_cusp(
    Triangulation   *manifold,
    int             cusp_index)
{
    Cusp    *cusp;

    for (cusp = manifold->cusp_list_begin.next;
         cusp != &manifold->cusp_list_end;
         cusp = cusp->next)

        if (cusp->index == cusp_index)

            return cusp;

    /*
     *  No Cusp with the given index was found.
     */
    uFatalError("find_cusp", "find_cusp");

    /*
     *  The C++ compiler would like a return value, even though
     *  we never return from the uFatalError() call.
     */
    return NULL;
}
#include "end_namespace.h"
