#include "kernel.h"

struct extra
{
        /*
         *      The first four Tetrahedra lopped off in the above
         *      algorithm will be called "outer vertex Tetrahedra".
         *      outer_vertex_tet[i] is a pointer to the outer vertex
         *      Tetrahedron at vertex i of the old Tetrahedron.
         */

        Tetrahedron     *outer_vertex_tet[4];

        /*
         *      The "inner vertex Tetrahedra" are the remaining
         *      Tetrahedra which are naturally associated to the
         *      ideal vertices of the old Tetrahedron.
         *      inner_vertex_tet[i] is a pointer to the new Tetrahedron
         *      which meets outer_vertex_tet[i] along it's i-th face.
         */

        Tetrahedron *inner_vertex_tet[4];

        /*
         *      The "edge Tetrahedra" are the 12 Tetrahedra which have
         *      precisely one edge contained within an edge of the
         *      old Tetrahedron.  edge_tet[i][j] is the Tetrahedron
         *      which has a face contained in face i of the old
         *      Tetrahedron, on the side (of face i) opposite vertex j.
         *      edge_tet[i][j] is defined iff i != j.
         */

        Tetrahedron *edge_tet[4][4];

        /*
         *      The "face Tetrahedra" are the 12 remaining Tetrahedra.
         *      face_tet[i][j] has a face contained in face i of the
         *      old Tetrahedron, on the side near vertex j (of the old
         *      Tetrahedron).  face_tet[i][j] is defined iff i != j.
         */

        Tetrahedron *face_tet[4][4];

};
