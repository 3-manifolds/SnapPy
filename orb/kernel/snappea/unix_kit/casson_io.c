/*
 *  cassion_io.c
 *
 */

#include "casson_io.h"
#include "ostream.h"

#include "kernel.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>

/* This should be in a header from kernel.h */
Boolean contains_flat_tetrahedra( Triangulation *manifold );

typedef struct CassonFormat CassonFormat;
typedef struct EdgeInfo EdgeInfo;
typedef struct TetEdgeInfo TetEdgeInfo;

static SolutionType string_to_solution_type(char *str);
static Boolean skip_blanks(char **str);
static Boolean fill_casson_struct(CassonFormat *cf, char **str);
static CassonFormat *read_casson_struct(char **str);
static Boolean verify_casson(CassonFormat *cf);
static void free_casson_format(CassonFormat *cf);
static Triangulation *casson_to_triangulation(CassonFormat *cf);

static const int vertex_at_faces[4][4] = {
    {9,2,3,1},
    {3,9,0,2},
    {1,3,9,0},
    {2,0,1,9}};

struct CassonFormat
{
    SolutionType    type;
    Boolean         vertices_known;
    int             num_tet;
    EdgeInfo        *head;
};

struct EdgeInfo
{
    int             index,
                    one_vertex,
                    other_vertex,
                    singular_index;
    double          singular_order,
                    e_inner_product,
                    v_inner_product1,
                    v_inner_product2;

    TetEdgeInfo     *head;
    EdgeInfo        *prev,
                    *next;
};

struct TetEdgeInfo
{
    int             tet_index,
	            f1, f2,
	            curves[8];
    double          dihedral_angle;
    TetEdgeInfo     *prev,
	            *next;
};

/* Ported from readCassonFormat in gui/organizer.cpp */
static SolutionType string_to_solution_type(
    char *str)
{
    #define _SOL_TYPE(t)            \
                                    \
        if (strcmp(str, #t) == 0) { \
            return t;               \
        }

    /* Branches similar to line 712... */
    _SOL_TYPE(geometric_solution);
    _SOL_TYPE(nongeometric_solution);
    _SOL_TYPE(flat_solution);
    _SOL_TYPE(degenerate_solution);
    _SOL_TYPE(other_solution);
    _SOL_TYPE(no_solution);
    _SOL_TYPE(step_failed);
    _SOL_TYPE(invalid_solution);

    #undef _SOL_TYPE

    if (strcmp(str, "partially_flat_solution") == 0) {
        return geometric_solution;
    }

    return not_attempted;
}

static Boolean skip_blanks(
    char **str)
{
    while (isspace(**str)) {
        if (**str == '\n') {
            (*str)++;
            return TRUE;
        }
        if (**str == '\0') {
            return TRUE;
        }
        (*str)++;
    }
    return FALSE;
}

/* Ported from readCassonFormat in gui/organizer.cpp */
static Boolean fill_casson_struct(
    CassonFormat *cf,
    char **str)
{
    cf->head = NULL;
    cf->num_tet = 0;
    cf->vertices_known = FALSE;
    cf->type = not_attempted;

    int consumed;

    {
        char solution_type[128];
        if (sscanf(
                *str, " SolutionType %127s%n", solution_type, &consumed) == 1)
        {
            /* Set cf->type based on section */
            /* see branches at 712 */
            cf->type = string_to_solution_type(solution_type);
            *str += consumed;
        }
    }

    if (sscanf(*str, " vertices_known%n", &consumed) == 0) {
        cf->vertices_known = TRUE;
        *str += consumed;
    }

    EdgeInfo *e = NULL;

    while (1)
    {
        if (e) {
            e = e->next = NEW_STRUCT(EdgeInfo);
        } else {
            e = cf->head = NEW_STRUCT(EdgeInfo);
        }
        e->next = NULL;
        e->head = NULL;

        if (sscanf(
                *str,
                " %d %d %lf%n",
                &e->index, &e->singular_index, &e->singular_order,
                &consumed) != 3) {
            return FALSE;
        }

        *str += consumed;

        e->index--;
        e->singular_index--;

        if (cf->vertices_known)
        {
            if (sscanf(
                    *str,
                    " %d %d%n",
                    &e->one_vertex, &e->other_vertex, &consumed) != 2)
            {
                return FALSE;
            }

            *str += consumed;

            if (e->one_vertex > 0) {
                e->one_vertex--;
            }
            if (e->other_vertex > 0) {
                e->other_vertex--;
            }
        }

        TetEdgeInfo * t = NULL;

        while (1)
        {
            if (t) {
                t = t->next = NEW_STRUCT(TetEdgeInfo);
            } else {
                t = e->head = NEW_STRUCT(TetEdgeInfo);
            }

            for(int i = 0; i < 8; i++ ) {
                t->curves[i] = 0;
            }

            char f1, f2;

            if (sscanf(
                    *str, "%d%c%c%n", &t->tet_index, &f1, &f2, &consumed) != 3) {
                return FALSE;
            }

            *str += consumed;

            if (t->tet_index > cf->num_tet) {
                cf->num_tet = t->tet_index;
            }

            t->tet_index--;

            if ('u' <= f1 && f1 <= 'x') {
                t->f1 = f1 - 'u';
            } else {
                return FALSE;
            }

            if ('u' <= f2 && f2 <= 'x') {
                t->f2 = f2 - 'u';
            } else {
                return FALSE;
            }

            if (skip_blanks(str)) {
                break;
            }
        }

        if (skip_blanks(str)) {
            break;
        }
    }

    if (cf->type != not_attempted)
    {
        /* Note that organizer.cpp returns early if the solution
         * type is not attempted.
         * This is not correct since we still need to parse the
         * peripheral curves.
	 * The peripheral curves are always there
	 * (every call to saveTriangulation says so).
	 */

        for (EdgeInfo * e = cf->head; e != NULL; e = e->next)
        {
            int index;
            if (sscanf(
                    *str,
                    " %d %lf %lf %lf%n",
                    &index,
                    &e->e_inner_product,
                    &e->v_inner_product1,
                    &e->v_inner_product2,
                    &consumed) != 4)
            {
                return FALSE;
            }
            *str += consumed;

            for (TetEdgeInfo * t = e->head; t != NULL; t = t->next)
            {
                if (sscanf(*str, " %lf%n", &t->dihedral_angle, &consumed) != 1)
                {
                    return FALSE;
                }
                *str += consumed;
            }

        }

        /* Line 873 in organizer.cpp */
    }

    if (cf->vertices_known)
    {
        for (EdgeInfo * e = cf->head; e != NULL; e = e->next)
        {
            int index;
            if (sscanf(*str, " %d%n", &index, &consumed) != 1) {
                return FALSE;
            }
            *str += consumed;

            for (TetEdgeInfo * t = e->head; t != NULL; t = t->next) {
                for (int i = 0; i < 8; i++) {
                    if (sscanf(*str, " %d%n", &t->curves[i], &consumed) != 1) {
                        return FALSE;
                    }
                    *str += consumed;
                }
            }
        }
    }

    return TRUE;
}

/* Ported from readCassonFormat in gui/organizer.cpp */
static CassonFormat *read_casson_struct(
    char **str)
{
    CassonFormat * cf = NEW_STRUCT(CassonFormat);
    if (fill_casson_struct(cf, str)) {
        return cf;
    }

    free_casson_format(cf);
    return NULL;
}

/* corresponds verify_casson_format in gui/organizer.cpp */
static Boolean verify_casson(CassonFormat *cf)
{
    if (cf == NULL) {
        return FALSE;
    }

    for (int i = 0; i < cf->num_tet; i++) {
	Boolean check[4][4];
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
		check[j][k] = j == k;

                EdgeInfo * ei = cf->head;

                if (ei == NULL) return FALSE;

                while (ei != NULL) {
                    TetEdgeInfo * tei = ei->head;

                    if (tei == NULL) return FALSE;

                    while (tei != NULL) {
                        if (tei->tet_index == i) {
                            if (check[tei->f1][tei->f2]) return TRUE;

                            check[tei->f1][tei->f2] = TRUE;
                            check[tei->f2][tei->f1] = TRUE;
                        }
                        tei = tei->next;
                    }
                    ei = ei->next;
                }

                for (int j = 0; j < 4; j++)
                    for (int k = 0; k < 4; k++)
                        if (check[j][k] == FALSE) return FALSE;
            }
        }
    }

    return TRUE;
}

/* Ported from freeCassonFormat in gui/organizer.cpp. */
void free_casson_format(CassonFormat *cf)
{
    if (cf == NULL)
        return;

    EdgeInfo * e1 = cf->head;

    while (e1 != NULL)
    {
        EdgeInfo * e2 = e1->next;
        TetEdgeInfo * t1 = e1->head;

        while (t1 != NULL)
        {
            TetEdgeInfo * t2 = t1->next;
            my_free(t1);
            t1 = t2;
        }

        my_free(e1);
        e1 = e2;
    }

    my_free(cf);
}

/* Ported from cassonToTriangulation in gui/organizer.cpp */
static Triangulation *casson_to_triangulation(CassonFormat *cf) {
    Triangulation * manifold = NEW_STRUCT(Triangulation);
    initialize_triangulation(manifold);

    manifold->num_tetrahedra = cf->num_tet;
    manifold->solution_type[complete] = not_attempted;
    manifold->solution_type[filled] = not_attempted;
    manifold->num_singular_arcs = 0;
    manifold->num_or_cusps = 0;
    manifold->num_nonor_cusps = 0;
    manifold->num_cusps = 0;

    Tetrahedron ** tet_array = NEW_ARRAY(manifold->num_tetrahedra, Tetrahedron *);
    for (int i = 0; i < manifold->num_tetrahedra; i++) {
        tet_array[i] = NEW_STRUCT(Tetrahedron);
        initialize_tetrahedron(tet_array[i]);
        INSERT_BEFORE(tet_array[i], &manifold->tet_list_end);
    }

    EdgeInfo * ei = cf->head;
    while (ei != NULL) {
        TetEdgeInfo * tei1 = ei->head;
        while (tei1 != NULL) {
	    TetEdgeInfo * tei2;
            if (tei1->next == NULL)
                tei2 = ei->head;
            else
                tei2 = tei1->next;

            int t1 = tei1->tet_index;
            int a1 = tei1->f1;
            int a2 = tei1->f2;
            int a3 = vertex_at_faces[a1][a2];
            int a4 = vertex_at_faces[a2][a1];

            int t2 = tei2->tet_index;
            int b1 = tei2->f1;
            int b2 = tei2->f2;
            int b3 = vertex_at_faces[b1][b2];
            int b4 = vertex_at_faces[b2][b1];

            tet_array[t1]
                ->dihedral_angle[ultimate][edge_between_faces[a1][a2]] =
                tei1->dihedral_angle;

            tet_array[t1]->neighbor[a1] = tet_array[t2];
            tet_array[t2]->neighbor[b2] = tet_array[t1];

            tet_array[t1]->gluing[a1] =
                CREATE_PERMUTATION(a1, b2, a2, b1, a3, b3, a4, b4);

            tet_array[t2]->gluing[b2] =
                CREATE_PERMUTATION(b1, a2, b2, a1, b3, a3, b4, a4);

            tei1 = tei1->next;
        }
        ei = ei->next;
    }

    create_edge_classes(manifold);
    orient_edge_classes(manifold);

    create_cusps(manifold);
    count_cusps(manifold);

    ei = cf->head;

    while (ei != NULL) {
        TetEdgeInfo * tei1 = ei->head;
        int t1 = tei1->tet_index;
        int a1 = tei1->f1;
        int a2 = tei1->f2;

        int index = edge_between_faces[a1][a2];
        EdgeClass *edge = tet_array[t1]->edge_class[index];

        edge->inner_product[ultimate] = ei->e_inner_product;
        edge->inner_product[penultimate] = ei->e_inner_product;

        tet_array[t1]->cusp[remaining_face[a1][a2]]->inner_product[ultimate] =
            ei->v_inner_product1;
        tet_array[t1]
            ->cusp[remaining_face[a1][a2]]
            ->inner_product[penultimate] = ei->v_inner_product1;

        tet_array[t1]->cusp[remaining_face[a2][a1]]->inner_product[ultimate] =
            ei->v_inner_product2;
        tet_array[t1]
            ->cusp[remaining_face[a2][a1]]
            ->inner_product[penultimate] = ei->v_inner_product2;

        if (ei->singular_index < 0) {
            edge->is_singular = FALSE;
            edge->singular_order = 1;
            edge->old_singular_order = 1;
            edge->singular_index = -1;
        } else {
            edge->is_singular = TRUE;
            manifold->num_singular_arcs++;
            edge->singular_order = ei->singular_order;
            edge->old_singular_order = ei->singular_order;
            edge->singular_index = ei->singular_index;
        }

        ei = ei->next;
    }

    if (cf->type != not_attempted)
        for (Tetrahedron * tet = manifold->tet_list_begin.next;
             tet != &manifold->tet_list_end; tet = tet->next) {
            Boolean neg = FALSE;

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    if (i != j) {
                        EdgeClass *edge =
                            tet->edge_class[edge_between_vertices[i][j]];
                        tet->Gram_matrix[i][j] = edge->inner_product[ultimate];

                        if (tet->dihedral_angle[ultimate]
                                               [edge_between_vertices[i][j]] <
                                -0.0001 ||
                            tet->dihedral_angle[ultimate]
                                               [edge_between_vertices[i][j]] >
                                PI + 0.0001)
                            neg = TRUE;
                    } else {
                        Cusp *cusp = tet->cusp[i];
                        tet->Gram_matrix[i][i] = cusp->inner_product[ultimate];
                    }

            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    tet->inverse_Gram_matrix[i][j] =
                        minor1(tet->Gram_matrix, i, j);

            tet->orientation_parameter[ultimate] =
                safe_sqrt(-gl4R_determinant(tet->Gram_matrix));
            if (neg) tet->orientation_parameter[ultimate] *= -1;

            for (int i = 0; i < 6; i++)
                if (cos(tet->dihedral_angle[ultimate][i]) *
                        cos(tet->dihedral_angle[ultimate][i]) <
                    1 / 2)
                    for (int j = 0; j < 2; j++)
                        tet->use_orientation_parameter[j][i] = FALSE;
                else
                    for (int j = 0; j < 2; j++)
                        tet->use_orientation_parameter[j][i] = TRUE;
        }

    identify_cusps(manifold);

    ei = cf->head;

    if (cf->vertices_known)
        while (ei != NULL) {
            TetEdgeInfo * tei1 = ei->head;
            int t1 = tei1->tet_index;
            int a1 = tei1->f1;
            int a2 = tei1->f2;

            tet_array[t1]->cusp[remaining_face[a1][a2]]->index = ei->one_vertex;
            tet_array[t1]->cusp[remaining_face[a2][a1]]->index =
                ei->other_vertex;

            while (tei1 != NULL) {
                t1 = tei1->tet_index;
                a1 = tei1->f1;
                a2 = tei1->f2;
                int top = remaining_face[a1][a2];
                int bottom = remaining_face[a2][a1];

                tet_array[t1]->curve[0][0][top][bottom] = tei1->curves[0];
                tet_array[t1]->curve[0][0][bottom][top] = tei1->curves[1];
                tet_array[t1]->curve[0][1][top][bottom] = tei1->curves[2];
                tet_array[t1]->curve[0][1][bottom][top] = tei1->curves[3];
                tet_array[t1]->curve[1][0][top][bottom] = tei1->curves[4];
                tet_array[t1]->curve[1][0][bottom][top] = tei1->curves[5];
                tet_array[t1]->curve[1][1][top][bottom] = tei1->curves[6];
                tet_array[t1]->curve[1][1][bottom][top] = tei1->curves[7];

                tei1 = tei1->next;
            }

            ei = ei->next;
        }

    orient(manifold);
    my_free(tet_array);

    manifold->solution_type[complete] = cf->type;
    manifold->solution_type[filled] = cf->type;

    if (manifold->solution_type[complete] == geometric_solution)
        my_tilts(manifold);
    peripheral_curves_as_needed(manifold);

    return manifold;
}

/* Ported from Organizer::readTriangulation in gui/organizer.cpp. */
Triangulation * read_casson_format(
    char ** str)
{
    CassonFormat * cf = read_casson_struct(str);
    if (cf == NULL) {
        return NULL;
    }
    Triangulation * t = NULL;
    if (verify_casson(cf)) {
        t = casson_to_triangulation(cf);
    }
    free_casson_format(cf);
    return t;
}

#define NL(f)   (f==0) ? 'u' : ((f==1) ? 'v' : ((f==2) ? 'w' : 'x'))

/* Ported from Console::saveTriangulation in console.cpp */
void write_casson_format_to_stream(
    OStream *stream,
    Triangulation *manifold,
    Boolean include_angular_error,
    Boolean include_geometric_structure_and_cusp_indices,
    Boolean include_peripheral_curves)
{
    Tetrahedron *tet;
    int index;

    for (tet = manifold->tet_list_begin.next, index = 1;
         tet != &manifold->tet_list_end;
	 tet = tet->next, index++)
    {
        tet->index = index;
    }

    if (include_geometric_structure_and_cusp_indices) {
        if (manifold->solution_type[complete] == geometric_solution) {
            if (contains_flat_tetrahedra(manifold) == TRUE)
                ostream_printf(stream,
                               "SolutionType partially_flat_solution\n");
            else
                ostream_printf(stream, "SolutionType geometric_solution\n");
        }

        if (manifold->solution_type[complete] == nongeometric_solution)
            ostream_printf(stream, "SolutionType nongeometric_solution\n");

        if (manifold->solution_type[complete] == not_attempted)
            ostream_printf(stream, "SolutionType not_attempted\n");

        if (manifold->solution_type[complete] == other_solution)
            ostream_printf(stream, "SolutionType other_solution\n");

        if (manifold->solution_type[complete] == step_failed)
            ostream_printf(stream, "SolutionType step_failed\n");

        if (manifold->solution_type[complete] == no_solution)
            ostream_printf(stream, "SolutionType no_solution\n");

        if (manifold->solution_type[complete] == invalid_solution)
            ostream_printf(stream, "SolutionType invalid_solution\n");

        if (manifold->solution_type[complete] == degenerate_solution)
            ostream_printf(stream, "SolutionType degenerate_solution\n");

        if (manifold->solution_type[complete] == flat_solution)
            ostream_printf(stream, "SolutionType flat_solution\n");

        ostream_printf(stream, "vertices_known\n\n");
    }

    EdgeClass *edge;
    for (edge = manifold->edge_list_begin.next, index = 1;
         edge != &manifold->edge_list_end;
	 edge = edge->next, index++)
    {
        ostream_printf(stream, "%3d %2d", index, edge->singular_index + 1);

        ostream_printf(stream, " %04.3f", edge->singular_order);

        PositionedTet ptet0;
        set_left_edge(edge, &ptet0);

        PositionedTet ptet = ptet0;

        if (include_angular_error) {
            double err =
                (edge->singular_order == 0) ? 0 : TWO_PI / edge->singular_order;

            do {
                err -= ptet.tet->dihedral_angle
                           [ultimate]
                           [edge_between_faces[ptet.near_face][ptet.left_face]];

                veer_left(&ptet);
            } while (!same_positioned_tet(&ptet, &ptet0));

            ostream_printf(stream, " %21.16f", err);
        }

        if (include_geometric_structure_and_cusp_indices)
	{
            if (ptet.tet->cusp[remaining_face[ptet.left_face][ptet.near_face]]
                    ->index > -1)
                ostream_printf(stream, " %2d",
                               ptet.tet->cusp[remaining_face[ptet.left_face]
                                                            [ptet.near_face]]
                                       ->index +
                                   1);
            else
                ostream_printf(
                    stream, " %2d",
                    ptet.tet
                        ->cusp[remaining_face[ptet.left_face][ptet.near_face]]
                        ->index);

            if (ptet.tet->cusp[remaining_face[ptet.near_face][ptet.left_face]]
                    ->index > -1)
                ostream_printf(stream, " %2d",
                               ptet.tet->cusp[remaining_face[ptet.near_face]
                                                            [ptet.left_face]]
                                       ->index +
                                   1);
            else
                ostream_printf(
                    stream, " %2d",
                    ptet.tet
                        ->cusp[remaining_face[ptet.near_face][ptet.left_face]]
                        ->index);
        }

        do
	{
            char c = NL(ptet.left_face);
            char d = NL(ptet.near_face);

            ostream_printf(stream, " %2d%c%c", ptet.tet->index, c, d);

            veer_left(&ptet);

        } while (!same_positioned_tet(&ptet, &ptet0));

        ostream_printf(stream, "\n");
    }

    if (include_geometric_structure_and_cusp_indices &&
        manifold->solution_type[complete] != not_attempted)
    {
        ostream_printf(stream, "\n");

        EdgeClass *edge;
        for (edge = manifold->edge_list_begin.next, index = 1;
             edge != &manifold->edge_list_end;
	     edge = edge->next, index++)
	{
            ostream_printf(stream, "%3d", index);

            PositionedTet ptet0;
            set_left_edge(edge, &ptet0);
            PositionedTet ptet = ptet0;

            ostream_printf(stream, " %21.16f", edge->inner_product[ultimate]);

            int top = remaining_face[ptet.left_face][ptet.near_face];
            ostream_printf(stream, " %21.16f",
                           ptet.tet->cusp[top]->inner_product[ultimate]);

            int bottom = remaining_face[ptet.near_face][ptet.left_face];
            ostream_printf(stream, " %21.16f",
                           ptet.tet->cusp[bottom]->inner_product[ultimate]);

            do {
                ostream_printf(
                    stream, " %21.16f",
                    ptet.tet->dihedral_angle
                        [ultimate]
                        [edge_between_faces[ptet.near_face][ptet.left_face]]);
                veer_left(&ptet);
            } while (!same_positioned_tet(&ptet, &ptet0));

            ostream_printf(stream, "\n");
        }
    }

    if (include_geometric_structure_and_cusp_indices &&
        include_peripheral_curves)
    {
        ostream_printf(stream, "\n");

        EdgeClass *edge;

        for (edge = manifold->edge_list_begin.next, index = 1;
             edge != &manifold->edge_list_end;
	     edge = edge->next, index++)
	{
            ostream_printf(stream, "%3d", index);
            PositionedTet ptet0;
            set_left_edge(edge, &ptet0);
            PositionedTet ptet = ptet0;

            do {
                int top = remaining_face[ptet.left_face][ptet.near_face];
                int bottom = remaining_face[ptet.near_face][ptet.left_face];

                ostream_printf(stream, "   %2d %2d %2d %2d %2d %2d %2d %2d",
                               ptet.tet->curve[0][0][top][bottom],
                               ptet.tet->curve[0][0][bottom][top],
                               ptet.tet->curve[0][1][top][bottom],
                               ptet.tet->curve[0][1][bottom][top],
                               ptet.tet->curve[1][0][top][bottom],
                               ptet.tet->curve[1][0][bottom][top],
                               ptet.tet->curve[1][1][top][bottom],
                               ptet.tet->curve[1][1][bottom][top]);
                veer_left(&ptet);
            } while (!same_positioned_tet(&ptet, &ptet0));

            ostream_printf(stream, "\n");
        }
    }
}

char *
write_casson_format_to_string(
    Triangulation * manifold,
    Boolean include_angular_error,
    Boolean include_geometric_structure_and_cusp_indices,
    Boolean include_peripheral_curves)
{
    OStream stream;
    string_stream_init(&stream);

    write_casson_format_to_stream(
        &stream,
        manifold,
        include_angular_error,
        include_geometric_structure_and_cusp_indices,
        include_peripheral_curves);

    return stream.buffer;
}
