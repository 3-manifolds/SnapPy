/**
 *  Symplectic Basis
 *
 *  Computes a symplectic basis of a triangulated knot with orientable torus cusps.
 *  This symplectic matrix extends the Neumann-Zagier matrix to one which is symplectic
 *  up to factors of 2, and which arises from the triangulation of the manifold.
 *
 *  See - https://arxiv.org/abs/2208.06969
 *
 */

#include <stdio.h>
#include <string.h>
#include "SnapPea.h"
#include "kernel.h"

#define ATLEAST_TWO(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))
#define TRI_TO_INDEX(tet_index, tet_vertex)     (4 * (tet_index) + (tet_vertex))

#define COPY_PATH_ENDPOINT(new, old)    {                                                       \
                                            (new)->vertex = (old)->vertex;                      \
                                            (new)->face = (old)->face;                          \
                                            (new)->tri = (old)->tri;                            \
                                            (new)->region_index = (old)->region_index;          \
                                            (new)->region = (old)->region;                      \
                                            (new)->node = (old)->node;                          \
                                            (new)->num_adj_curves = (old)->num_adj_curves;      \
                                        }

#define COPY_PATH_NODE(new, old)        {                                                           \
                                            (new)->next = NULL;                                     \
                                            (new)->prev = NULL;                                     \
                                            (new)->next_face = (old)->next_face;                    \
                                            (new)->prev_face = (old)->prev_face;                    \
                                            (new)->inside_vertex = (old)->inside_vertex;            \
                                            (new)->cusp_region_index = (old)->cusp_region_index;    \
                                            (new)->tri = (old)->tri;                                \
                                        }

enum pos {
    START,
    FINISH
};

static int debug = 0;

/**
 * Queue
 */

typedef struct Queue {
    int                         front;
    int                         rear;
    int                         len;
    int                         size;
    int                         *array;
} Queue ;

/**
 * Graph
 */

typedef struct EdgeNode {
    int                         y;
    struct EdgeNode             *next;
    struct EdgeNode             *prev;
} EdgeNode;

typedef struct Graph {
    EdgeNode                    *edge_list_begin;        /** header node of doubly linked list */
    EdgeNode                    *edge_list_end;          /** tailer node ... */
    int                         *degree;                 /** degree of each vertex */
    int                         *color;                  /** color a tree bipartite */
    int                         num_vertices;            /** number of vertices in the graph */
    Boolean                     directed;                /** is the graph directed */
} Graph;

/**
 * The End multi graph is a quotient of certain Heegard surface by a map which
 * collapses each boundary component of the manifold to a point and collapses
 * each annulus around an edge of the triangulation to an edge.
 */

typedef struct CuspEndPoint {
    int                         cusp_index;
    int                         edge_class[2];
    struct CuspEndPoint         *next;
    struct CuspEndPoint         *prev;
} CuspEndPoint;

typedef struct EndMultiGraph {
    int                         e0;                      /** edge connecting vertices of the same color */
    int                         num_edge_classes;
    int                         num_cusps;
    int                         **edges;                 /** edge_class[u][v] is the edge class of the edge u->v */
    Boolean                     *edge_classes;           /** which edge classes are in the multigraph */
    Graph                       *multi_graph;            /** tree with extra edge of cusps */
} EndMultiGraph;

/**
 * Dual Curves
 *
 * Each oscillating curve contributes combinatorial holonomy, we store this in
 * curve[4][4] in a similar way to the curve[4][4] attribute of a Tetrahedron.
 * An array of size num_edge_classes is attached to each Tetrahedron.
 * tet->extra[edge_class]->curve[v][f] is the intersection number of
 * the oscillating curve associated to edge_class with the face 'f' of the
 * cusp triangle at vertex 'v' of tet.
 */

struct extra {
    int                         curve[4][4];            /** oscillating curve holonomy for a cusp triangle */
};

/**
 * Path End Points
 *
 * Path endpoints can have different states of initialisation. As a convention 
 * if the pointers tri and region are NULL then endpoint is not initialised. 
 * If tri is not NULL and region is NULL then the endpoint is initialised but 
 * the region is not known either because it has not been choosen or we have 
 * split along the curve. In either of the previous cases the tri pointer is 
 * still valid. If tri is not NULL and region is not NULL then the region 
 * pointer is valid and tri = region->tri.
 *
 * Oscillating curves consist of a path through the end multi graph, on each
 * cusp of the path, we have a curve component, which contains a path through
 * the cusp region graph on the cusp, and two path end points, one at either
 * end of the curve.
 */

typedef struct PathEndPoint {
    FaceIndex                   face;                   /** face containg the short rectangle carrying the curve */
    VertexIndex                 vertex;                 /** vertex we dive through the manifold along */
    int                         region_index;           /** index of the region the endpoint lies in */
    int                         num_adj_curves;         /** where the curve dives into the manifold */
    struct PathNode             *node;                  /** pointer to the path node which connects to the endpoint */
    struct CuspRegion           *region;                /** pointer to the region the endpoint lies in */
    struct CuspTriangle         *tri;                   /** pointer to the cusp triangle the endpoint lies in */
} PathEndPoint;

typedef struct PathNode { 
    int                         cusp_region_index;
    FaceIndex                   next_face;               /** face the path crosses to the next node */
    FaceIndex                   prev_face;               /** face the path crosses to the prev node */
    VertexIndex                 inside_vertex;           /** inside vertex of the path */
    struct CuspTriangle         *tri;                    /** cusp triangle the node lies in */
    struct PathNode             *next;                   /** next node in doubly linked list */
    struct PathNode             *prev;
} PathNode;

typedef struct CurveComponent {
    int                         edge_class[2];          /** edge classes at path end points */
    int                         cusp_index;             /** which cusp does the curve lie in */
    PathNode                    path_begin;             /** header node of doubbly linked list */
    PathNode                    path_end;               /** tailer node of ... */
    PathEndPoint                endpoints[2];           /** path end points */
    struct CurveComponent       *next;                  /** next curve component in doubly linked list */
    struct CurveComponent       *prev;                  /** prev ... */
} CurveComponent;

typedef struct OscillatingCurves {
    int                         num_curves;
    int                         *edge_class;
    CurveComponent              *curve_begin;          /** array of doubly linked lists of dual curves */
    CurveComponent              *curve_end;            /** array of ... */
} OscillatingCurves;

/**
 * Cusp Triangulation
 *
 * CuspTriangle stores information about a triangle in the cusp triangulation. 
 * The homology curves bound a fundamental domain, and cusp regions store the 
 * information for intersection of this domain with each cusp triangle. When
 * we add oscillating curves, these regions are divided further.
*/

typedef struct CuspVertex {
    int                         edge_class;
    int                         edge_index;
    EdgeClass                   *edge;
    VertexIndex                 v1;
    VertexIndex                 v2;
} CuspVertex;

typedef struct CuspTriangle {
    Tetrahedron                 *tet;                   /** tetrahedron the triangle comes from */
    Cusp                        *cusp;                  /** cusp the triangle lies in */
    int                         tet_index;              /** tet->index */
    VertexIndex                 tet_vertex;             /** vertex the triangle comes from */
    int                         num_curves;             /** number of curves on the triangle */
    CuspVertex                  vertices[4];            /** information about each vertex */
    struct CuspTriangle         *neighbours[4];         /** triangle neighbouring a face */
    struct CuspTriangle         *next;                  /** next cusp triangle on doubly linked list */
    struct CuspTriangle         *prev;                  /** prev cusp triangle on doubly linkled list */
} CuspTriangle;

typedef struct CuspRegion {
    CuspTriangle                *tri;                   /** cusp triangle the region lies on */
    int                         tet_index;              /** tri->tetIndex */
    VertexIndex                 tet_vertex;             /** tri->tet_vertex */
    int                         index;                  /** index of the cusp region */
    int                         curve[4][4];            /** looking at face, number of curves between the region and vertex */
    Boolean                     adj_cusp_triangle[4];   /** does the region meet this edge of the cusp triangle */
    Boolean                     dive[4][4];             /** can we dive along the face into this vertex */
    int                         num_adj_curves[4][4];   /** stores the number of curves between a region and a face */
    int                         temp_adj_curves[4][4];  /** store the adj curve until pathfinding is complete */
    struct CuspRegion           *adj_cusp_regions[4];   /** index of the adjacent regions */
    struct CuspRegion           *next;                  /** next cusp region in doubly linked list */
    struct CuspRegion           *prev;                  /** prev cusp region in doubly linked list */
} CuspRegion;

typedef struct CuspStructure {
    int                         intersect_tet_index;    /** index of the intersection triangle */
    VertexIndex                 intersect_tet_vertex;   /** vertex of the intersection triangle */
    int                         num_edge_classes;       /** number of edge classes in the cusp */
    int                         num_cusp_triangles;     /** number of cusp triangle in the cusp */
    int                         num_cusp_regions;       /** number of cusp regions in the cusp */
    Triangulation               *manifold;              /** manifold */
    Cusp                        *cusp;                  /** which manifold cusp does the struct lie in */
    Graph                       *cusp_region_graph;     /** dual graph of the cusp region */
    CuspRegion                  **regions;              /** regions in the cusp region graph */
    CuspTriangle                cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle                cusp_triangle_end;      /** tailer node of ... */
    CuspRegion                  *cusp_region_begin;     /** array of header nodes for cusp regions, index by cusp tri */
    CuspRegion                  *cusp_region_end;       /** array of tailer nodes for ...*/
} CuspStructure;

/**
 * Queue Data Structure
 */

Queue                   *init_queue(int);
Queue                   *enqueue(Queue *, int);
int                     dequeue(Queue *);
Queue                   *resize_queue(Queue *);
Boolean                 empty_queue(Queue *);
void                    free_queue(Queue *);

/**
 * Graph
 */

Graph                   *init_graph(int, Boolean);
void                    free_graph(Graph *);
int                     insert_edge(Graph *, int, int, Boolean);
void                    delete_edge(Graph *, int, int, Boolean);
Boolean                 edge_exists(Graph *, int, int);

/**
 * Breadth First Search
 */

void                    init_search(Graph *, Boolean *, Boolean *, int *);
void                    bfs(Graph *, int, Boolean *, Boolean *, int *);
void                    find_path(int, int, int *, EdgeNode *, EdgeNode *);
Boolean                 cycle_exists(Graph *, int, Boolean *, Boolean *, int *, int *, int *);
int                     **ford_fulkerson(Graph *, int, int);
int                     augment_path(Graph *, int **, Boolean *, int, int, int);
int                     bfs_target_list(Graph *, int, int *, int, Boolean *, Boolean *, int *);
Boolean                 contains(int *, int, int);
void                    free_edge_node(EdgeNode *, EdgeNode *);

/**
 * Symplectic Basis
 */

int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     *combinatorial_holonomy(Triangulation *, int);
void                    oscillating_curves(Triangulation *, Boolean *);

/**
 * Initialisation Functions
 */

CuspStructure           *init_cusp_structure(Triangulation *, Cusp *);
void                    free_cusp_structure(CuspStructure **, int, int);
void                    init_cusp_triangulation(Triangulation *, CuspStructure *);
void                    init_cusp_region(CuspStructure *);
int                     init_intersect_cusp_region(CuspStructure *, CuspTriangle *, int);
int                     init_intersect_vertex_two_zero_flows(CuspStructure *,
			    CuspTriangle *, int);
int                     init_normal_cusp_region(CuspStructure *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspStructure *, CuspTriangle *,
			    const int [4], const Boolean [4], int);
void                    init_train_line(CuspStructure *);
CurveComponent          *init_curve_component(int, int, int);
OscillatingCurves       *init_oscillating_curves(Triangulation *, const Boolean *);
void                    free_oscillating_curves(OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, CuspStructure *);

/**
 * Cusp Functions
 */

int                     net_flow_around_vertex(CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *,
					    CuspTriangle *, int);
void                    update_adj_region_data(CuspStructure *);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    construct_cusp_region_dual_graph(CuspStructure *);
void                    log_structs(Triangulation *, CuspStructure **,
				    OscillatingCurves *, const char *);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(CuspStructure **, OscillatingCurves *, EndMultiGraph *);
void                    do_one_oscillating_curve(CuspStructure **, OscillatingCurves *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int, int);
CurveComponent          *setup_first_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
CurveComponent          *setup_last_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
void                    do_curve_component_to_new_edge_class(CuspStructure *, CurveComponent *);
void                    find_single_endpoint(CuspStructure *, PathEndPoint *, int, int);
void                    find_single_matching_endpoint(CuspStructure *, PathEndPoint *, PathEndPoint *);

void                    graph_path_to_dual_curve(CuspStructure *, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    endpoint_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *, PathEndPoint *, int);
void                    interior_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *);

void                    split_cusp_regions_along_path(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_path_len_one(CuspStructure *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_region_path_interior(CuspRegion *, CuspRegion *, PathNode *, int);
void                    split_cusp_region_path_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, PathNode *, int);

void                    update_adj_curve_along_path(CuspStructure **, OscillatingCurves *, int, Boolean);
void                    update_adj_curve_at_endpoint(PathEndPoint *, CurveComponent *, int);
void                    update_adj_curve_on_cusp(CuspStructure *);
void                    update_path_holonomy(CurveComponent *, int);

/**
 * End Multi Graph
 */

EndMultiGraph           *init_end_multi_graph(Triangulation *);
void                    free_end_multi_graph(EndMultiGraph *);
Graph                   *spanning_tree(Graph *, int, int *);
int                     **find_end_multi_graph_edge_classes(EndMultiGraph *, Triangulation *);
int                     find_edge_class(Triangulation *, int, int);
void                    cusp_graph(Triangulation *, Graph *);
void                    color_graph(Graph *);
int                     find_same_color_edge(Triangulation *, EndMultiGraph *, Graph *);
int                     find_path_len(int, int, int *, int);
void                    find_multi_graph_path(Triangulation *, EndMultiGraph *, CuspEndPoint *, CuspEndPoint *, int);
void                    graph_path_to_cusp_path(EndMultiGraph *, EdgeNode *, EdgeNode *, CuspEndPoint *, CuspEndPoint *, int);
void                    find_edge_ends(Graph *, Triangulation *, int, int *, int *);

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};

// -------------------------------------------------

// Queue 

Queue *init_queue(int size) {
    Queue *q    = NEW_STRUCT( Queue );

    q->front    = 0;
    q->rear     = -1;
    q->len      = 0;
    q->size     = MAX(size, 256);
    q->array    = NEW_ARRAY(q->size, int);

    return q;
}

Queue *enqueue(Queue *q, int i) {
    // Queue is full
    if ( q->size == q->len ) {
        q = resize_queue(q);
        q = enqueue(q, i);
    } else {
        q->rear = (q->rear + 1) % q->size;
        q->array[q->rear] = i;
        q->len++;
    }

    return q;
}

int dequeue(Queue *q) {
    // User to verify queue is not empty
    int i = q->array[q->front];

    q->front = (q->front + 1) % q->size;
    q->len--;

    return i;
}

Boolean empty_queue(Queue *q) {
    if (q->len > 0)
        return FALSE;

    return TRUE;
}

Queue *resize_queue(Queue *q) {
    int i;
    Queue *p = init_queue(2 * q->size);

    // Copy elements to new array
    while (!empty_queue(q)) {
        i = dequeue(q);
        enqueue(p, i);
    }

    free_queue(q);
    return p;
}

void free_queue(Queue *q) {
    my_free(q->array);
    my_free(q);
}

// Graph

/*
 * Initialise the arrays of the graph 'g' to their default values
 */

Graph *init_graph(int max_vertices, Boolean directed) {
    int i;
    Graph *g = NEW_STRUCT(Graph);

    g->num_vertices = max_vertices;
    g->directed = directed;

    g->edge_list_begin      = NEW_ARRAY(max_vertices, EdgeNode);
    g->edge_list_end        = NEW_ARRAY(max_vertices, EdgeNode);
    g->degree               = NEW_ARRAY(max_vertices, int);
    g->color                = NEW_ARRAY(max_vertices, int);

    for (i = 0; i < max_vertices; i++) {
        g->degree[i] = 0;
        g->color[i] = -1;

        g->edge_list_begin[i].next     = &g->edge_list_end[i];
        g->edge_list_begin[i].prev     = NULL;
        g->edge_list_end[i].next       = NULL;
        g->edge_list_end[i].prev       = &g->edge_list_begin[i];
    }

    return g;
}

void free_graph(Graph *g) {
    if (g == NULL)
        return;

    for (int i = 0; i < g->num_vertices; i++) {
        free_edge_node(&g->edge_list_begin[i], &g->edge_list_end[i]);
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->color);
    my_free(g);
}

/*
 * Insert an edge into the graph 'g' from vertex x to y.
 */

int insert_edge(Graph *g, int x, int y, Boolean directed) {
    // Ignore edge if it already exists
    if (edge_exists(g, x, y))
        return x;

    EdgeNode *p = NEW_STRUCT( EdgeNode);
    INSERT_AFTER(p, &g->edge_list_begin[x]);
    p->y = y;
    g->degree[x]++;

    if (!directed) {
        insert_edge(g, y, x, TRUE);
    }

    return x;
}

/*
 * Remove the edge from vertex x to vertex y
 */

void delete_edge(Graph *g, int vertex_x, int vertex_y, Boolean directed) {
    EdgeNode *node;

    for (node = g->edge_list_begin[vertex_x].next;
        node != &g->edge_list_end[vertex_x] && node->y != vertex_y;
        node = node->next);

    if (node == &g->edge_list_end[vertex_x])
        return;

    REMOVE_NODE(node)
    my_free(node);

    if (!directed) {
        delete_edge(g, vertex_y, vertex_x, TRUE);
    }
}

/*
 * Check if an edge already exists in the graph
 */

Boolean edge_exists(Graph *g, int v1, int v2) {
    EdgeNode *node = &g->edge_list_begin[v1];

    while ((node = node->next)->next != NULL) {
        if (node->y == v2) {
            return TRUE;
        }
    }

    return FALSE;
}

// ---------------------------------------------------------------

// Breadth First Search

/*
 * Initialise default values for bfs arrays
 */

void init_search(Graph *g, Boolean *processed, Boolean *discovered, int *parent) {
    int i;

    for (i = 0; i < g->num_vertices; i ++) {
        processed[i] = FALSE;
        discovered[i] = FALSE;
        parent[i] = -1;
    }
}

/*
 * Graph search algorithm starting at vertex 'start'.
 */

void bfs(Graph *g, int start, Boolean *processed, Boolean *discovered, int *parent) {
    Queue *q = init_queue(10);
    int v, y;
    EdgeNode *p;

    enqueue(q, start);
    discovered[start] = TRUE;

    while (!empty_queue(q)) {
        v = dequeue(q);
        processed[v] = TRUE;
        p = &g->edge_list_begin[v];

        while ((p = p->next)->next != NULL) {
            y = p->y;

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
}

/*
 * Recover the path through the graph from the parents array and store
 * in the doubly linked list node_begin -> ... -> node_end.
 */

void find_path(int start, int end, int *parents, EdgeNode *node_begin, EdgeNode *node_end) {
    int u;

    if (start != end && parents[end] == -1) {
        uFatalError("find_path", "symplectic_basis");
    }

    u = end;
    while (u != start) {
        EdgeNode *new_node = NEW_STRUCT(EdgeNode);
        new_node->y = u;
        INSERT_AFTER(new_node, node_begin);

        u = parents[u];
    };

    EdgeNode *new_node = NEW_STRUCT(EdgeNode);
    new_node->y = start;
    INSERT_AFTER(new_node, node_begin);
}

void free_edge_node(EdgeNode *node_begin, EdgeNode *node_end) {
    EdgeNode *node;

    while (node_begin->next != node_end) {
        node = node_begin->next;
        REMOVE_NODE(node);
        my_free(node);
    }
}

// ---------------------------------------------------

// Symplectic Basis

/*
 * Allocates arrays for symplectic basis and gluing equations.
 * get_gluing_equations find oscillating curves on the manifold.
 * Constructs return array using gluing_equations_for_edge_class
 * and combinatorial_holonomy
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols, int log) {
    int i, j, k;
    debug = log;
    Boolean *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, Boolean);
    Tetrahedron *tet;

    peripheral_curves(manifold);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        if (tet->extra != NULL)
            uFatalError("oscillating_curves", "symplectic_basis");

        tet->extra = NEW_ARRAY(manifold->num_tetrahedra, Extra);

        for (i = 0; i < manifold->num_tetrahedra; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    tet->extra[i].curve[j][k] = 0;
    }

    // Dual Edge Curves Gamma_i -> symplectic equations
    oscillating_curves(manifold, edge_classes);

    // Construct return array
    *num_rows = 2 * (manifold->num_tetrahedra - manifold->num_cusps);
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i]) {
            continue;
        }

        eqns[2 * j]     = gluing_equations_for_edge_class(manifold, i);
        eqns[2 * j + 1] = combinatorial_holonomy(manifold, i);
        j++;
    }

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        my_free(tet->extra);
        tet->extra = NULL;
    }
    my_free(edge_classes);

    *num_cols = 3 * manifold->num_tetrahedra;
    return eqns;
}

/*
 * Copy of get_gluings_equations.c get_gluing_equations() which finds 
 * the edge gluings equations for a given edge index. Used instead 
 * of get_gluing_equations to ensure we have the correct edge index 
 * and simplify memory management since we don't need all the rows of 
 * the gluing equations matrix.
 */

int *gluing_equations_for_edge_class(Triangulation *manifold, int edgeClass) {
    int *eqns, i, T;
    EdgeClass *edge;
    PositionedTet ptet0, ptet;

    T = manifold->num_tetrahedra;
    eqns = NEW_ARRAY(3 * T, int);

    for (i = 0; i < 3 * T; i++)
        eqns[i] = 0;

    /*
     *  Build edge equations.
     */

    for (edge = manifold->edge_list_begin.next; edge != &manifold->edge_list_end; edge = edge->next) {
        if (edge->index == edgeClass)
            break;
    }

    set_left_edge(edge, &ptet0);
    ptet = ptet0;
    do {
        eqns[3 * ptet.tet->index + edge3_between_faces[ptet.near_face][ptet.left_face]]++;
        veer_left(&ptet);
    } while (same_positioned_tet(&ptet, &ptet0) == FALSE);

    return eqns;
}

/*
 * Construct the symplectic equations from the oscillating curves
 */

int *combinatorial_holonomy(Triangulation *manifold, int edge_class) {
    int v, f, ff;
    int *eqns = NEW_ARRAY(3 * manifold->num_tetrahedra, int);
    Tetrahedron *tet;

    for (int i = 0; i < 3 * manifold->num_tetrahedra; i++) {
        eqns[i] = 0;
    }

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which tet vertex
        for (v = 0; v < 4; v++) {
            // which face
            for (f = 0; f < 4; f++) {
                if (f == v)
                    continue;

                ff = (int) remaining_face[v][f];

                eqns[3 * tet->index + edge3_between_faces[f][ff]]
                    += FLOW(tet->extra[edge_class].curve[v][f], tet->extra[edge_class].curve[v][ff]);
            }
        }
    }

    return eqns;
}

/*
 * Initialise cusp structure on each cusp, construct train lines, construct
 * oscillating curves and store the intersection numbers of each curve with the
 * cusp triangles it enters in tet->extra[edge_class]->curve, in the same fashion
 * as the peripheral curves.
 */

void oscillating_curves(Triangulation *manifold, Boolean *edge_classes) {
    int i;
    label_triangulation_edges(manifold);

    CuspStructure **cusps         = NEW_ARRAY(manifold->num_cusps, CuspStructure *);
    EndMultiGraph *multi_graph    = init_end_multi_graph(manifold);
    Cusp *cusp;

    for (i = 0; i < multi_graph->num_edge_classes; i++)
        edge_classes[i] = multi_graph->edge_classes[i] == TRUE ? FALSE : TRUE;

    edge_classes[multi_graph->e0] = FALSE;

    OscillatingCurves *curves   = init_oscillating_curves(manifold, edge_classes);

    for (i = 0; i < manifold->num_cusps; i++) {
        for (cusp = manifold->cusp_list_begin.next; cusp != &manifold->cusp_list_end && cusp->index != i; cusp = cusp->next);

        if (cusp == &manifold->cusp_list_end)
            uFatalError("oscillating_curves", "symplectic_basis");

        cusps[i] = init_cusp_structure(manifold, cusp);
    }

    if (debug) {
        printf("\n");
        printf("Struct Initialisation\n");
        printf("\n");

        log_structs(manifold, cusps, NULL, "gluing");
        log_structs(manifold, cusps, NULL, "homology");
        log_structs(manifold, cusps, NULL, "edge_indices");
        log_structs(manifold, cusps, NULL, "inside_edge");
        log_structs(manifold, cusps, NULL, "cusp_regions");
    }

    do_oscillating_curves(cusps, curves, multi_graph);

    if (debug) {
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("%d, ", cusps[i]->num_cusp_regions);
        }
        printf("\n");
    }

    free_end_multi_graph(multi_graph);
    free_oscillating_curves(curves);
    free_cusp_structure(cusps, manifold->num_cusps, manifold->num_tetrahedra);
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

// ------------------------------------

/*
 * Initialisation Functions
 */

CuspStructure *init_cusp_structure(Triangulation *manifold, Cusp *cusp) {
    CuspStructure *boundary = NEW_STRUCT(CuspStructure);

    // Invalid cusp topology
    if (cusp->topology == Klein_cusp)
        uFatalError("init_cusp_structure", "symplectic_basis");

    boundary->manifold              = manifold;
    boundary->cusp                  = cusp;
    boundary->num_edge_classes      = manifold->num_tetrahedra;
    boundary->num_cusp_triangles    = 0;
    boundary->num_cusp_regions      = 0;
    boundary->regions    = NULL;

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);

    boundary->cusp_region_graph = NULL;
    construct_cusp_region_dual_graph(boundary);

    return boundary;
}

void free_cusp_structure(CuspStructure **cusps, int num_cusps, int num_edge_classes) {
    int cusp_index;
    CuspTriangle *tri;
    CuspRegion *region;
    CuspStructure *cusp;

    for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
        cusp = cusps[cusp_index];
        // free graph
        free_graph(cusp->cusp_region_graph);
        my_free(cusp->regions);

        // free cusp regions
        for (int i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
            while (cusp->cusp_region_begin[i].next != &cusp->cusp_region_end[i]) {
                region = cusp->cusp_region_begin[i].next;
                REMOVE_NODE(region);
                my_free(region);
            }
        }

        my_free(cusp->cusp_region_begin);
        my_free(cusp->cusp_region_end);

        // free cusp triangle
        while (cusp->cusp_triangle_begin.next != &cusp->cusp_triangle_end) {
            tri = cusp->cusp_triangle_begin.next;
            REMOVE_NODE(tri);
            my_free(tri);
        }

        my_free(cusp);
    }

    my_free(cusps);
}

/*
 * Construct the cusp triangle doubly linked list which consists of the
 * triangles in the cusp triangulation
 */

void init_cusp_triangulation(Triangulation *manifold, CuspStructure *cusp) {
    int index = 0;
    VertexIndex vertex;
    FaceIndex face;
    Tetrahedron *tet;
    CuspTriangle *tri;

    // Allocate Cusp Triangulation Header and Tail Null nodes
    cusp->cusp_triangle_begin.next      = &cusp->cusp_triangle_end;
    cusp->cusp_triangle_begin.prev      = NULL;
    cusp->cusp_triangle_end.next        = NULL;
    cusp->cusp_triangle_end.prev        = &cusp->cusp_triangle_begin;

    // which tetrahedron are we on
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // while vertex are we on
        for (vertex = 0; vertex < 4; vertex++) {
            // is this vertex on the right cusp
            if (tet->cusp[vertex] != cusp->cusp) {
                continue;
            }

            tri = NEW_STRUCT( CuspTriangle );
            INSERT_BEFORE(tri, &cusp->cusp_triangle_end);
            index++;

            tri->tet = tet;
            tri->cusp = tet->cusp[vertex];
            tri->tet_index = tri->tet->index;
            tri->tet_vertex = vertex;

            tri->num_curves = net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][0])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][1])
                              + net_flow_around_vertex(tri, edgesThreeToFour[tri->tet_vertex][2]);

            for (face = 0; face < 4; face ++) {
                if (tri->tet_vertex == face)
                    continue;

                tri->vertices[face].v1              = tri->tet_vertex;
                tri->vertices[face].v2              = face;
                tri->vertices[face].edge            = tri->tet->edge_class[
                        edge_between_vertices[tri->vertices[face].v1][tri->vertices[face].v2]];
                tri->vertices[face].edge_class      = tri->vertices[face].edge->index;
                tri->vertices[face].edge_index      = -1;
            }
        }
    }

    // which cusp triangle
    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        // which vertex
        for (face = 0; face < 4; face++) {
            if (face == tri->tet_vertex)
                continue;

            tri->neighbours[face] = find_cusp_triangle(&cusp->cusp_triangle_begin, &cusp->cusp_triangle_end, tri, face);
        }
    }

    label_cusp_vertex_indices(&cusp->cusp_triangle_begin, &cusp->cusp_triangle_end, cusp->num_edge_classes);
    cusp->num_cusp_triangles = index;
}

/*
 * Initialise the cusp region doubly linked list to cotain the regions bounded 
 * by the meridian and longitude curves.
 */

void init_cusp_region(CuspStructure *cusp) {
    int index;
    CuspTriangle *tri;

    // Header and tailer nodes.
    cusp->cusp_region_begin = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);
    cusp->cusp_region_end   = NEW_ARRAY(4 * cusp->manifold->num_tetrahedra, CuspRegion);

    for (index = 0; index < 4 * cusp->manifold->num_tetrahedra; index++) {
        cusp->cusp_region_begin[index].next    = &cusp->cusp_region_end[index];
        cusp->cusp_region_begin[index].prev    = NULL;
        cusp->cusp_region_end[index].next      = NULL;
        cusp->cusp_region_end[index].prev      = &cusp->cusp_region_begin[index];
    }

    index = 0;
    for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
        // Intersection vertex doesn't have a center
        if (tri->tet_index == cusp->intersect_tet_index && tri->tet_vertex == cusp->intersect_tet_vertex) {
            index = init_intersect_cusp_region(cusp, tri, index);
            continue;
        }

        index = init_normal_cusp_region(cusp, tri, index);
    }

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

/*
 * Assume peripheral_curves() has been called, and as a result the only curves 
 * on the intersection triangle are those which intersect, and they give a 
 * valid intersection.
 */

int init_intersect_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, curve_index, vertex, v1, v2, v3;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 1; i < net_flow_around_vertex(tri, vertex); i++) {
            for (curve_index = 0; curve_index < 2; curve_index++) {
                distance[v1]                    = i;
                distance[v2]                    = MIN(distance[v1], 2 * net_flow_around_vertex(tri, vertex) - distance[v1])
                                                + net_flow_around_vertex(tri, v2) + net_flow_around_vertex(tri, v1);
                distance[vertex]                = net_flow_around_vertex(tri, vertex)
                                                - distance[v1] + net_flow_around_vertex(tri, v1);
                distance[tri->tet_vertex]       = -1;

                adj_triangle[v1]                = 1;
                adj_triangle[v2]                = 0;
                adj_triangle[vertex]            = 0;
                adj_triangle[tri->tet_vertex]   = -1;

                set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
                index++;

                // Swap vertices
                v1 = (int) remaining_face[vertex][tri->tet_vertex];
                v2 = (int) remaining_face[tri->tet_vertex][vertex];
            }
        }

        // Region in the middle of face vertex
        if (net_flow_around_vertex(tri, v1) && net_flow_around_vertex(tri, v2)) {
            distance[v1]                    = net_flow_around_vertex(tri, v2);
            distance[v2]                    = net_flow_around_vertex(tri, v1);
            distance[vertex]                = MIN(net_flow_around_vertex(tri, v1) + distance[v1],
                                              net_flow_around_vertex(tri, v2) + distance[v2])
                                                      + net_flow_around_vertex(tri, vertex);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[v1]                = 0;
            adj_triangle[v2]                = 0;
            adj_triangle[vertex]            = 1;
            adj_triangle[tri->tet_vertex]   = -1;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }
    }

    // Region of distance 0 to vertex
    v1 = edgesThreeToFour[tri->tet_vertex][0];
    v2 = edgesThreeToFour[tri->tet_vertex][1];
    v3 = edgesThreeToFour[tri->tet_vertex][2];

    // Edge Case: Two vertices with 0 flow
    if (ATLEAST_TWO(!net_flow_around_vertex(tri, v1),
                    !net_flow_around_vertex(tri, v2),
                    !net_flow_around_vertex(tri, v3)))
        return init_intersect_vertex_two_zero_flows(cusp, tri, index);

    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        distance[vertex]               = 0;
        distance[v1]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v1);
        distance[v2]                   = net_flow_around_vertex(tri, vertex) + net_flow_around_vertex(tri, v2);
        distance[tri->tet_vertex]      = -1;

        adj_triangle[vertex]           = 0;
        adj_triangle[v1]               = 1;
        adj_triangle[v2]               = 1;
        adj_triangle[tri->tet_vertex]  = 0;

        set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
        index++;
    }

    return index;
}

int init_intersect_vertex_two_zero_flows(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int vertex, v1, v2, v3, distance[4];
    Boolean adj_triangle[4];

    v1 = (int) edgesThreeToFour[tri->tet_vertex][0];
    v2 = (int) edgesThreeToFour[tri->tet_vertex][1];
    v3 = (int) edgesThreeToFour[tri->tet_vertex][2];

    distance[v1]                   = net_flow_around_vertex(tri, v1);
    distance[v2]                   = net_flow_around_vertex(tri, v2);
    distance[v3]                   = net_flow_around_vertex(tri, v3);
    distance[tri->tet_vertex]      = -1;

    adj_triangle[v1]               = 1;
    adj_triangle[v2]               = 1;
    adj_triangle[v3]               = 1;
    adj_triangle[tri->tet_vertex]  = -1;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;

    // Find vertex with non-zero flow
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        if (net_flow_around_vertex(tri, vertex)) {
            v1 = vertex;
            v2 = (int) remaining_face[tri->tet_vertex][v1];
            v3 = (int) remaining_face[v1][tri->tet_vertex];
            break;
        }
    }
    distance[v1]                    = 0;
    distance[v2]                    = net_flow_around_vertex(tri, v1);
    distance[v3]                    = net_flow_around_vertex(tri, v1);
    distance[tri->tet_vertex]       = -1;

    adj_triangle[v1]                = 0;
    adj_triangle[v2]                = 1;
    adj_triangle[v3]                = 1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);

    return index + 1;
}

int init_normal_cusp_region(CuspStructure *cusp, CuspTriangle *tri, int index) {
    int i, vertex, v1, v2;
    int distance[4];
    Boolean adj_triangle[4];

    // which vertex are we inside the flow of
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex) {
            continue;
        }

        v1 = (int) remaining_face[tri->tet_vertex][vertex];
        v2 = (int) remaining_face[vertex][tri->tet_vertex];

        for (i = 0; i < net_flow_around_vertex(tri, vertex); i++) {
            distance[vertex]                = i;
            distance[v1]                    = net_flow_around_vertex(tri, v1)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[v2]                    = net_flow_around_vertex(tri, v2)
                                            + (net_flow_around_vertex(tri, vertex) - distance[vertex]);
            distance[tri->tet_vertex]       = -1;

            adj_triangle[vertex]            = 0;
            adj_triangle[v1]                = 1;
            adj_triangle[v2]                = 1;
            adj_triangle[tri->tet_vertex]   = 0;

            set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
            index++;
        }

    }

    // center region
    for (vertex = 0; vertex < 4; vertex++) {
        if (vertex == tri->tet_vertex)
            continue;

        distance[vertex]        = net_flow_around_vertex(tri, vertex);
        adj_triangle[vertex]    = 1;
    }

    distance[tri->tet_vertex]       = -1;
    adj_triangle[tri->tet_vertex]   = 0;

    set_cusp_region_data(cusp, tri, distance, adj_triangle, index);
    index++;
    return index;
}

/*
 * Helper function to init_cusp_regions which allocates the attributes of the 
 * cusp region
 */

void set_cusp_region_data(CuspStructure *cusp, CuspTriangle *tri, const int distance[4],
                          const Boolean adj_cusp_triangle[4], int index) {
    int i, j, v1, v2, v3;
    CuspRegion *region = NEW_STRUCT( CuspRegion );
    INSERT_BEFORE(region, &cusp->cusp_region_end[TRI_TO_INDEX(tri->tet_index, tri->tet_vertex)]);

    region->tri             = tri;
    region->tet_index       = region->tri->tet_index;
    region->tet_vertex      = region->tri->tet_vertex;
    region->index           = index;

    // default values
    for (i = 0; i < 4; i++) {
        region->adj_cusp_triangle[i] = FALSE;
        region->adj_cusp_regions[i]  = NULL;

        for (j = 0; j < 4; j++) {
            region->curve[i][j]             = -1;
            region->dive[i][j]              = 0;
            region->num_adj_curves[i][j]    = 0;
            region->temp_adj_curves[i][j]   = 0;
        }
    }

    for (i = 0; i < 3; i++) {
        v1 = edgesThreeToFour[tri->tet_vertex][i];
        v2 = edgesThreeToFour[tri->tet_vertex][(i + 1) % 3];
        v3 = edgesThreeToFour[tri->tet_vertex][(i + 2) % 3];

        region->curve[v2][v1]   = distance[v1];
        region->curve[v3][v1]   = distance[v1];
        region->dive[v2][v1]    = distance[v1] ? FALSE : TRUE;
        region->dive[v3][v1]    = distance[v1] ? FALSE : TRUE;

        region->adj_cusp_triangle[v1] = adj_cusp_triangle[v1];
    }
}

CurveComponent *init_curve_component(int edge_class_start, int edge_class_finish, int cusp_index) {
    int i;

    CurveComponent *path = NEW_STRUCT(CurveComponent );

    path->path_begin.next     = &path->path_end;
    path->path_begin.prev     = NULL;
    path->path_end.next       = NULL;
    path->path_end.prev       = &path->path_begin;

    path->edge_class[START]     = edge_class_start;
    path->edge_class[FINISH]    = edge_class_finish;
    path->cusp_index            = cusp_index;

    for (i = 0; i < 2; i++) {
        path->endpoints[i].tri              = NULL;
        path->endpoints[i].region           = NULL;
        path->endpoints[i].num_adj_curves   = 0;
    }

    return path;
}

/*
 * Initialise dual curve doubly linked list which stores the oscillating curves
 * on the cusp
 */

OscillatingCurves *init_oscillating_curves(Triangulation *manifold, const Boolean *edge_classes) {
    int i, j;
    OscillatingCurves *curves = NEW_STRUCT(OscillatingCurves );

    curves->num_curves = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++)
        if (edge_classes[i])
            curves->num_curves++;

    curves->curve_begin               = NEW_ARRAY(curves->num_curves, CurveComponent );
    curves->curve_end                 = NEW_ARRAY(curves->num_curves, CurveComponent );
    curves->edge_class                = NEW_ARRAY(curves->num_curves, int);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i])
            continue;

        curves->edge_class[j] = i;
        j++;
    }

    // which curve
    for (i = 0; i < curves->num_curves; i++) {
        curves->curve_begin[i].next    = &curves->curve_end[i];
        curves->curve_begin[i].prev    = NULL;
        curves->curve_end[i].next      = NULL;
        curves->curve_end[i].prev      = &curves->curve_begin[i];
    }

    return curves;
}

void free_oscillating_curves(OscillatingCurves *curves) {
    int i;
    CurveComponent *path;
    PathNode *path_node;

    for (i = 0; i < curves->num_curves; i++) {
        while (curves->curve_begin[i].next != &curves->curve_end[i]) {
            path = curves->curve_begin[i].next;
            REMOVE_NODE(path);

            while (path->path_begin.next != &path->path_end) {
                path_node = path->path_begin.next;
                REMOVE_NODE(path_node);
                my_free(path_node);
            }

            my_free(path);
        }
    }

    my_free(curves->curve_begin);
    my_free(curves->curve_end);
    my_free(curves->edge_class);
    my_free(curves);
}

// ----------------------------------------------------

/*
 * Cusp Functions
 */

/*
 * peripheral_curves.c places a meridian and longitude curve on each cusp. It
 * starts at a base triangle, the intersection point, and searches outwards.
 * Note it does not visit a cusp triangle more than once. So we find a cusp
 * triangle which contains both a meridian and longitude (this should be the
 * same intersection triangle that peripheral_curves sets since it is the same
 * search process) and assert this is the intersection triangle. Currently
 * init_cusp_regions assumes the intersection triangle only contains curves
 * which intersect. This is because we need some information about the curves
 * to construct the cusp regions.
 */

void find_intersection_triangle(Triangulation *manifold, CuspStructure *boundary) {
    FaceIndex   face;
    Cusp *cusp = boundary->cusp;
    int n;

    for (cusp->basepoint_tet = manifold->tet_list_begin.next;
         cusp->basepoint_tet != &manifold->tet_list_end;
         cusp->basepoint_tet = cusp->basepoint_tet->next)

        for (cusp->basepoint_vertex = 0;
             cusp->basepoint_vertex < 4;
             cusp->basepoint_vertex++)
        {
            if (cusp->basepoint_tet->cusp[cusp->basepoint_vertex] != cusp)
                continue;

            for (face = 0; face < 4; face++)
            {
                if (face == cusp->basepoint_vertex)
                    continue;

                for (n = 0; n < 2; n++) {
                    cusp->basepoint_orientation = ORIENTATION(n);

                    if (cusp->basepoint_tet->curve
                        [M]
                        [cusp->basepoint_orientation]
                        [cusp->basepoint_vertex]
                        [face] != 0
                        && cusp->basepoint_tet->curve
                           [L]
                           [cusp->basepoint_orientation]
                           [cusp->basepoint_vertex]
                           [face] != 0) {
                        /*
                         *  We found the basepoint!
                         */

                        boundary->intersect_tet_index  = cusp->basepoint_tet->index;
                        boundary->intersect_tet_vertex = cusp->basepoint_vertex;
                        return;
                    }


                }
            }
        }
}

/*
 * Calculate the number of curves passing around a vertex in the cusp 
 * triangulation.
 */

int net_flow_around_vertex(CuspTriangle *tri, int vertex) {
    int mflow, lflow, retval;

    // Contribution from meridian curves
    mflow = FLOW(tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[M][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    // Contribution from longitudinal curves
    lflow = FLOW(tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[tri->tet_vertex][vertex]],
                 tri->tet->curve[L][right_handed][tri->tet_vertex][remaining_face[vertex][tri->tet_vertex]]);

    retval = ABS(mflow) + ABS(lflow);
    return retval;
}

/*
 * Returns a pointer to the cusp triangle which is the neighbour of tri across 
 * face 'face'.
 */

CuspTriangle *find_cusp_triangle(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end,
        CuspTriangle *tri, int face) {
    int tet_index, tet_vertex;
    CuspTriangle *pTri;

    tet_index = tri->tet->neighbor[face]->index;
    tet_vertex = EVALUATE(tri->tet->gluing[face], tri->tet_vertex);

    for (pTri = cusp_triangle_begin->next; pTri != cusp_triangle_end; pTri = pTri->next) {
        if (pTri->tet_index == tet_index && pTri->tet_vertex == tet_vertex)
            return pTri;
    }

    // Didn't find a neighbour
    return NULL;
}

/*
 * Give each edge of the triangulation an index to identify the cusp vertices
 */

void label_triangulation_edges(Triangulation *manifold) {
    int i = 0;
    EdgeClass *edge = &manifold->edge_list_begin;

    while ((edge = edge->next)->next != NULL)
        edge->index = i++;

    // incorrect number of edge classes
    if (i != manifold->num_tetrahedra)
        uFatalError("label_triangulation_edges", "symplectic_basis");
}

/*
 * Each edge class of the manifold appears as two vertices in the cusp
 * triangulation. We iterate over the cusp triangulation, walking around each
 * vertex to give it the same index.
 */

void label_cusp_vertex_indices(CuspTriangle *cusp_triangle_begin, CuspTriangle *cusp_triangle_end, int numEdgeClasses) {
    int i, vertex;
    CuspTriangle *tri;

    int *current_index = NEW_ARRAY(numEdgeClasses, int);

    for (i = 0; i < numEdgeClasses; i++)
        current_index[i] = 0;

    for (tri = cusp_triangle_begin->next; tri != cusp_triangle_end; tri = tri->next) {
        for (vertex = 0; vertex < 4; vertex++) {
            if (vertex == tri->tet_vertex || tri->vertices[vertex].edge_index != -1)
                continue;

            walk_around_cusp_vertex(tri, vertex, current_index[tri->vertices[vertex].edge_class]);
            current_index[tri->vertices[vertex].edge_class]++;
        }
    }

    my_free(current_index);
}

/*
 * Walk around vertex cusp_vertex of triangle *tri and set edge_index to index.
 */

void walk_around_cusp_vertex(CuspTriangle *tri, int cusp_vertex, int index) {
    int gluing_vertex, outside_vertex, old_gluing_vertex, old_cusp_vertex, old_outside_vertex;
    gluing_vertex = (int) remaining_face[cusp_vertex][tri->tet_vertex];
    outside_vertex = (int) remaining_face[tri->tet_vertex][cusp_vertex];

    while (tri->vertices[cusp_vertex].edge_index == -1) {
        tri->vertices[cusp_vertex].edge_index = index;

        // Move to the next cusp triangle
        old_cusp_vertex         = cusp_vertex;
        old_gluing_vertex       = gluing_vertex;
        old_outside_vertex      = outside_vertex;

        cusp_vertex             = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_cusp_vertex);
        gluing_vertex           = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_outside_vertex);
        outside_vertex          = EVALUATE(tri->tet->gluing[old_gluing_vertex], old_gluing_vertex);
        tri                     = tri->neighbours[old_gluing_vertex];
    }
}

/*
 * Calculate which regions are located across cusp edges and store the result
 * in the adj_cusp_regions attribute
 */

void update_adj_region_data(CuspStructure *cusp) {
    CuspTriangle *adj_triangle;
    CuspRegion *region;
    FaceIndex f;
    int i, adj_index;

    // Add adjacent region info
    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            for (f = 0; f < 4; f++) {
                if (!region->adj_cusp_triangle[f] || region->tet_vertex == f) {
                    region->adj_cusp_regions[f] = NULL;
                    continue;
                }

                adj_triangle = region->tri->neighbours[f];
                adj_index = TRI_TO_INDEX(adj_triangle->tet_index, adj_triangle->tet_vertex);
                region->adj_cusp_regions[f] = find_adj_region(&cusp->cusp_region_begin[adj_index],
                                                              &cusp->cusp_region_end[adj_index],
                                                              region, f);
            }
        }
    }
}

/*
 * Find the cusp region which is adjacent to x across face.
 */

CuspRegion *find_adj_region(CuspRegion *cusp_region_begin, CuspRegion *cusp_region_end,
                            CuspRegion *x, int face) {
    int v1, v2, y_vertex1, y_vertex2, y_face, distance_v1, distance_v2, tet_index, tet_vertex;
    Boolean adj_face;
    CuspTriangle *tri = x->tri;
    CuspRegion *region;

    v1 = (int) remaining_face[tri->tet_vertex][face];
    v2 = (int) remaining_face[face][tri->tet_vertex];

    y_vertex1    = EVALUATE(tri->tet->gluing[face], v1);
    y_vertex2    = EVALUATE(tri->tet->gluing[face], v2);
    y_face       = EVALUATE(tri->tet->gluing[face], face);

    // Check current adj region first
    if (x->adj_cusp_regions[face] != NULL) {
        distance_v1      = (x->curve[face][v1] == x->adj_cusp_regions[face]->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == x->adj_cusp_regions[face]->curve[y_face][y_vertex2]);
        adj_face         = x->adj_cusp_regions[face]->adj_cusp_triangle[y_face];

        if (distance_v1 && distance_v2 && adj_face)
            return x->adj_cusp_regions[face];
    }

    /*
     * We search through the regions in reverse as the new regions
     * are added to the end of the doubly linked list
     */
    for (region = cusp_region_end->prev; region != cusp_region_begin; region = region->prev) {
        tet_index    = (tri->neighbours[face]->tet_index == region->tet_index);
        tet_vertex   = (tri->neighbours[face]->tet_vertex == region->tet_vertex);

        if (!tet_index || !tet_vertex)
            continue;

        distance_v1      = (x->curve[face][v1] == region->curve[y_face][y_vertex1]);
        distance_v2      = (x->curve[face][v2] == region->curve[y_face][y_vertex2]);
        adj_face         = region->adj_cusp_triangle[y_face];

        // missing distance
        if (region->curve[y_face][y_vertex1] == -1 || region->curve[y_face][y_vertex2] == -1)
            uFatalError("find_adj_region", "symplectic_basis");

        if (distance_v1 && distance_v2 && adj_face)
            return region;
    }

    // We didn't find a cusp region
    //uFatalError("find_cusp_region", "symplectic_basis");
    return NULL;
}

/*
 * region1 splits into region1 and region2, set them up to be split
 */

void copy_region(CuspRegion *region1, CuspRegion *region2) {
    int i, j;

    if (region1 == NULL || region2 == NULL || region1->tri == NULL)
        uFatalError("copy_region", "symplectic_basis");

    region2->tri            = region1->tri;
    region2->tet_index      = region1->tet_index;
    region2->tet_vertex     = region1->tet_vertex;

    for (i = 0; i < 4; i++) {
        region2->adj_cusp_triangle[i]   = region1->adj_cusp_triangle[i];
        region2->adj_cusp_regions[i]    = NULL;

        for (j = 0; j < 4; j++) {
            region2->curve[i][j]            = region1->curve[i][j];
            region2->dive[i][j]             = FALSE;
            region2->num_adj_curves[i][j]   = region1->num_adj_curves[i][j];
            region2->temp_adj_curves[i][j]  = region1->temp_adj_curves[i][j];
        }
    }
}

/*
 * Construct the graph with edges coming from adjacent regions, using
 * region->index to label each vertex.
 */

void construct_cusp_region_dual_graph(CuspStructure *cusp) {
    int i, face;
    CuspRegion *region;

    Graph *graph1 = init_graph(cusp->num_cusp_regions, FALSE);

    int *visited = NEW_ARRAY(graph1->num_vertices, int);

    my_free(cusp->regions);
    cusp->regions = NEW_ARRAY(graph1->num_vertices, CuspRegion *);

    for (i = 0; i < graph1->num_vertices; i++) {
        visited[i] = FALSE;
        cusp->regions[i] = NULL;
    }

    // Walk around the cusp triangulation inserting edges
    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            if (visited[region->index])
                continue;

            for (face = 0; face < 4; face++) {
                if (!region->adj_cusp_triangle[face])
                    continue;

                // Missing adj region data
                if (region->adj_cusp_regions[face] == NULL)
                    uFatalError("construct_cusp_region_dual_graph", "symplectic_basis");

                insert_edge(graph1, region->index, region->adj_cusp_regions[face]->index, graph1->directed);
                cusp->regions[region->index] = region;
            }

            visited[region->index] = 1;
        }
    }

    free_graph(cusp->cusp_region_graph);
    my_free(visited);

    cusp->cusp_region_graph = graph1;
}

/*
 * Types: gluing, train_lines, cusp_regions, homology, edge_indices,
 * dual_curves, inside_edge, graph, endpoints
 */

void log_structs(Triangulation *manifold, CuspStructure **cusps,
		 OscillatingCurves *curves, const char *type) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    CuspTriangle *tri;
    CuspRegion *region;
    EdgeNode *edge_node;
    PathNode *path_node;
    CurveComponent *path;
    Graph *g;
    CuspStructure *cusp;

    if (strcmp(type, "gluing") == 0) {
        printf("Triangle gluing info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);
            cusp = cusps[i];

            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                for (j = 0; j < 4; j++) {
                    if (j == tri->tet_vertex)
                        continue;

                    x_vertex1 = (int) remaining_face[tri->tet_vertex][j];
                    x_vertex2 = (int) remaining_face[j][tri->tet_vertex];
                    y_vertex1 = EVALUATE(tri->tet->gluing[j], x_vertex1);
                    y_vertex2 = EVALUATE(tri->tet->gluing[j], x_vertex2);

                    printf("    (Tet Index: %d, Tet Vertex: %d) Cusp Edge %d glues to "
                           "(Tet Index: %d, Tet Vertex: %d) Cusp Edge %d. (%d -> %d, %d -> %d)\n",
                           tri->tet_index,               // Tet Index
                           tri->tet_vertex,                // Tet Vertex
                           j,      // Cusp Edge
                           tri->tet->neighbor[j]->index,                              // Tet Index
                           EVALUATE(tri->tet->gluing[j], tri->tet_vertex),             // Tet Vertex
                           EVALUATE(tri->tet->gluing[j], j),   // Cusp Edge
                           x_vertex1, y_vertex1,
                           x_vertex2, y_vertex2
                    );
                }
            }
        }
    } else if (strcmp(type, "cusp_regions") == 0) {
        printf("Cusp Region info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (j = 0; j < 4 * cusp->manifold->num_tetrahedra; j++) {
                printf("    Cusp Triangle (Tet Index %d Tet Vertex %d)\n", j / 4, j % 4);
                for (region = cusp->cusp_region_begin[j].next;
                     region != &cusp->cusp_region_end[j]; region = region->next) {
                    v1 = edgesThreeToFour[region->tet_vertex][0];
                    v2 = edgesThreeToFour[region->tet_vertex][1];
                    v3 = edgesThreeToFour[region->tet_vertex][2];

                    printf("    Region %d (Tet Index: %d, Tet Vertex: %d) (Adj Tri: %d, %d, %d) (Adj Regions: %d, %d, %d) "
                           " (Curves: [%d %d] [%d %d] [%d %d]) (Adj Curves: [%d %d] [%d %d] [%d %d]) (Dive: [%d %d] [%d %d] [%d %d])\n",
                           region->index, region->tet_index, region->tet_vertex,
                           region->adj_cusp_triangle[v1], region->adj_cusp_triangle[v2], region->adj_cusp_triangle[v3],
                           region->adj_cusp_regions[v1] == NULL ? -1 : region->adj_cusp_regions[v1]->index,
                           region->adj_cusp_regions[v2] == NULL ? -1 : region->adj_cusp_regions[v2]->index,
                           region->adj_cusp_regions[v3] == NULL ? -1 : region->adj_cusp_regions[v3]->index,
                           region->curve[v2][v1], region->curve[v3][v1],
                           region->curve[v1][v2], region->curve[v3][v2],
                           region->curve[v1][v3], region->curve[v2][v3],
                           region->num_adj_curves[v2][v1], region->num_adj_curves[v3][v1],
                           region->num_adj_curves[v1][v2], region->num_adj_curves[v3][v2],
                           region->num_adj_curves[v1][v3], region->num_adj_curves[v2][v3],
                           region->dive[v2][v1], region->dive[v3][v1],
                           region->dive[v1][v2], region->dive[v3][v2],
                           region->dive[v1][v3], region->dive[v2][v3]
                    );
                }
            }
        }

    } else if (strcmp(type, "homology") == 0) {
        printf("Homology info\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            cusp = cusps[i];

            printf("Boundary %d\n", i);
            printf("Intersect Tet Index %d, Intersect Tet Vertex %d\n", cusp->intersect_tet_index, cusp->intersect_tet_vertex);
            printf("    Meridian\n");

            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[M][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[M][right_handed][tri->tet_vertex][3]
                );
            }
            printf("    Longitude\n");
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("        (Tet Index: %d, Tet Vertex: %d) %d %d %d %d\n",
                       tri->tet_index,
                       tri->tet_vertex,
                       tri->tet->curve[L][right_handed][tri->tet_vertex][0],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][1],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][2],
                       tri->tet->curve[L][right_handed][tri->tet_vertex][3]
                );
            }
        }

    } else if (strcmp(type, "edge_indices") == 0) {
        printf("Edge classes\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                v1 = edgesThreeToFour[tri->tet_vertex][0];
                v2 = edgesThreeToFour[tri->tet_vertex][1];
                v3 = edgesThreeToFour[tri->tet_vertex][2];

                printf("    (Tet Index: %d, Tet Vertex: %d) Vertex %d: (%d %d), "
                       "Vertex %d: (%d %d), Vertex %d: (%d %d)\n",
                       tri->tet_index, tri->tet_vertex,
                       v1, tri->vertices[v1].edge_class, tri->vertices[v1].edge_index,
                       v2, tri->vertices[v2].edge_class, tri->vertices[v2].edge_index,
                       v3, tri->vertices[v3].edge_class, tri->vertices[v3].edge_index
                );
            }
        }
    } else if (strcmp(type, "dual_curves") == 0) {
        printf("Oscillating curve paths\n");

        // which dual curve
        for (i = 0; i < curves->num_curves; i++) {
            j = 0;

            printf("Dual Curve %d\n", i);
            // which curve component
            for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
                printf("    Part %d: \n", j);

                for (path_node = path->path_begin.next;
                     path_node != &path->path_end;
                     path_node = path_node->next)
                    printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n",
                           path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex,
                           path_node->next_face, path_node->prev_face, path_node->inside_vertex
                    );
                j++;
            }
        }
    } else if (strcmp(type, "inside_edge") == 0) {
        printf("Inside edge info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
                printf("    (Tet Index: %d, Tet Vertex: %d) Edge label (%d, %d, %d)\n",
                       tri->tet_index,               // Tet Index
                       tri->tet_vertex,                // Tet Vertex
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][1]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][2]],
                       edge3_between_faces[edgesThreeToFour[tri->tet_vertex][0]][edgesThreeToFour[tri->tet_vertex][1]]
                );
            }
        }
    } else if (strcmp(type, "graph") == 0) {
        printf("Graph info\n");

        for (i = 0; i < manifold->num_cusps; i++) {
            cusp = cusps[i];

            printf("Boundary %d\n", i);
            g = cusp->cusp_region_graph;
            for (j = 0; j < 4 * manifold->num_tetrahedra; j++) {
                for (region = cusp->cusp_region_begin[j].next;
                     region != &cusp->cusp_region_end[j];
                     region = region->next) {
                    if (region->index >= g->num_vertices) {
                        continue;
                    }

                    printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                           region->index, region->tet_index, region->tet_vertex
                    );
                    for (edge_node = g->edge_list_begin[region->index].next;
                         edge_node != &g->edge_list_end[region->index];
                         edge_node = edge_node->next)
                        printf("%d ", edge_node->y);

                    printf("\n");
                }
            }
        }
    } else if (strcmp(type, "endpoints") == 0) {
        printf("EndPoint Info\n");

        // which curve
        for (i = 0; i < curves->num_curves; i++) {
            printf("Dual Curve %d\n", i);

            j = 0;
            // which component
            for (path = curves->curve_begin[i].next; path != &curves->curve_end[i]; path = path->next) {
                printf("    Part %d Cusp %d\n", j, path->endpoints[0].tri->tet->cusp[path->endpoints[0].tri->tet_vertex]->index);
                for (k = 0; k < 2; k++) {
                    if (k == 0)
                        printf("        Start: ");
                    else
                        printf("        End:   ");

                    x_vertex1 = (int) remaining_face[path->endpoints[k].tri->tet_vertex][path->endpoints[k].vertex];
                    x_vertex2 = (int) remaining_face[path->endpoints[k].vertex][path->endpoints[k].tri->tet_vertex];

                    printf("Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d) Adj Curves %d\n",
                           path->endpoints[k].region_index, path->endpoints[k].tri->tet_index,
                           path->endpoints[k].tri->tet_vertex, path->endpoints[k].face, path->endpoints[k].vertex,
                           path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_class,
                           path->endpoints[k].tri->vertices[path->endpoints[k].vertex].edge_index,
                           path->endpoints[k].num_adj_curves);
                }

                j++;
            }
        }
    } else {
        printf("Unknown type: %s\n", type);
    }
    printf("-------------------------------\n");
}

// ------------------------------------

/*
 * Find oscillating curves. Each curve is made up of an even number of
 * components, with each component contained in a cusp, and connecting
 * two cusp vertices. Each oscillating curve is associated to an edge
 * of the triangulation, the rest of the edges come from the end multi
 * graph.
 *
 * The result is stored in tet->extra[edge_class].curve[f][v] array
 * on each tetrahedron.
 */

void do_oscillating_curves(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    CuspEndPoint cusp_path_begin, cusp_path_end, *temp_cusp;
    int i;

    cusp_path_begin.next = &cusp_path_end;
    cusp_path_begin.prev = NULL;
    cusp_path_end.next   = NULL;
    cusp_path_end.prev   = &cusp_path_begin;

    for (i = 0; i < curves->num_curves; i++) {
        find_multi_graph_path(cusps[0]->manifold, multi_graph,
                              &cusp_path_begin, &cusp_path_end, curves->edge_class[i]);
        do_one_oscillating_curve(cusps, curves, multi_graph, &cusp_path_begin, &cusp_path_end,
                                 curves->edge_class[i], i);

        while (cusp_path_begin.next != &cusp_path_end) {
            temp_cusp = cusp_path_begin.next;
            REMOVE_NODE(temp_cusp);
            my_free(temp_cusp);
        }

        if (debug) {
            printf("\n");
            printf("Oscillating Curve %d\n", i);
            printf("\n");
            printf("-------------------------------\n");

            log_structs(cusps[0]->manifold, cusps, curves, "dual_curves");
            log_structs(cusps[0]->manifold, cusps, curves, "endpoints");
            log_structs(cusps[0]->manifold, cusps, curves, "cusp_regions");
            log_structs(cusps[0]->manifold, cusps, curves, "graph");
        }
    }
}

/*
 * Construct a curve dual to the edge class 'edge_class'. The first and last 
 * components connect to edge_class which is not in the end multi graph so 
 * we need to find a new curve. Any intermediate components, if they exist, will
 * make use of the train lines, as they consist of curves between edge classes 
 * in the end multi graph and thus is a segment of the train line.
 */

void do_one_oscillating_curve(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph,
                              CuspEndPoint *cusp_path_begin, CuspEndPoint *cusp_path_end,
                              int edge_class, int curve_index) {
    int orientation = START;
    CuspEndPoint *endpoint = cusp_path_begin->next;
    CurveComponent *path,
                   *curve_begin = &curves->curve_begin[curve_index],
                   *curve_end = &curves->curve_end[curve_index];

    curve_begin->edge_class[FINISH] = edge_class;
    curve_end->edge_class[START]    = edge_class;

    path = setup_first_curve_component(cusps[endpoint->cusp_index], multi_graph, endpoint,
                                       curve_begin, curve_end);
    do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    update_path_holonomy(path, edge_class);

    // interior curve components (not used for knots)
    for (endpoint = endpoint->next; endpoint->next != cusp_path_end; endpoint = endpoint->next) {
        // not implemented
        uFatalError("do_one_oscillating_curve", "symplectic_basis");
    }

    orientation = (orientation == START ? FINISH : START);

    path = setup_last_curve_component(cusps[endpoint->cusp_index], multi_graph, endpoint,
                                      curve_begin, curve_end);
    do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    update_path_holonomy(path, edge_class);

    update_adj_curve_along_path(cusps, curves, curve_index,
                                (Boolean) (cusp_path_begin->next->next->next != cusp_path_end));
}

/*
 * Initalise the first curve component of an oscillating curve.
 * Set edge classes and find path endpoints.
 */

CurveComponent *setup_first_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                            CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START], endpoint->edge_class[FINISH], endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end);

    construct_cusp_region_dual_graph(cusp);
    find_single_endpoint(cusp, &path->endpoints[START], path->edge_class[START], START);
    find_single_endpoint(cusp, &path->endpoints[FINISH], path->edge_class[FINISH], START);
    return path;
}

/*
 * Initalise the last curve component of an oscillating curve.
 * Set edge classes and find path endpoints.
 */

CurveComponent *setup_last_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                           CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START], endpoint->edge_class[FINISH], endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end);

    construct_cusp_region_dual_graph(cusp);
    find_single_matching_endpoint(cusp, &curves_begin->next->endpoints[START], &path->endpoints[START]);
    find_single_matching_endpoint(cusp, &path->prev->endpoints[FINISH], &path->endpoints[FINISH]);

    return path;
}

/*
 * Construct an oscillating curve component, which is either the
 * first or last component of an oscillating curve.
 */

void do_curve_component_to_new_edge_class(CuspStructure *cusp, CurveComponent *curve) {
    int *parent;
    Boolean *processed, *discovered;
    EdgeNode node_begin, node_end;

    processed   = NEW_ARRAY(cusp->cusp_region_graph->num_vertices, Boolean);
    discovered  = NEW_ARRAY(cusp->cusp_region_graph->num_vertices, Boolean);
    parent      = NEW_ARRAY(cusp->cusp_region_graph->num_vertices, int);

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    // Find curve using bfs
    init_search(cusp->cusp_region_graph, processed, discovered, parent);
    bfs(cusp->cusp_region_graph, curve->endpoints[START].region_index, processed, discovered, parent);

    find_path(curve->endpoints[START].region_index, curve->endpoints[FINISH].region_index,
              parent, &node_begin, &node_end);
    graph_path_to_dual_curve(cusp, &node_begin, &node_end,
                             &curve->path_begin, &curve->path_end,
                             &curve->endpoints[START], &curve->endpoints[FINISH]);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the curve
    split_cusp_regions_along_path(cusp, &curve->path_begin, &curve->path_end,
                                  &curve->endpoints[START], &curve->endpoints[FINISH]);

    free_edge_node(&node_begin, &node_end);
}


/*
 * Find a cusp region which can dive along a face into a vertex of
 * the cusp triangle which corresponds to 'edge_class' and 'edge_index',
 * and store the result in path_endpoint.
 */

void find_single_endpoint(CuspStructure *cusp, PathEndPoint *path_endpoint, int edge_class, int edge_index) {
    int i;
    VertexIndex vertex;
    FaceIndex face1, face2, face;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < cusp->num_cusp_triangles; i++) {
        for (region = cusp->cusp_region_begin[i].next;
             region != &cusp->cusp_region_end[i];
             region = region->next) {

            // which vertex to dive through
            for (vertex = 0; vertex < 4; vertex++) {
                if (vertex == region->tet_vertex)
                    continue;

                if (region->tri->vertices[vertex].edge_class != edge_class)
                    continue;

                if (region->tri->vertices[vertex].edge_index != edge_index)
                    continue;

                face1 = remaining_face[region->tet_vertex][vertex];
                face2 = remaining_face[vertex][region->tet_vertex];

                if (region->dive[face1][vertex])
                    face = face1;
                else if (region->dive[face2][vertex])
                    face = face2;
                else
                    continue;

                path_endpoint->region           = region;
                path_endpoint->tri              = region->tri;
                path_endpoint->vertex           = vertex;
                path_endpoint->face             = face;
                path_endpoint->region_index     = region->index;
                path_endpoint->num_adj_curves   = region->num_adj_curves[path_endpoint->face][path_endpoint->vertex];

                return ;
            }
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_endpoints", "symplectic_basis");
}

/*
 * Find a cusp region which can dive into a vertex of the cusp triangle
 * corresponding 'edge_class' and 'edge_index', while matching path_endpoint1.
 *
 * See 'region_index', 'region_vertex', 'region_dive', 'region_curve' for the
 * conditions for a matching endpoint.
 */

void find_single_matching_endpoint(CuspStructure *cusp, PathEndPoint *path_endpoint1, PathEndPoint *path_endpoint2) {
    int i;
    Boolean region_index, region_vertex, region_dive, region_curve;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < cusp->num_cusp_triangles; i++) {
        for (region = cusp->cusp_region_begin[i].next;
             region != &cusp->cusp_region_end[i];
             region = region->next) {

            // are we in the matching endpoint
            region_index    = (Boolean) (region->tet_index != path_endpoint1->tri->tet_index);
            region_vertex   = (Boolean) (region->tet_vertex != path_endpoint1->vertex);
            region_dive     = (Boolean) !region->dive[path_endpoint1->face][path_endpoint1->tri->tet_vertex];
            region_curve    = (Boolean) (region->num_adj_curves[path_endpoint1->face][path_endpoint1->tri->tet_vertex]
                    != path_endpoint1->num_adj_curves);

            if (region_index || region_vertex || region_dive || region_curve)
                continue;

            path_endpoint2->region          = region;
            path_endpoint2->tri             = region->tri;
            path_endpoint2->vertex          = path_endpoint1->tri->tet_vertex;
            path_endpoint2->face            = path_endpoint1->face;
            path_endpoint2->region_index    = region->index;
            path_endpoint2->num_adj_curves  = region->num_adj_curves[path_endpoint2->face][path_endpoint2->vertex];

            return ;
        }
    }

    // didn't find valid path endpoints
    uFatalError("find_single_matching_endpoints", "symplectic_basis");
}

/*
 * After finding a path, each node contains the index of the region it lies in. 
 * Update path info calculates the face the path crosses to get to the next node 
 * and the vertex it cuts off to simplify combinatorial holonomy calculation.
 */

void graph_path_to_dual_curve(CuspStructure *cusp, EdgeNode *node_begin, EdgeNode *node_end, PathNode *path_begin,
                              PathNode *path_end, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    FaceIndex face;
    EdgeNode *edge_node;
    PathNode *path_node;
    CuspRegion *region;

    // path len 0
    if (node_begin->next == node_end)
        return;

    edge_node = node_begin->next;
    // path len 1
    if (edge_node->next == node_end) {
        for (face = 0; face < 4; face++)
            if (cusp->regions[edge_node->y]->tet_vertex != face &&
                start_endpoint->vertex != face &&
                finish_endpoint->vertex != face)
                break;

        region = cusp->regions[edge_node->y];

        path_node = NEW_STRUCT( PathNode );
        INSERT_BEFORE(path_node, path_end);
        path_node->next_face = finish_endpoint->face;
        path_node->prev_face = start_endpoint->face;
        path_node->cusp_region_index = edge_node->y;
        path_node->tri = region->tri;
        path_node->inside_vertex = face;
        return;
    }

    // Set Header node
    endpoint_edge_node_to_path_node(cusp->regions[edge_node->y], path_end, edge_node,
                                    start_endpoint, START);

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next)
        interior_edge_node_to_path_node(cusp->regions[edge_node->y], path_end, edge_node);

    // Set Tail node
    endpoint_edge_node_to_path_node(cusp->regions[edge_node->y], path_end, edge_node,
                                    finish_endpoint, FINISH);
}

void endpoint_edge_node_to_path_node(CuspRegion *region, PathNode *path_end, EdgeNode *edge_node,
                                     PathEndPoint *path_endpoint, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        path_node->next_face = -1;
        for (face = 0; face < 4; face++) {
            if (face == region->tet_vertex || !region->adj_cusp_triangle[face] || path_node->next_face != -1)
                continue;

        if (region->adj_cusp_regions[face]->index == edge_node->next->y)
            path_node->next_face = face;
        }

        // next node isn't in an adjacent region
        if (path_node->next_face == -1)
            uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

        path_node->prev_face = path_endpoint->face;

        if (path_node->next_face == path_endpoint->vertex) {
            if (path_endpoint->face == vertex1)
                path_node->inside_vertex = vertex2;
            else
                path_node->inside_vertex = vertex1;
        } else if (path_node->next_face == path_endpoint->face) {
            path_node->inside_vertex = -1;
        } else {
            path_node->inside_vertex = path_endpoint->vertex;
        }
    } else {
        path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face],
                                        path_end->prev->next_face);
        path_node->next_face = path_endpoint->face;

        if (path_node->prev_face == path_endpoint->vertex) {
            if (path_endpoint->face == vertex1)
                path_node->inside_vertex = vertex2;
            else
                path_node->inside_vertex = vertex1;
        } else if (path_node->prev_face == path_endpoint->face) {
            path_node->inside_vertex = -1;
        } else {
            path_node->inside_vertex = path_endpoint->vertex;
        }
    }

    INSERT_BEFORE(path_node, path_end);
}

/*
 * node lies in 'region', find the vertex which the subpath 
 * node->prev->y --> node->y --> node->next->y cuts off of the cusp triangle 
 * >tri.
 */

void interior_edge_node_to_path_node(CuspRegion *region, PathNode *path_end, EdgeNode *edge_node) {
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face],
                                    path_end->prev->next_face);

    vertex1 = remaining_face[path_node->tri->tet_vertex][path_node->prev_face];
    vertex2 = remaining_face[path_node->prev_face][path_node->tri->tet_vertex];

    if (region->adj_cusp_triangle[vertex1] && region->adj_cusp_regions[vertex1]->index == edge_node->next->y) {
        path_node->next_face = vertex1;
        path_node->inside_vertex = vertex2;
    } else if (region->adj_cusp_triangle[vertex2] && region->adj_cusp_regions[vertex2]->index == edge_node->next->y) {
        path_node->next_face = vertex2;
        path_node->inside_vertex = vertex1;
    } else
        uFatalError("interior_edge_node_to_path_node", "symplectic_basis");

    INSERT_BEFORE(path_node, path_end);
}

/*
 * The oscillating curve splits the region it passes through into two regions. 
 * Split each region in two and update attributes
 */

void split_cusp_regions_along_path(CuspStructure *cusp, PathNode *path_begin, PathNode *path_end,
                                   PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, region_index;
    PathNode *node;
    CuspRegion *region;

    // empty path
    if (path_begin->next == path_end)
        return ;

    // path of len 1
    if (path_begin->next->next == path_end) {
        split_path_len_one(cusp, path_begin->next, start_endpoint, finish_endpoint);
        return;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and region to the left of the curve.
     */
    node = path_begin->next;
    region = cusp->regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   region, start_endpoint, node, START);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], region,
                                    node, start_endpoint, index, START);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        region = cusp->regions[node->cusp_region_index];
        region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin[region_index],
                                           &cusp->cusp_region_end[region_index], region, node);
        split_cusp_region_path_interior(&cusp->cusp_region_end[region_index], region, node, index);
        index++;
    }

    // update last region
    region = cusp->regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(region->tet_index, region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   region, finish_endpoint, node, FINISH);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], region,
                                    node, finish_endpoint, index, FINISH);
    index++;

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

void split_path_len_one(CuspStructure *cusp, PathNode *node, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, region_index;
    FaceIndex face;
    CuspRegion *new_region, *old_region, *region;

    new_region = NEW_STRUCT(CuspRegion);
    old_region = cusp->regions[node->cusp_region_index];
    region_index = TRI_TO_INDEX(old_region->tet_index, old_region->tet_vertex);
    INSERT_BEFORE(new_region, &cusp->cusp_region_end[region_index])
    copy_region(old_region, new_region);

    face = node->inside_vertex;

    new_region->index = index;
    new_region->adj_cusp_triangle[start_endpoint->vertex]                   = FALSE;
    new_region->adj_cusp_triangle[finish_endpoint->vertex]                  = FALSE;
    new_region->dive[face][start_endpoint->vertex]                          = TRUE;
    new_region->dive[face][finish_endpoint->vertex]                         = TRUE;
    new_region->dive[start_endpoint->vertex][finish_endpoint->vertex]       = (Boolean) (face != finish_endpoint->face);
    new_region->dive[finish_endpoint->vertex][start_endpoint->vertex]       = (Boolean) (face != start_endpoint->face);
    new_region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
    new_region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;

    old_region->adj_cusp_triangle[face]             = FALSE;
    old_region->dive[face][start_endpoint->vertex]  = (Boolean) (face == start_endpoint->face);
    old_region->dive[face][finish_endpoint->vertex] = (Boolean) (face == finish_endpoint->face);
    old_region->temp_adj_curves[face][start_endpoint->vertex]++;
    old_region->temp_adj_curves[face][finish_endpoint->vertex]++;

    // update other cusp regions
    for (region = cusp->cusp_region_begin[region_index].next;
         region != &cusp->cusp_region_end[region_index];
         region = region->next) {

        if (new_region->tet_index != region->tet_index || new_region->tet_vertex != region->tet_vertex)
            continue;

        if (region == new_region || region == old_region)
            continue;

        if (region->adj_cusp_triangle[start_endpoint->vertex] || region->adj_cusp_triangle[finish_endpoint->vertex]) {
            region->temp_adj_curves[face][finish_endpoint->vertex]++;
            region->temp_adj_curves[face][start_endpoint->vertex]++;

        } else {
            region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
            region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;
        }
    }

    update_adj_region_data(cusp);
    cusp->num_cusp_regions++;
}

/*
 * Set the new and old region data. Draw a picture to see how the attributes 
 * change in each case
 */

void split_cusp_region_path_interior(CuspRegion *region_end, CuspRegion *region, PathNode *node, int index) {
    int v1, v2;
    CuspRegion *new_region = NEW_STRUCT( CuspRegion );

    v1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    v2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    /*
     * new_region becomes the cusp region closest to the inside vertex and
     * region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    // Update new region
    new_region->curve[v1][v2]++;
    new_region->curve[v2][v1]++;
    new_region->dive[v1][node->inside_vertex]           = region->dive[v1][node->inside_vertex];
    new_region->dive[v2][node->inside_vertex]           = region->dive[v2][node->inside_vertex];
    new_region->adj_cusp_triangle[node->inside_vertex]  = FALSE;

    // Update region
    region->curve[v1][node->inside_vertex]++;
    region->curve[v2][node->inside_vertex]++;
    region->dive[v1][node->inside_vertex]           = FALSE;
    region->dive[v2][node->inside_vertex]           = FALSE;

    INSERT_BEFORE(new_region, region_end);
}

void split_cusp_region_path_endpoint(CuspRegion *region_end, CuspRegion *region, PathNode *path_node, 
                                     PathEndPoint *path_endpoint, int index, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;
    path_endpoint->region = NULL;

    if (pos == START) {
        face = path_node->next_face;
    } else {
        face = path_node->prev_face;
    }

    if (face == path_endpoint->vertex) {
        // curve passes through the face opposite the vertex it dives through
        new_region->curve[path_endpoint->vertex][vertex2]++;
        new_region->temp_adj_curves[vertex1][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]      = (Boolean) (path_endpoint->face == vertex1);
        new_region->dive[vertex2][path_endpoint->vertex]      = region->dive[vertex2][path_endpoint->vertex];
        new_region->dive[vertex2][vertex1]                    = region->dive[vertex2][vertex1];
        new_region->dive[path_endpoint->vertex][vertex1]      = region->dive[path_endpoint->vertex][vertex1];
        new_region->adj_cusp_triangle[vertex1]                = FALSE;

        region->curve[path_endpoint->vertex][vertex1]++;
        region->temp_adj_curves[vertex2][path_endpoint->vertex]++;
        region->dive[vertex2][path_endpoint->vertex]         = (Boolean) (path_endpoint->face == vertex2);
        region->dive[vertex2][vertex1]                       = FALSE;
        region->dive[path_endpoint->vertex][vertex1]         = FALSE;
        region->adj_cusp_triangle[vertex2]                   = FALSE;
    } else if (face == path_endpoint->face) {
        // curve passes through the face that carries it
        new_region->curve[path_endpoint->face][path_endpoint->face == vertex1 ? vertex2 : vertex1]++;
        new_region->temp_adj_curves[face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex]
                        = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex]                                 = FALSE;
        new_region->adj_cusp_triangle[path_endpoint->face == vertex1 ? vertex2 : vertex1]    = FALSE;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
    } else {
        // Curve goes around the vertex
        new_region->curve[face][path_endpoint->face]++;
        new_region->temp_adj_curves[path_endpoint->face][path_endpoint->vertex]++;
        new_region->dive[vertex1][path_endpoint->vertex]              = region->dive[vertex1][path_endpoint->vertex];
        new_region->dive[vertex2][path_endpoint->vertex]              = region->dive[vertex2][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->face]            = FALSE;
        new_region->adj_cusp_triangle[path_endpoint->vertex]          = FALSE;

        region->curve[face][path_endpoint->vertex]++;
        region->temp_adj_curves[face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face == vertex1 ? vertex2 : vertex1][path_endpoint->vertex] = FALSE;
    }

    INSERT_BEFORE(new_region, region_end);
}

/*
 * After splitting each region the path travels through, the attributes for 
 * other regions in the same cusp triangle is now out of date. Update cusp 
 * triangles for nodes in the interior of the path.
 */

void update_cusp_triangle_path_interior(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
                                        CuspRegion *region, PathNode *node) {
    int face1, face2;
    CuspRegion *current_region;

    face1 = (int) remaining_face[region->tet_vertex][node->inside_vertex];
    face2 = (int) remaining_face[node->inside_vertex][region->tet_vertex];

    for (current_region = cusp_region_start->next;
         current_region != cusp_region_end;
         current_region = current_region->next) {

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (current_region->curve[face1][node->inside_vertex] > region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][node->inside_vertex]++;
        }
        else if (current_region->curve[face1][node->inside_vertex] < region->curve[face1][node->inside_vertex]) {
            current_region->curve[face1][face2]++;
        }

        if (current_region->curve[face2][node->inside_vertex] > region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][node->inside_vertex]++;
        }
        else if (current_region->curve[face2][node->inside_vertex] < region->curve[face2][node->inside_vertex]) {
            current_region->curve[face2][face1]++;
        }
    }
}

/*
 * After splitting each curveRegion the path travels through, the attributes 
 * for other regions in the same cusp triangle is now out of date. Update cusp 
 * triangles for nodes at the end of the path.
 */

void update_cusp_triangle_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end, CuspRegion *region,
                                    PathEndPoint *path_endpoint, PathNode *node, int pos) {
    FaceIndex face, face1, face2;
    CuspRegion *current_region;

    face1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    face2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    if (pos == START) {
        face = node->next_face;
    } else {
        face = node->prev_face;
    }

    for (current_region = cusp_region_start->next;
         current_region != cusp_region_end;
         current_region = current_region->next) {
        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (face == path_endpoint->vertex) {
            // curve passes through the face opposite the vertex it dives through
            if (!current_region->adj_cusp_triangle[face]) {
                if (!current_region->adj_cusp_triangle[face1]) {
                    current_region->temp_adj_curves[face1][path_endpoint->vertex]++;
                } else if (!current_region->adj_cusp_triangle[face2]) {
                    current_region->temp_adj_curves[face2][path_endpoint->vertex]++;
                } else {
                    uFatalError("update_cusp_triangle_endpoints", "symplectic_basis");
                }
            } else if (current_region->curve[path_endpoint->vertex][face1] > region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face1]++;
                current_region->temp_adj_curves[face2][path_endpoint->vertex]++;
            } else if (current_region->curve[path_endpoint->vertex][face1] < region->curve[path_endpoint->vertex][face1]) {
                current_region->curve[face][face2]++;
                current_region->temp_adj_curves[face1][path_endpoint->vertex]++;
            }

            continue;
        }

        if (!current_region->adj_cusp_triangle[face]) {
            current_region->temp_adj_curves[face][path_endpoint->vertex]++;
            continue;
        }

        // Curve goes around the vertex or passes through the face that carries it
        if (current_region->curve[face][path_endpoint->vertex] > region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][path_endpoint->vertex]++;
            current_region->temp_adj_curves[face][path_endpoint->vertex]++;

        } else if (current_region->curve[face][path_endpoint->vertex] < region->curve[face][path_endpoint->vertex]) {
            current_region->curve[face][face == face1 ? face2 : face1]++;
            current_region->temp_adj_curves[face == face1 ? face2 : face1][path_endpoint->vertex]++;
        }
    }
}

void update_adj_curve_along_path(CuspStructure **cusps, OscillatingCurves *curves, int curve_index, Boolean train_line) {
    CurveComponent *curve,
                   *current_begin = &curves->curve_begin[curve_index],
                   *current_end   = &curves->curve_end[curve_index];

    // Update regions curve data along the current curve
    for (curve = current_begin->next; curve != current_end; curve = curve->next)
        update_adj_curve_on_cusp(cusps[curve->cusp_index]);

    // update endpoint curve data
    for (int i = 0; i < curve_index; i++) {
        // which oscillating curve

        for (curve = curves->curve_begin[i].next; curve != &curves->curve_end[i]; curve = curve->next) {
            // which component of the curve

            for (int j = 0; j < 2; j++) {
                // which end point
                update_adj_curve_at_endpoint(&curve->endpoints[j], current_begin->next, START);
                update_adj_curve_at_endpoint(&curve->endpoints[j], current_end->prev, FINISH);
            }
        }
    }
}

/*
 * curve_begin and curve_end are header and tailer nodes of a doubly linked list of path
 * components for a new path. Update the path_endpoint->num_adj_curves attribute to account for this
 * new curve.
 */

void update_adj_curve_at_endpoint(PathEndPoint *path_endpoint, CurveComponent *path, int pos) {
    PathEndPoint *curve_end_point;

    curve_end_point = &path->endpoints[pos];

    if (curve_end_point->tri->tet_index != path_endpoint->tri->tet_index ||
        curve_end_point->tri->tet_vertex != path_endpoint->tri->tet_vertex ||
        curve_end_point->face != path_endpoint->face ||
        curve_end_point->vertex != path_endpoint->vertex)
        return;

    path_endpoint->num_adj_curves++;
}

/*
 * Move the temp adj curves into the current num of adj curves.
 */

void update_adj_curve_on_cusp(CuspStructure *cusp) {
    int i, j, k;
    CuspRegion *region;

    for (i = 0; i < 4 * cusp->manifold->num_tetrahedra; i++) {
        for (region = cusp->cusp_region_begin[i].next; region != &cusp->cusp_region_end[i]; region = region->next) {
            // which cusp region
            for (j = 0; j < 4; j++) {
                for (k = 0; k < 4; k++) {
                    region->num_adj_curves[j][k] += region->temp_adj_curves[j][k];
                    region->temp_adj_curves[j][k] = 0;
                }
            }
        }
    }
}

void update_path_holonomy(CurveComponent *path, int edge_class) {
    PathNode *path_node;

    for (path_node = path->path_begin.next; path_node != &path->path_end; path_node = path_node->next) {
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->next_face]++;
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->prev_face]--;
    }
}

// -------------------------------------------------------

/*
 * End Multi Graph
 *
 * The end multi graph is a graph with vertices for the cusps of M and
 * edges for each edge of the triangulation. We also refer to the spanning
 * tree of this graph as the end multi graph. The end multi graph structure also
 * keeps track of a special E_0 edge, which is used to construct a path of even
 * length through the graph.
 */

EndMultiGraph *init_end_multi_graph(Triangulation *manifold) {
    int i, j;
    int *parent;
    EndMultiGraph *multi_graph = NEW_STRUCT( EndMultiGraph );
    
    multi_graph->num_cusps = manifold->num_cusps;
    multi_graph->num_edge_classes = manifold->num_tetrahedra;

    Graph *g = init_graph(multi_graph->num_cusps, FALSE);
    cusp_graph(manifold, g);

    parent = NEW_ARRAY(g->num_vertices, int);

    multi_graph->multi_graph = spanning_tree(g, 0, parent);
    color_graph(multi_graph->multi_graph);

    multi_graph->edges = find_end_multi_graph_edge_classes(multi_graph, manifold);
    multi_graph->e0 = find_same_color_edge(manifold, multi_graph, g);

    multi_graph->edge_classes = NEW_ARRAY(multi_graph->num_edge_classes, Boolean);
    for (i = 0; i < multi_graph->num_edge_classes; i++) {
        multi_graph->edge_classes[i] = FALSE;
    }

    for (i = 0; i < multi_graph->num_cusps; i++) { 
        for (j = 0; j < multi_graph->num_cusps; j++) { 
            if (multi_graph->edges[i][j] == -1)
                continue;

            multi_graph->edge_classes[multi_graph->edges[i][j]] = TRUE;
        }
    }

    free_graph(g);
    my_free(parent);
    return multi_graph;
}

void free_end_multi_graph(EndMultiGraph *multi_graph) {
    int i;

    free_graph(multi_graph->multi_graph);
    
    for (i = 0; i < multi_graph->num_cusps; i++)
        my_free(multi_graph->edges[i]);

    my_free(multi_graph->edge_classes);
    my_free(multi_graph->edges);
    my_free(multi_graph);
}

void cusp_graph(Triangulation *manifold, Graph *g) {
    int vertex1, vertex2;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (vertex1 = 0; vertex1 < 4; vertex1++) {
            // which vertex of the cusp triangle at vertex1
            for (vertex2 = 0; vertex2 < 4; vertex2++) {
                if (vertex1 == vertex2)
                    continue;

                insert_edge(g, tet->cusp[vertex1]->index, tet->cusp[vertex2]->index, g->directed);
            }
        }
    }
}

/*
 * Find a spanning tree of graph1
 */

Graph *spanning_tree(Graph *graph1, int start, int *parent) {
    int i;

    Boolean *processed = NEW_ARRAY(graph1->num_vertices, Boolean);
    Boolean *discovered = NEW_ARRAY(graph1->num_vertices, Boolean);

    Graph *graph2 = init_graph(graph1->num_vertices, graph1->directed);

    // Find path using bfs
    init_search(graph1, processed, discovered, parent);
    bfs(graph1, start, processed, discovered, parent);

    for (i = 0; i < graph1->num_vertices; i++) {
        if (parent[i] == -1)
            continue;

        insert_edge(graph2, i, parent[i], graph2->directed);
    }

    my_free(processed);
    my_free(discovered);

    return graph2;
}

/*
 * Assign an edge class to each edge of the graph g and return an array of 
 * Booleans indicating if an edge class is in the graph.
 */

int **find_end_multi_graph_edge_classes(EndMultiGraph *multi_graph, Triangulation *manifold) {
    int i, j, edge_class, **cusps;
    EdgeNode *edge_node;
    Graph *g = multi_graph->multi_graph;

    cusps = NEW_ARRAY(multi_graph->num_cusps, int *);

    for (i = 0; i < multi_graph->num_cusps; i++) {
        cusps[i] = NEW_ARRAY(multi_graph->num_cusps, int);

        for (j = 0; j < multi_graph->num_cusps; j++)
            cusps[i][j] = -1;
    }

    for (i = 0; i < g->num_vertices; i++) {
        for (edge_node = g->edge_list_begin[i].next; edge_node != &g->edge_list_end[i]; edge_node = edge_node->next) {
            edge_class = find_edge_class(manifold, i, edge_node->y);
            cusps[i][edge_node->y] = edge_class; 
            cusps[edge_node->y][i] = edge_class; 
        }
    }
    
    return cusps;
}

/*
 * Find an edge class whose edge connects cusp1 and cusp2
 */

int find_edge_class(Triangulation *manifold, int cusp1, int cusp2) {
    int v1, v2;
    EdgeClass *edge;
    Tetrahedron *tet;

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        for (v1 = 0; v1 < 4; v1++) {
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                if (tet->cusp[v1]->index != cusp1 || tet->cusp[v2]->index != cusp2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                return edge->index;
            }
        }
    }

    uFatalError("find_edge_class", "symplectic_basis");
    return 0;
}

void color_graph(Graph *g) {
    int color = 0, v;
    Queue *q = init_queue(g->num_vertices);
    EdgeNode *node;

    g->color[0] = color;
    q = enqueue(q, 0);

    while (!empty_queue(q)) {
        v = dequeue(q);
        color = g->color[v];

        for (node = g->edge_list_begin[v].next; node != &g->edge_list_end[v]; node = node->next) {
            // graph is not bipartite
            if (g->color[node->y] == color)
                uFatalError("color_graph", "symplectic_basis");

            if (g->color[node->y] != -1)
                continue;

            g->color[node->y] = !color;
            q = enqueue(q, node->y);
        }
    }

    free_queue(q);
}

/*
 * g1 is the colored spanning tree of g2, return the edge class of the edge in 
 * g2 which connects vertices in g1 of the same color
 */

int find_same_color_edge(Triangulation *manifold, EndMultiGraph *multi_graph, Graph *g2) {
    int cusp;
    EdgeNode *node;
    Graph *g1 = multi_graph->multi_graph;

    for (cusp = 0; cusp < g2->num_vertices; cusp++) {
        for (node = g2->edge_list_begin[cusp].next; node != &g2->edge_list_end[cusp]; node = node->next) {
            if (g1->color[cusp] == g1->color[node->y] && multi_graph->edges[cusp][node->y] == -1) 
                // we found an edge
                return find_edge_class(manifold, cusp, node->y);
        }
    }

    // we didn't find an edge connecting vertices of the same color
    uFatalError("find_same_color_edge", "symplectic_basis");
    return -1;
}

/*
 * Find the length of a path between start and end
 */

int find_path_len(int start, int end, int *parents, int path_length) {
    if ((start == end) || (end == -1)) {
        return path_length;
    } else {
        return find_path_len(start, parents[end], parents, path_length + 1);
    }
}

/*
 * Find a path through the end multi graph, starting at 'edge_class' and which
 * an odd length, since this corresponds to an even number of oscillating curve
 * components. The path is stored in the doubly linked list cusp_path_begin ->
 * cusp_path_end.
 */

void find_multi_graph_path(Triangulation *manifold, EndMultiGraph *multi_graph, CuspEndPoint *cusp_path_begin,
                           CuspEndPoint *cusp_path_end, int edge_class) {
    Graph *g = multi_graph->multi_graph;
    Boolean *processed     = NEW_ARRAY(g->num_vertices, Boolean);
    Boolean *discovered    = NEW_ARRAY(g->num_vertices, Boolean);
    int *parent         = NEW_ARRAY(g->num_vertices, int);
    int start, end, startE0, endE0, path_len = 0;
    EdgeNode node_begin, node_end;

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    find_edge_ends(g, manifold, edge_class, &start, &end);
    find_edge_ends(g, manifold, multi_graph->e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    path_len = find_path_len(start, end, parent, path_len);

    if (path_len % 2 == 1) {
        find_path(start, end, parent, &node_begin, &node_end);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, &node_begin, &node_end);

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end.prev, &node_end);
    }

    graph_path_to_cusp_path(multi_graph, &node_begin, &node_end, cusp_path_begin, cusp_path_end, edge_class);

    free_edge_node(&node_begin, &node_end);
    my_free(parent);
    my_free(discovered);
    my_free(processed);
}

/*
 * Converts the EdgeNode path through the cusps in the end multigraph to 
 * a CuspEndPoint path which contains the edge classes on each cusp. 
 * A CuspEndPoint corresponding to one section of an oscillating curve, and
 * constructing such a section for all CuspEndPoints gives the whole curve.
 *
 * node_begin -> node_end is a doubly linked list through the end multi graph
 * and the result path is stored in the doubly linked list cusp_path_begin ->
 * cusp_path_end.
 */

void graph_path_to_cusp_path(EndMultiGraph *multi_graph, EdgeNode *node_begin, EdgeNode *node_end,
                             CuspEndPoint *cusp_path_begin, CuspEndPoint *cusp_path_end, int edge_class) {
    int cusp, prev_edge_class;
    EdgeNode *node;
    CuspEndPoint *endpoint;

    prev_edge_class = edge_class;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        endpoint = NEW_STRUCT( CuspEndPoint );
        INSERT_BEFORE(endpoint, cusp_path_end);

        endpoint->cusp_index = cusp;
        endpoint->edge_class[START] = prev_edge_class;
        endpoint->edge_class[FINISH] = multi_graph->edges[node->y][node->next->y];

        if (endpoint->edge_class[FINISH] == -1)
            endpoint->edge_class[FINISH] = multi_graph->e0;
        
        prev_edge_class = endpoint->edge_class[FINISH];
    }

    endpoint = NEW_STRUCT( CuspEndPoint );
    INSERT_BEFORE(endpoint, cusp_path_end);

    endpoint->cusp_index = node->y;
    endpoint->edge_class[START] = prev_edge_class;
    endpoint->edge_class[FINISH] = edge_class;
}

void find_edge_ends(Graph *g, Triangulation *manifold, int edge_class, int *start, int *end) {
    int v1, v2;
    Tetrahedron *tet;
    EdgeClass *edge;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        // which vertex
        for (v1 = 0; v1 < 4; v1++) {
            // which vertex of the cusp triangle at v1
            for (v2 = 0; v2 < 4; v2++) {
                if (v1 == v2)
                    continue;

                edge = tet->edge_class[edge_between_vertices[v1][v2]];
                if (edge->index != edge_class)
                    continue;

                *start = tet->cusp[v1]->index;
                *end   = tet->cusp[v2]->index;
                return;
            }
        }
    }


    // didn't find the edge class in the graph
    uFatalError("find_edge_ends", "symplectic_basis");
}
