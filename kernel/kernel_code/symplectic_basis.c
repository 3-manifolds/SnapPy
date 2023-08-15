/**
 *  Symplectic Basis
 *
 *  Computes a symplectic basis of a triangulated knot or link exterior with
 *  orientable torus cusps. This symplectic matrix extends the Neumann-Zagier
 *  matrix to one which is symplectic up to factors of 2, and which arises
 *  from the triangulation of the manifold.
 *
 *  See - https://arxiv.org/abs/2208.06969
 *
 */

#include <stdio.h>
#include <string.h>
#include "SnapPea.h"
#include "kernel.h"
#include "addl_code.h"

#define at_least_two(a, b, c)                    ((a) && (b)) || ((a) && (c)) || ((b) && (c))
#define tri_to_index(tet_index, tet_vertex)      (4 * (tet_index) + (tet_vertex))

enum pos {
    START,
    FINISH
};

static int debug;

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
    struct CuspRegion           **regions;               /** list of regions in the graph */
    int                         *degree;                 /** degree of each vertex */
    int                         *color;                  /** color a tree bipartite */
    int                         num_vertices;            /** number of vertices in the graph */
    Boolean                     directed;                /** is the graph directed */
} Graph;

typedef struct CuspEndPoint {
    int                         cusp_index;
    int                         edge_class[2];
    struct CuspEndPoint         *next;
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
    PathNode                    curves_begin;           /** header node of doubbly linked list */
    PathNode                    curves_end;             /** tailer node of ... */
    PathEndPoint                endpoints[2];           /** path end points */
    struct CurveComponent       *next;                  /** next curve component in doubly linked list */
    struct CurveComponent       *prev;                  /** prev ... */
} CurveComponent;

typedef struct OscillatingCurves {
    int                         num_curves;
    int                         *edge_class;
    CurveComponent              *curve_begin;          /** array of doubly linked lists of dual curves */
    CurveComponent              *curve_end;            /** array of doubly linkek lists of dual curves */
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
    Graph                       *dual_graph;            /** dual graph of the cusp region */
    CuspTriangle                cusp_triangle_begin;    /** header node of doubly linked list of cusp triangles */
    CuspTriangle                cusp_triangle_end;      /** tailer node of ... */
    CuspRegion                  *cusp_region_begin;     /** array of header nodes for cusp regions, index by cusp tri */
    CuspRegion                  *cusp_region_end;       /** array of tailer nodes for ...*/
    PathNode                    train_line_path_begin;  /** header node of doubly linked list of train line node */
    PathNode                    train_line_path_end;    /** tailer node of ... */
    PathEndPoint                *train_line_endpoint;   /** train line endpoints indexed by edge class */
    PathEndPoint                extra_endpoint_e0;      /** extra endpoint in case E_0 lies in one cusp */
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

/**
 * Symplectic Basis
 */

int                     *gluing_equations_for_edge_class(Triangulation *, int);
int                     **get_symplectic_equations(Triangulation *, Boolean *, int *, int);

CuspStructure           *init_cusp_structure(Triangulation *, Cusp *);
void                    free_cusp_structure(CuspStructure **, int, int);
void                    init_cusp_triangulation(Triangulation *, CuspStructure *);
int                     net_flow_around_vertex(CuspTriangle *, int);
void                    label_triangulation_edges(Triangulation *);
CuspTriangle            *find_cusp_triangle(CuspTriangle *, CuspTriangle *, CuspTriangle *, int);
void                    label_cusp_vertex_indices(CuspTriangle *, CuspTriangle *, int);
void                    walk_around_cusp_vertex(CuspTriangle *, int, int);
void                    init_cusp_region(CuspStructure *);
int                     init_intersect_cusp_region(CuspStructure *, CuspTriangle *, int);
int                     init_intersect_vertex_two_zero_flows(CuspStructure *, CuspTriangle *, int);
int                     init_normal_cusp_region(CuspStructure *, CuspTriangle *, int);
void                    set_cusp_region_data(CuspStructure *, CuspTriangle *, const int [4], const Boolean [4], int);
void                    update_adj_region_data(CuspStructure *);
CuspRegion              *find_adj_region(CuspRegion *, CuspRegion *, CuspRegion *, int);
void                    init_train_line(CuspStructure *);
CurveComponent          *init_curve_component(int, int, int);
OscillatingCurves       *init_oscillating_curves(Triangulation *, const Boolean *);
void                    free_oscillating_curves(OscillatingCurves *);
void                    find_intersection_triangle(Triangulation *, CuspStructure *);
void                    construct_cusp_region_dual_graph(CuspStructure *);
void                    log_structs(Triangulation *, CuspStructure **, OscillatingCurves *, char *);

/**
 * Train lines
 */

void                    do_manifold_train_lines(CuspStructure **, EndMultiGraph *);
int                     *find_tet_index_for_edge_classes(Triangulation *, EndMultiGraph *);
void                    find_edge_class_edges(CuspStructure **, EndMultiGraph *);
void                    find_edge_class_edges_on_cusp(CuspStructure *, EndMultiGraph *, const Boolean *, const int *);
void                    find_e0_edges_on_cusp(CuspStructure **, EndMultiGraph *, const int *);
Boolean                 *update_edge_classes_on_cusp(CuspStructure **, Boolean **, int, int, int);
void                    find_primary_train_line(CuspStructure **, CuspStructure *, EndMultiGraph *);
void                    do_initial_train_line_segment_on_cusp(CuspStructure *, PathEndPoint *, PathEndPoint *);
void                    do_train_line_segment_on_cusp(CuspStructure *, PathEndPoint *, PathEndPoint *);
int                     next_valid_endpoint_index(PathEndPoint *, int, int);
void                    tri_endpoint_to_region_endpoint(CuspStructure *, PathEndPoint *);
void                    graph_path_to_path_node(Graph *, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_regions_along_train_line_segment(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    update_cusp_triangle_train_line_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int);
void                    split_cusp_region_train_line_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    cycle_path(Graph *, EdgeNode *, EdgeNode *, Boolean *, Boolean *, int *, int, int, int, int, int);

/**
 * Construct Oscillating Curves and calculate holonomy
 */

void                    do_oscillating_curves(CuspStructure **, OscillatingCurves *, EndMultiGraph *);
void                    do_one_oscillating_curve(CuspStructure **, OscillatingCurves *, EndMultiGraph *, int, int);
void                    copy_path_endpoint(PathEndPoint *, PathEndPoint *);
PathNode                *copy_path_node(PathNode *);
CurveComponent          *setup_train_line_component(CuspStructure *, EndMultiGraph *, CurveComponent *, CurveComponent *, CuspEndPoint *, int);
void                    do_curve_component_on_train_line(CuspStructure *, CurveComponent *);
CurveComponent          *setup_first_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
CurveComponent          *setup_last_curve_component(CuspStructure *, EndMultiGraph *, CuspEndPoint *, CurveComponent *, CurveComponent *);
void                    do_curve_component_to_new_edge_class(CuspStructure *, CurveComponent *);
void                    find_single_endpoint(Graph *, PathEndPoint *, int, int);
void                    find_single_matching_endpoint(Graph *, PathEndPoint *, PathEndPoint *, int, int);
void                    find_train_line_endpoint(CuspStructure *, PathEndPoint *, int, int, int, Boolean);
void                    graph_path_to_dual_curve(Graph *g, EdgeNode *, EdgeNode *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    endpoint_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *, PathEndPoint *, int);
void                    interior_edge_node_to_path_node(CuspRegion *, PathNode *, EdgeNode *);
void                    split_cusp_regions_along_path(CuspStructure *, PathNode *, PathNode *, PathEndPoint *, PathEndPoint *);
void                    split_cusp_region_path_interior(CuspRegion *, CuspRegion *, PathNode *, int);
void                    split_cusp_region_path_endpoint(CuspRegion *, CuspRegion *, PathNode *, PathEndPoint *, int, int);
void                    update_cusp_triangle_path_interior(CuspRegion *, CuspRegion *, CuspRegion *, PathNode *);
void                    update_cusp_triangle_endpoints(CuspRegion *, CuspRegion *, CuspRegion *, PathEndPoint *, PathNode *, int);
void                    copy_region(CuspRegion *, CuspRegion *);
void                    update_adj_curve_along_path(CuspStructure **, OscillatingCurves *, int);
void                    update_adj_curve_at_endpoint(PathEndPoint *, CurveComponent *, int);
void                    update_adj_curve_on_cusp(CuspStructure *);
void                    update_path_holonomy(CurveComponent *, int);
void                    calculate_holonomy(Triangulation *, int **, int);

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
CuspEndPoint            *find_multi_graph_path(EndMultiGraph *, Triangulation *, int, int *);
CuspEndPoint            *graph_path_to_cusp_path(EndMultiGraph *, EdgeNode *, EdgeNode *, int);
void                    find_edge_ends(Graph *, Triangulation *, int, int *, int *);
Boolean                 *edge_classes_on_cusp(EndMultiGraph *, int);

int edgesThreeToFour[4][3] = {{1, 2, 3},
                              {0, 2, 3},
                              {0, 1, 3},
                              {0, 1, 2}};

// -------------------------------------------------

// Data Structures

// Queue Data Structure

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

// Graph Data Structure

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
    g->regions              = NEW_ARRAY(max_vertices, CuspRegion *);

    for (i = 0; i < max_vertices; i++) {
        g->degree[i] = 0;
        g->color[i] = -1;

        g->edge_list_begin[i].next     = &g->edge_list_end[i];
        g->edge_list_begin[i].prev     = NULL;
        g->edge_list_end[i].next       = NULL;
        g->edge_list_end[i].prev       = &g->edge_list_begin[i];
        g->regions[i]                  = NULL;
    }

    return g;
}

void free_graph(Graph *g) {
    int i;
    EdgeNode *edge;

    if (g == NULL)
        return;

    for (i = 0; i < g->num_vertices; i++) {
        while (g->edge_list_begin[i].next != &g->edge_list_end[i]) {
            edge = g->edge_list_begin[i].next;
            REMOVE_NODE(edge);
            my_free(edge);
        }
    }

    my_free(g->edge_list_begin);
    my_free(g->edge_list_end);
    my_free(g->degree);
    my_free(g->color);
    my_free(g->regions);
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

/*
 * Find a cycle
 */

Boolean cycle_exists(Graph *g, int start, Boolean *processed, Boolean *discovered,
                     int *parent, int *cycle_start, int *cycle_end) {
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

            if (processed[y] && y != parent[v]) {
                free_queue(q);
                *cycle_start = y;
                *cycle_end = v;
                return TRUE;
            }

            if (!discovered[y]) {
                q = enqueue(q, y);
                discovered[y] = TRUE;
                parent[y] = v;
            }
        }
    }

    free_queue(q);
    *cycle_start = 0;
    *cycle_end = 0;
    return FALSE;
}

int **ford_fulkerson(Graph *g, int source, int sink) {
    int i, j, augment;
    int **residual_network = NEW_ARRAY(g->num_vertices, int *);
    EdgeNode *node;
    Boolean *visited = NEW_ARRAY(g->num_vertices, Boolean);

    for (i = 0; i < g->num_vertices; i++) {
        residual_network[i] = NEW_ARRAY(g->num_vertices, int);

        for (j = 0; j < g->num_vertices; j++)
            residual_network[i][j] = -1;
    }


    for (i = 0; i < g->num_vertices; i++) {
        for (node = g->edge_list_begin[i].next; node != &g->edge_list_end[i]; node = node->next) {
            residual_network[i][node->y] = 1;
            residual_network[node->y][i] = 0;
        }
    }

    augment = 1;
    while (augment > 0) {
        for (i = 0; i < g->num_vertices; i++)
            visited[i] = FALSE;

        augment = augment_path(g, residual_network, visited, source, sink, INT_MAX);
    }

    my_free(visited);
    return residual_network;
}

int augment_path(Graph *g, int **residual_network, Boolean *visited, int u, int t, int bottleneck) {
    int residual, augment, v;
    EdgeNode *node;

    if (u == t)
        return bottleneck;

    visited[u] = TRUE;

    for (node = g->edge_list_begin[u].next; node != &g->edge_list_end[u]; node = node->next) {
        v = node->y;
        residual = residual_network[u][v];

        if (residual > 0 && !visited[v]) {
            augment = augment_path(g, residual_network, visited, v, t, MIN(bottleneck, residual));

            if (augment > 0) {
                residual_network[u][v] -= augment;
                residual_network[v][u] += augment;
                return augment;
            }
        }
    }

    return 0;
}

// ---------------------------------------------------

// Symplectic Basis

/*
 * Allocates arrays for symplectic basis and gluing equations
 * Calls the gluing_equations_for_edge_class and get_symplectic_equations 
 * functions. 
 * Constructs return array
 */

int** get_symplectic_basis(Triangulation *manifold, int *num_rows, int *num_cols, int log) {
    int i, j;
    Boolean *edge_classes = NEW_ARRAY(manifold->num_tetrahedra, Boolean);
    peripheral_curves(manifold);

    // Dual Edge Curves Gamma_i -> symplectic equations
    int **symp_eqns, symp_num_rows;
    debug = log;

    // Get Symplectic Equations
    symp_eqns = get_symplectic_equations(manifold, edge_classes, &symp_num_rows, *num_cols);

    // Construct return array
    *num_rows = 2 * (manifold->num_tetrahedra - manifold->num_cusps);
    int **eqns = NEW_ARRAY(*num_rows, int *);

    j = 0;
    for (i = 0; i < manifold->num_tetrahedra; i++) {
        if (!edge_classes[i]) {
            my_free(symp_eqns[i]);
            continue;
        }

        eqns[2 * j] = gluing_equations_for_edge_class(manifold, i);
        eqns[2 * j + 1] = symp_eqns[i];
        j++;
    }

    my_free(edge_classes);
    my_free(symp_eqns);

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
 * Initialise cusp structure on each cusp.
 * Construct train lines. 
 * Construct oscillating curves. 
 * Construct Neumann-Zagier matrix.
 */

int **get_symplectic_equations(Triangulation *manifold, Boolean *edge_classes, int *num_rows, int num_cols) {
    int i, j, k;
    label_triangulation_edges(manifold);

    CuspStructure **cusps                   = NEW_ARRAY(manifold->num_cusps, CuspStructure *);
    EndMultiGraph *end_multi_graph          = init_end_multi_graph(manifold);
    Cusp *cusp;

    for (i = 0; i < end_multi_graph->num_edge_classes; i++)
        edge_classes[i] = end_multi_graph->edge_classes[i] == TRUE ? FALSE : TRUE;

    edge_classes[end_multi_graph->e0] = FALSE;

    OscillatingCurves *dual_curves   = init_oscillating_curves(manifold, edge_classes);

    for (i = 0; i < manifold->num_cusps; i++) {
        for (cusp = manifold->cusp_list_begin.next; cusp != &manifold->cusp_list_end && cusp->index != i; cusp = cusp->next);

        if (cusp == &manifold->cusp_list_end)
            uFatalError("get_symplectic_equations", "symplectic_basis");

        cusps[i] = init_cusp_structure(manifold, cusp);
    }

    Tetrahedron *tet;
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        if (tet->extra != NULL)
            uFatalError("get_symplectic_equations", "symplectic_basis");

        tet->extra = NEW_ARRAY(manifold->num_tetrahedra, Extra);

        for (i = 0; i < manifold->num_tetrahedra; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    tet->extra[i].curve[j][k] = 0;
    }

    if (debug) {
        printf("\n");
        printf("Struct Initialisation\n");
        printf("\n");

        log_structs(manifold, cusps, dual_curves, "gluing");       
        log_structs(manifold, cusps, dual_curves, "homology");    
        log_structs(manifold, cusps, dual_curves, "edge_indices");     
        log_structs(manifold, cusps, dual_curves, "inside_edge");       
        log_structs(manifold, cusps, dual_curves, "cusp_regions");       
    }

    // Allocate Symplectic Equations Array
    *num_rows = dual_curves->num_curves;
    int **symp_eqns = NEW_ARRAY(manifold->num_tetrahedra, int *);

    for (i = 0; i < manifold->num_tetrahedra; i ++) {
        symp_eqns[i] = NEW_ARRAY(3 * manifold->num_tetrahedra, int);

        for (j = 0; j < 3 * manifold->num_tetrahedra; j++)
            symp_eqns[i][j] = 0;
    }

    do_manifold_train_lines(cusps, end_multi_graph);
    do_oscillating_curves(cusps, dual_curves, end_multi_graph);
    calculate_holonomy(manifold, symp_eqns, manifold->num_tetrahedra);

    if (debug) {
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("%d, ", cusps[i]->num_cusp_regions);
        }
        printf("\n");
    }

    free_end_multi_graph(end_multi_graph);
    free_oscillating_curves(dual_curves);
    free_cusp_structure(cusps, manifold->num_cusps, manifold->num_tetrahedra);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        my_free(tet->extra);
        tet->extra = NULL;
    }

    return symp_eqns;
}

void free_symplectic_basis(int **eqns, int num_rows) {
    int i;

    for (i = 0; i < num_rows; i++)
        my_free(eqns[i]);
    my_free(eqns);
}

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

    find_intersection_triangle(manifold, boundary);
    init_cusp_triangulation(manifold, boundary);
    init_cusp_region(boundary);
    init_train_line(boundary);

    boundary->dual_graph = NULL;
    construct_cusp_region_dual_graph(boundary);

    return boundary;
}

void free_cusp_structure(CuspStructure **cusps, int num_cusps, int num_edge_classes) {
    int cusp_index;
    CuspTriangle *tri;
    CuspRegion *region;
    PathNode *path_node;
    CuspStructure *cusp;

    for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
        cusp = cusps[cusp_index];
        // free graph
        free_graph(cusp->dual_graph);

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

        // free train line path
        while (cusp->train_line_path_begin.next != &cusp->train_line_path_end) {
            path_node = cusp->train_line_path_begin.next;
            REMOVE_NODE(path_node);
            my_free(path_node);
        }

        my_free(cusp->train_line_endpoint);
        my_free(cusp);
    }

    my_free(cusps);
}

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
                                              net_flow_around_vertex(tri, v2) + distance[v2]) + net_flow_around_vertex(tri, vertex);
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
    if (at_least_two(!net_flow_around_vertex(tri, v1), !net_flow_around_vertex(tri, v2), !net_flow_around_vertex(tri, v3)))
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
    INSERT_BEFORE(region, &cusp->cusp_region_end[tri_to_index(tri->tet_index, tri->tet_vertex)]);

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
                adj_index = tri_to_index(adj_triangle->tet_index, adj_triangle->tet_vertex);
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

void init_train_line(CuspStructure *cusp) {
    int edge_class;

    cusp->train_line_path_begin.next    = &cusp->train_line_path_end;
    cusp->train_line_path_begin.prev    = NULL;
    cusp->train_line_path_end.next      = NULL;
    cusp->train_line_path_end.prev      = &cusp->train_line_path_begin;

    cusp->train_line_endpoint = NEW_ARRAY(cusp->manifold->num_tetrahedra, PathEndPoint);

    for (edge_class = 0; edge_class < cusp->manifold->num_tetrahedra; edge_class++) {
        cusp->train_line_endpoint[edge_class].tri               = NULL;
        cusp->train_line_endpoint[edge_class].region            = NULL;
        cusp->train_line_endpoint[edge_class].num_adj_curves    = 0;
    }

    cusp->extra_endpoint_e0.tri                 = NULL;
    cusp->extra_endpoint_e0.region              = NULL;
    cusp->extra_endpoint_e0.node                = NULL;
    cusp->extra_endpoint_e0.num_adj_curves      = 0;
}

CurveComponent *init_curve_component(int edge_class_start, int edge_class_finish, int cusp_index) {
    int i;

    CurveComponent *path = NEW_STRUCT(CurveComponent );

    path->curves_begin.next     = &path->curves_end;
    path->curves_begin.prev     = NULL;
    path->curves_end.next       = NULL;
    path->curves_end.prev       = &path->curves_begin;

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

            while (path->curves_begin.next != &path->curves_end) {
                path_node = path->curves_begin.next;
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

/*
 * Construct the graph dual to the cusp regions, using region->index to label 
 * each vertex, and adding edges using region->adj_cusp_regions[].
 */

void construct_cusp_region_dual_graph(CuspStructure *cusp) {
    int i, face;
    CuspRegion *region;

    Graph *graph1 = init_graph(cusp->num_cusp_regions, FALSE);

    int *visited = NEW_ARRAY(graph1->num_vertices, int);

    for (i = 0; i < graph1->num_vertices; i++)
        visited[i] = FALSE;

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
                graph1->regions[region->index] = region;
            }

            visited[region->index] = 1;
        }
    }

    free_graph(cusp->dual_graph);
    my_free(visited);

    cusp->dual_graph = graph1;
}

/*
 * Types: gluing, train_lines, cusp_regions, homology, edge_indices, 
 * dual_curves, inside_edge, graph, endpoints
 */

void log_structs(Triangulation *manifold, CuspStructure **cusps, OscillatingCurves *curves, char *type) {
    int i, j, k, x_vertex1, x_vertex2, y_vertex1, y_vertex2, v1, v2, v3;

    CuspTriangle *tri;
    CuspRegion *region;
    EdgeNode *edge_node;
    PathNode *path_node;
    CurveComponent *path;
    Graph *g;
    CuspStructure *cusp;
    PathEndPoint *endpoint;

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
    } else if (strcmp(type, "train_lines") == 0) {
        printf("Train Lines\n");
        for (i = 0; i < manifold->num_cusps; i++) {
            printf("Boundary %d\n", i);

            cusp = cusps[i];
            printf("    Train Line Path: \n");

            for (path_node = cusp->train_line_path_begin.next; path_node != &cusp->train_line_path_end; path_node = path_node->next) {
                printf("        Node %d: (Tet Index %d, Tet Vertex %d) Next Face: %d, Prev Face: %d, Inside Vertex: %d\n", 
                     path_node->cusp_region_index, path_node->tri->tet_index, path_node->tri->tet_vertex, 
                     path_node->next_face, path_node->prev_face, path_node->inside_vertex
                     );
            }

            printf("    Train Line Endpoints\n");
            for (j = 0; j < cusp->num_edge_classes; j++) {
                if (cusp->train_line_endpoint[j].tri == NULL)
                    continue;

                endpoint = &cusp->train_line_endpoint[j];
                printf("        Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d)\n",
                       endpoint->region_index, endpoint->tri->tet_index,
                       endpoint->tri->tet_vertex, endpoint->face, endpoint->vertex,
                       endpoint->tri->vertices[endpoint->vertex].edge_class,
                       endpoint->tri->vertices[endpoint->vertex].edge_index);
            }

            if (cusp->extra_endpoint_e0.tri != NULL) {
                endpoint = &cusp->extra_endpoint_e0;

                printf("        Region %d (Tet Index %d, Tet Vertex %d) Face %d Vertex %d Edge Class (%d, %d)\n",
                       endpoint->region_index, endpoint->tri->tet_index,
                       endpoint->tri->tet_vertex, endpoint->face, endpoint->vertex,
                       endpoint->tri->vertices[endpoint->vertex].edge_class,
                       endpoint->tri->vertices[endpoint->vertex].edge_index);
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

                for (path_node = path->curves_begin.next; 
                    path_node != &path->curves_end; 
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
            g = cusp->dual_graph;
            for (j = 0; j < g->num_vertices; j++) {
                if (g->regions[j] == NULL)
                        continue;

                printf("    Vertex %d (Tet Index: %d, Tet Vertex: %d): ",
                       j,
                       g->regions[j]->tet_index,
                       g->regions[j]->tet_vertex
                );
                for (edge_node = g->edge_list_begin[j].next;
                    edge_node != &g->edge_list_end[j];
                    edge_node = edge_node->next)
                    printf("%d ", edge_node->y);

                printf("\n");
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

// ------------------------------------

/*
 * Train lines 
 *
 * Each cusp contains a number of edge classes which are split into two partitions, 
 * those in the end multi graph and those not in the end multi graph.
 *
 * On each cusp we want a curve (which we call a train line) which passes through
 * regions in the cusp which can dive to each edge class in the end multi graph.
 * This allows us to draw a curve between any two edge classes in the end multi
 * graph by walking along sections of the train line. The resulting curves have
 * intersection number 0 since the combinatorial holonomy on their intersection
 * is 0.
 *
 * Firstly we choose a collection of cusp edges for the train tracks. These will
 * be where the train line will dive into the manifold. We start with cusp 0, 
 * choose edges for this cusp and add the corresponding edges on the cusp 
 * at the other end of each edge class and the cusp index to a queue. Then iterate
 * until the queue is empty, only looking for edges if they have not been found 
 * yet. The edges are stored as path endpoints. A path endpoint which has a non-null
 * pointer for the tri attribute and a null pointer for the region attribute is an
 * edge. To find a curve and split we convert this path endpoint to one with a
 * valid region.
 *
 * On each cusp, find a path endpoint for the first edge and the second edge. Draw
 * a curve between them and split. Find a path endpoint for the third edge and 
 * find a path from the second endpoint to the third which does not backtrack in 
 * the first step, and split appropriately. Repeat until all edges are done.
 *
 */

void do_manifold_train_lines(CuspStructure **cusps, EndMultiGraph *multi_graph) {
    int cusp_index;

    find_edge_class_edges(cusps, multi_graph);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        find_primary_train_line(cusps, cusps[cusp_index], multi_graph);
        update_adj_curve_on_cusp(cusps[cusp_index]);
    }

    if (debug) {
        printf("\n");
        printf("Manifold Train Lines\n");
        printf("\n");
        printf("-------------------------------\n");

        log_structs(cusps[0]->manifold, cusps, NULL, "train_lines");
        log_structs(cusps[0]->manifold, cusps, NULL, "cusp_regions");
        log_structs(cusps[0]->manifold, cusps, NULL, "graph");
    }
}

/*
 * Since each manifold contains num_tetrahedra edge classes and each tetraheron 
 * contains atleast one edge class, we can pick an eage class from the edge 
 * classes on a tetrahedron with no two tetrahedra having the same edge class. 
 * This is done to simplify the process of constructing the train line. When we 
 * find endpoints for the train line, for a given edge class we ensure it lies
 * on the tetrahedron assigned to said edge class. This means on a given cusp, 
 * no two endpoints lie on the same cusp triangle.
 */

int *find_tet_index_for_edge_classes(Triangulation *manifold, EndMultiGraph *multi_graph) {
    int i, j, num_edge_classes = manifold->num_tetrahedra;
    int edge_source = 2 * num_edge_classes, tet_sink = 2 * num_edge_classes + 1;
    int *edge_class_to_tet_index = NEW_ARRAY(num_edge_classes, int);
    int **residual_network;
    Graph *g = init_graph(2 * num_edge_classes + 2, TRUE);
    Tetrahedron *tet;

    for (i = 0; i < num_edge_classes; i++)
        edge_class_to_tet_index[i] = -1;

    /*
     * We build a graph with a vertex for each edge class and tetrahedron,
     * and edges for edge class which lie in a given tetrahedron.
     * A bipartite matching then gives the assignment required.
     */

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {
        for (i = 0; i < 6; i++) {
            if (multi_graph->edge_classes[tet->edge_class[i]->index] || tet->edge_class[i]->index == multi_graph->e0)
                insert_edge(g, tet->edge_class[i]->index, num_edge_classes + tet->index, g->directed);
        }
    }

    /*
     * Convert the graph to a maximum flow problem
     */
    for (i = 0; i < num_edge_classes; i++)
        insert_edge(g, edge_source, i, g->directed);

    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next)
        insert_edge(g, num_edge_classes + tet->index, tet_sink, g->directed);

    residual_network = ford_fulkerson(g, edge_source, tet_sink);

    for (i = 0; i < num_edge_classes; i++) {
        for (j = num_edge_classes; j < 2 * num_edge_classes; j++)
            if (residual_network[j][i] == 1)
                edge_class_to_tet_index[i] = j - num_edge_classes;

        if ((multi_graph->edge_classes[i] || i == multi_graph->e0) && edge_class_to_tet_index[i] == -1) {
            uFatalError("find_tet_index_for_edge_classes", "symplectic_basis");
        }
    }

    for (i = 0; i < 2 * num_edge_classes + 2; i++)
        my_free(residual_network[i]);

    free_graph(g);
    my_free(residual_network);
    return edge_class_to_tet_index;
}

/*
 * Assign a cusp triangle, face and vertex to each PathEndPoint of the train 
 * line. This is done in a breadth first search fashion from the first cusp, 
 * adding cusps to the search queue after diving through them. We save E_0 for 
 * last since it takes into consideration the edge index, whereas all other edge 
 * classes in the end multi graph connect two disctint cusps. 
 */

void find_edge_class_edges(CuspStructure **cusps, EndMultiGraph *multi_graph) {
    int edge_class, cusp_index, other_cusp_index;
    int *edge_class_to_tet_index = find_tet_index_for_edge_classes(cusps[0]->manifold, multi_graph);
    Boolean found_edge_class;
    Boolean *visited_cusps, **edge_classes = NEW_ARRAY(multi_graph->num_cusps, Boolean *);
    Queue *queue = init_queue(multi_graph->num_cusps);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        edge_classes[cusp_index] = edge_classes_on_cusp(multi_graph, cusp_index);
    }

    enqueue(queue, 0);

    while (!empty_queue(queue)) {
        cusp_index = dequeue(queue);

        found_edge_class = FALSE;
        for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
            if (!edge_classes[cusp_index][edge_class])
                found_edge_class = TRUE;
        }

        if (!found_edge_class)
            continue;

        // assign edges to edge classes
        find_edge_class_edges_on_cusp(cusps[cusp_index],
                                      multi_graph,
                                      edge_classes[cusp_index],
                                      edge_class_to_tet_index);

        // update dive edges classes
        visited_cusps = update_edge_classes_on_cusp(cusps,
                                                    edge_classes,
                                                    multi_graph->num_cusps,
                                                    multi_graph->num_edge_classes,
                                                    cusp_index);

        for (other_cusp_index = 0; other_cusp_index < multi_graph->num_cusps; other_cusp_index++) {
            if (!visited_cusps[other_cusp_index])
                continue;

            enqueue(queue, other_cusp_index);
        }

        my_free(visited_cusps);
    }

    find_e0_edges_on_cusp(cusps, multi_graph, edge_class_to_tet_index);

    for (cusp_index = 0; cusp_index < multi_graph->num_cusps; cusp_index++) {
        my_free(edge_classes[cusp_index]);
    }

    my_free(edge_classes);
    my_free(edge_class_to_tet_index);
    free_queue(queue);
}


void find_edge_class_edges_on_cusp(CuspStructure *cusp, EndMultiGraph *multi_graph,
                                   const Boolean *edge_classes, const int *edge_class_to_tet_index) {
    int edge_class;
    VertexIndex v1, v2;
    FaceIndex face;
    CuspTriangle *tri;
    CuspVertex *vertex1, *vertex2;
    Boolean found;

    for (edge_class = 0; edge_class < multi_graph->num_edge_classes; edge_class++) {
        if (!edge_classes[edge_class])
            continue;

        if (edge_class_to_tet_index[edge_class] == -1)
            uFatalError("find_edge_class_edges_on_cusp", "symplectic_basis");

        found = FALSE;

        // find a cusp edge incident to the edge class 
        for (tri = cusp->cusp_triangle_begin.next; tri != &cusp->cusp_triangle_end; tri = tri->next) {
            if (found || edge_class_to_tet_index[edge_class] != tri->tet_index)
                continue;

            for (face = 0; face < 4; face++) {
                if (face == tri->tet_vertex || found)
                    continue;

                v1 = remaining_face[tri->tet_vertex][face];
                v2 = remaining_face[face][tri->tet_vertex];

                vertex1 = &tri->vertices[v1];
                vertex2 = &tri->vertices[v2];

                if (vertex1->edge_class == edge_class) {
                    cusp->train_line_endpoint[edge_class].tri = tri;
                    cusp->train_line_endpoint[edge_class].face = face;
                    cusp->train_line_endpoint[edge_class].vertex = v1;
                    found = TRUE;

                } else if (vertex2->edge_class == edge_class) {
                    cusp->train_line_endpoint[edge_class].tri = tri;
                    cusp->train_line_endpoint[edge_class].face = face;
                    cusp->train_line_endpoint[edge_class].vertex = v2;
                    found = TRUE;

                }
            }
        }
    }
}

void find_e0_edges_on_cusp(CuspStructure **cusps, EndMultiGraph *multi_graph, const int *edge_class_to_tet_index) {
    int cusp1, cusp2;
    VertexIndex vertex;
    CuspTriangle *tri;
    EdgeClass *edge;
    PathEndPoint *endpoint;
    Triangulation *manifold = cusps[0]->manifold;

    for (edge = manifold->edge_list_begin.next; edge != &manifold->edge_list_end; edge = edge->next) {
        if (edge->index != multi_graph->e0)
            continue;

        cusp1 = edge->incident_tet->cusp[one_vertex_at_edge[edge->incident_edge_index]]->index;
        cusp2 = edge->incident_tet->cusp[other_vertex_at_edge[edge->incident_edge_index]]->index;
    }

    for (tri = cusps[cusp1]->cusp_triangle_begin.next; tri != &cusps[cusp1]->cusp_triangle_end; tri = tri->next) {
        if (tri->tet_index != edge_class_to_tet_index[multi_graph->e0])
            continue;

        for (vertex = 0; vertex < 4; vertex++) {
            if (tri->tet_vertex == vertex)
                continue;

            if (tri->vertices[vertex].edge_class != multi_graph->e0 || tri->vertices[vertex].edge_index != START)
                continue;

            cusps[cusp1]->train_line_endpoint[multi_graph->e0].vertex = vertex;
            cusps[cusp1]->train_line_endpoint[multi_graph->e0].face = remaining_face[tri->tet_vertex][vertex];
            cusps[cusp1]->train_line_endpoint[multi_graph->e0].tri = tri;
        }
    }

    endpoint = &cusps[cusp1]->train_line_endpoint[multi_graph->e0];

    for (tri = cusps[cusp2]->cusp_triangle_begin.next; tri != &cusps[cusp2]->cusp_triangle_end; tri = tri->next) {
        if (tri->tet_index != edge_class_to_tet_index[multi_graph->e0] || tri->tet_vertex != endpoint->vertex)
            continue;

        vertex = endpoint->tri->tet_vertex;

        if (tri->vertices[vertex].edge_class != multi_graph->e0
            || (tri->vertices[vertex].edge_index != FINISH && cusp1 == cusp2))
            continue;

        if (cusp1 == cusp2) {
            cusps[cusp2]->extra_endpoint_e0.vertex = vertex;
            cusps[cusp2]->extra_endpoint_e0.face = endpoint->face;
            cusps[cusp2]->extra_endpoint_e0.tri = tri;
        } else {
            cusps[cusp2]->train_line_endpoint[multi_graph->e0].vertex = vertex;
            cusps[cusp2]->train_line_endpoint[multi_graph->e0].face = endpoint->face;
            cusps[cusp2]->train_line_endpoint[multi_graph->e0].tri = tri;
        }
    }
}

/*
 * Each edge we choose to add to the list of edges in find_edge_class_edges_on_cusp
 * has a corresponding edge on another cusp, which represents diving through the
 * manifold along that edge. Find these corresponding edges and set them in the
 * edge_begin and edges_end arrays, so the final train lines are consistent.
 */

Boolean *update_edge_classes_on_cusp(CuspStructure **cusps, Boolean **edge_classes,
                                     int num_cusps, int num_edge_classes, int current_cusp_index) {
    int cusp_index, other_cusp_index, edge_class;
    VertexIndex v1, v2, vertex;
    CuspStructure *cusp = cusps[current_cusp_index];
    Boolean *visited_cusp = NEW_ARRAY(num_edge_classes, Boolean);
    PathEndPoint *endpoint;
    CuspTriangle *tri, *adj_tri;

    for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
        visited_cusp[cusp_index] = FALSE;
    }

    for (edge_class = 0; edge_class < num_edge_classes; edge_class++) {
        if (!edge_classes[current_cusp_index][edge_class])
            continue;

        endpoint = &cusps[current_cusp_index]->train_line_endpoint[edge_class];

        for (cusp_index = 0; cusp_index < num_cusps; cusp_index++) {
            if (edge_classes[cusp_index][edge_class] && cusp_index != cusp->cusp->index)
                other_cusp_index = cusp_index;
        }

        v1 = remaining_face[endpoint->tri->tet_vertex][endpoint->face];
        v2 = remaining_face[endpoint->face][endpoint->tri->tet_vertex];

        if (endpoint->tri->vertices[v1].edge_class == edge_class)
            vertex = v1;
        else if (endpoint->tri->vertices[v2].edge_class == edge_class)
            vertex = v2;
        else 
            continue;
        

        adj_tri = NULL;

        for (tri = cusps[other_cusp_index]->cusp_triangle_begin.next;
             tri != &cusps[other_cusp_index]->cusp_triangle_end;
             tri = tri->next) {
            if (tri->tet_vertex != vertex || tri->tet_index != endpoint->tri->tet_index)
                continue;

            adj_tri = tri;
        }

        if (adj_tri == NULL)
            uFatalError("update_edge_classes_on_cusp", "symplectic_basis");

        edge_classes[other_cusp_index][edge_class] = FALSE;
        visited_cusp[other_cusp_index] = TRUE;

        cusps[other_cusp_index]->train_line_endpoint[edge_class].tri = adj_tri;
        cusps[other_cusp_index]->train_line_endpoint[edge_class].vertex = endpoint->tri->tet_vertex;
        cusps[other_cusp_index]->train_line_endpoint[edge_class].face = endpoint->face;
    }

    return visited_cusp;
}

/*
 * Use the regions on either side of the target edges to find a curve
 * through a cusp which passes along each target edge.
 */

void find_primary_train_line(CuspStructure **cusps, CuspStructure *cusp, EndMultiGraph *multi_graph) {
    int endpoint_start_index, endpoint_finish_index;
    PathEndPoint *start, *finish;

    endpoint_start_index = next_valid_endpoint_index(cusp->train_line_endpoint,
                                                     multi_graph->num_edge_classes,
                                                     -1);
    endpoint_finish_index = next_valid_endpoint_index(cusp->train_line_endpoint,
                                                      multi_graph->num_edge_classes,
                                                      endpoint_start_index);

    if (endpoint_start_index == -1 || (cusp->extra_endpoint_e0.tri == NULL && endpoint_finish_index == -1))
        return;

    if (cusp->extra_endpoint_e0.tri == NULL) {
        start = &cusp->train_line_endpoint[endpoint_start_index];
        finish = &cusp->train_line_endpoint[endpoint_finish_index];
    } else {
        start = &cusp->extra_endpoint_e0;
        finish = &cusp->train_line_endpoint[endpoint_start_index];

        endpoint_finish_index = endpoint_start_index;
    }

    tri_endpoint_to_region_endpoint(cusp, start);
    tri_endpoint_to_region_endpoint(cusp, finish);

    do_initial_train_line_segment_on_cusp(cusp, start, finish);

    endpoint_start_index = endpoint_finish_index;
    endpoint_finish_index = next_valid_endpoint_index(cusp->train_line_endpoint,
                                                      multi_graph->num_edge_classes,
                                                      endpoint_start_index);

    while (endpoint_finish_index != -1) {
        start = &cusp->train_line_endpoint[endpoint_start_index];
        finish = &cusp->train_line_endpoint[endpoint_finish_index];

        tri_endpoint_to_region_endpoint(cusp, finish);

        do_train_line_segment_on_cusp(cusp, start, finish);

        // update endpoint indices
        endpoint_start_index = endpoint_finish_index;
        endpoint_finish_index = next_valid_endpoint_index(cusp->train_line_endpoint,
                                                          multi_graph->num_edge_classes,
                                                          endpoint_start_index);
    }
}

/*
 * Construct the first segment of a train line. Essentially the same process
 * as do_curve_component_to_new_edge_class but stores the result in the cusp train
 * line.
 */

void do_initial_train_line_segment_on_cusp(CuspStructure *cusp, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    EdgeNode *node_begin, *node_end, *edge_node;
    Boolean *processed;
    Boolean *discovered;
    int *parent;

    node_begin  = NEW_STRUCT( EdgeNode );
    node_end    = NEW_STRUCT( EdgeNode );

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    construct_cusp_region_dual_graph(cusp);

    processed = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    init_search(cusp->dual_graph, processed, discovered, parent);
    bfs(cusp->dual_graph, start_endpoint->region_index, processed, discovered, parent);
    find_path(start_endpoint->region_index, finish_endpoint->region_index, parent, node_begin, node_end);

    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // split along curve
    graph_path_to_dual_curve(cusp->dual_graph, node_begin, node_end,
                             &cusp->train_line_path_begin, &cusp->train_line_path_end,
                             start_endpoint, finish_endpoint);
    split_cusp_regions_along_path(cusp, &cusp->train_line_path_begin,
                                  &cusp->train_line_path_end, start_endpoint, finish_endpoint);

    start_endpoint->node = cusp->train_line_path_begin.next;
    finish_endpoint->node = cusp->train_line_path_end.prev;

    while (node_begin->next != node_end) {
        edge_node = node_begin->next;
        REMOVE_NODE(edge_node);
        my_free(edge_node);
    }

    my_free(node_begin);
    my_free(node_end);
}

/*
 * Construct the next train line segment after the first. The start endpoint
 * is already set, so we set the region the start endpoint is in to visited
 * before starting the breadth first search and instead start from the region
 * adjacent across the face of the cusp triangle we dive along.
 */

void do_train_line_segment_on_cusp(CuspStructure *cusp, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    EdgeNode node_begin, node_end, *edge_node;
    PathNode *start_node;
    CuspRegion *region;
    Boolean cycle, *processed, *discovered;
    int start, start_index, cycle_start, cycle_end, visited, *parent;

    node_begin.next = &node_end;
    node_begin.prev = NULL;
    node_end.next   = NULL;
    node_end.prev   = &node_begin;

    construct_cusp_region_dual_graph(cusp);
    processed = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    init_search(cusp->dual_graph, processed, discovered, parent);

    start_index = tri_to_index(start_endpoint->tri->tet_index, start_endpoint->tri->tet_vertex);
    for (region = cusp->cusp_region_begin[start_index].next; region != &cusp->cusp_region_end[start_index]; region = region->next) {
        if (region->tet_index != start_endpoint->tri->tet_index || region->tet_vertex != start_endpoint->tri->tet_vertex)
            continue;

        if (!region->adj_cusp_triangle[start_endpoint->face] || !region->dive[start_endpoint->face][start_endpoint->vertex])
            continue;

        if (start_endpoint->face == cusp->train_line_path_end.prev->prev_face
            && region->curve[start_endpoint->face][start_endpoint->vertex] != 1)
            continue;

        start_endpoint->region_index = region->index;
        start_endpoint->region = region;
    }

    if (start_endpoint->region == NULL)
        uFatalError("do_train_line_segment_on_cusp", "symplectic_basis");

    /*
     * We require curves run between distinct sides of each cusp triangle
     * it enters. Hence, we need to remove the edge of the dual graph
     * corresponding to the last curve segment we drew. This edge will be
     * added back when the dual graph is reconstructed.
     */

    if (start_endpoint->face == cusp->train_line_path_end.prev->prev_face) {
        // curve dives along the face it passes through

        start = start_endpoint->region_index;
        visited = start_endpoint->region->adj_cusp_regions[start_endpoint->face]->index;
    } else {
        // curve dives through the vertex opposite the face it passes through or
        // curve travells around the vertex it dives through

        start = start_endpoint->region->adj_cusp_regions[start_endpoint->face]->index;
        visited = start_endpoint->region_index;
    }

    delete_edge(cusp->dual_graph, visited, start, cusp->dual_graph->directed);

    bfs(cusp->dual_graph, start, processed, discovered, parent);

    if (parent[finish_endpoint->region_index] == -1 && start != finish_endpoint->region_index) {
        /*
         * The finish endpoint is not in the subgraph we created by removing the edge
         * (visited, start). Assume there exists a cycle in this subgraph, we use this to
         * 'turn the curve around' and use the edge (visited, start).
         */

        init_search(cusp->dual_graph, processed, discovered, parent);
        cycle = cycle_exists(cusp->dual_graph, start, processed, discovered, parent, &cycle_start, &cycle_end);

        if (cycle == FALSE)
            // nothing we can do, train line does not work
            uFatalError("do_train_line_segment_on_cusp", "symplectic_basis");

        // reset parent array
        init_search(cusp->dual_graph, processed, discovered, parent);
        bfs(cusp->dual_graph, start, processed, discovered, parent);

        cycle_path(cusp->dual_graph, &node_begin, &node_end, processed, discovered, parent,
                   start, visited, finish_endpoint->region_index, cycle_start, cycle_end);
    } else {
        find_path(start, finish_endpoint->region_index, parent, &node_begin, &node_end);
    }

    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // split along curve
    start_node = cusp->train_line_path_end.prev;
    graph_path_to_path_node(cusp->dual_graph, &node_begin, &node_end,
                            &cusp->train_line_path_begin, &cusp->train_line_path_end, 
                            start_endpoint, finish_endpoint);
    split_cusp_regions_along_train_line_segment(cusp, start_node, &cusp->train_line_path_end, start_endpoint, finish_endpoint);

    finish_endpoint->node = cusp->train_line_path_end.prev;

    while (node_begin.next != &node_end) {
        edge_node = node_begin.next;
        REMOVE_NODE(edge_node);
        my_free(edge_node);
    }
}



int next_valid_endpoint_index(PathEndPoint *endpoints, int num_endpoints, int current_index) {
    int index;

    for (index = current_index + 1; index < num_endpoints; index++) {
        if (endpoints[index].tri == NULL)
            continue;

        return index;
    }

    return -1;
}

/*
 * Find a valid region for a path endpoint
 */

void tri_endpoint_to_region_endpoint(CuspStructure *cusp, PathEndPoint *endpoint) {
    CuspRegion *region;
    int index;

    index = tri_to_index(endpoint->tri->tet_index, endpoint->tri->tet_vertex);
    for (region = cusp->cusp_region_begin[index].next; region != &cusp->cusp_region_end[index]; region = region->next) {
        if (region->tet_index != endpoint->tri->tet_index || region->tet_vertex != endpoint->tri->tet_vertex)
            continue;

        if (!region->adj_cusp_triangle[endpoint->face] || !region->dive[endpoint->face][endpoint->vertex])
            continue;

        endpoint->region = region;
        endpoint->region_index = region->index;
    }

    if (endpoint->region == NULL)
        uFatalError("tri_endpoint_to_region_endpoint", "symplectic_basis");
}

void graph_path_to_path_node(Graph *g, EdgeNode *node_begin, EdgeNode *node_end, PathNode *path_begin,
                             PathNode *path_end, PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    FaceIndex face;
    EdgeNode *edge_node, *node;
    PathNode *path_node;
    CuspRegion *region;

    if (node_begin->next == node_end) {
        // path len 0
        return;
    } else if (node_begin->next->next == node_end) {
        // path len 1
        region = g->regions[node_begin->next->y];
        path_end->prev->next_face = start_endpoint->face;

        path_node = NEW_STRUCT( PathNode );
        INSERT_BEFORE(path_node, path_end);
        path_node->next_face = finish_endpoint->face;
        path_node->prev_face = EVALUATE(start_endpoint->tri->tet->gluing[start_endpoint->face], start_endpoint->face);
        path_node->cusp_region_index = node_begin->next->y;
        path_node->tri = region->tri;

        for (face = 0; face < 4; face++)
            if (region->tet_vertex != face &&
                path_node->next_face != face &&
                path_node->prev_face != face)
                break;

        path_node->inside_vertex = face;
        return;
    }

    // Set Header node
    path_end->prev->next_face = -1;

    // Add in a node for the start pos when the start endpoint is not in the same cusp tri as the first node.
    region = g->regions[node_begin->next->y];

    if (region->tet_index != start_endpoint->tri->tet_index || region->tet_vertex != start_endpoint->tri->tet_vertex) {
        node = NEW_STRUCT( EdgeNode );
        INSERT_AFTER(node, node_begin);
        node->y = start_endpoint->region_index;
    }

    for (face = 0; face < 4; face++) {
        if (!start_endpoint->region->adj_cusp_triangle[face])
            continue;

        if (start_endpoint->region->adj_cusp_regions[face]->index != g->regions[node_begin->next->next->y]->index)
            continue;

        path_end->prev->next_face = face;
    }

    if (path_end->prev->next_face == -1)
        uFatalError("graph_path_to_path_node", "symplectic_basis");

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next)
        interior_edge_node_to_path_node(g->regions[edge_node->y], path_end, edge_node);

    // Set Tail node
    endpoint_edge_node_to_path_node(g->regions[edge_node->y], path_end, edge_node, finish_endpoint, FINISH);
}

void split_cusp_regions_along_train_line_segment(CuspStructure *cusp, PathNode *path_begin, PathNode *path_end, 
                                                 PathEndPoint *start_endpoint, PathEndPoint *finish_endpoint) {
    int index = cusp->num_cusp_regions, split_type, region_index;
    PathNode *node;
    CuspRegion *p_region;
    Graph *g = cusp->dual_graph;

    if (path_begin->tri->tet_index == start_endpoint->tri->tet_index 
        && path_begin->tri->tet_vertex == start_endpoint->tri->tet_vertex) {
        node = path_begin;
    } else if (path_begin->next->tri->tet_index == start_endpoint->tri->tet_index 
               && path_begin->next->tri->tet_vertex == start_endpoint->tri->tet_vertex) {
        node = path_begin->prev;
    } else {
        uFatalError("split_cusp_regions_along_train_line_segment", "symplectic_basis");
        return;
    }

    if (node->next == path_end) {
        // empty path
        return ;
    }

    if (start_endpoint->face == node->prev_face) {
        // curve dives along the face it passes through
        split_type = 0;
    } else if (start_endpoint->vertex == node->prev_face) {
        // curve dives through the vertex opposite the face it passes through
        split_type = 1;
    } else {
        // curve travells around the vertex it dives through
        split_type = 2;
    }

    p_region = g->regions[start_endpoint->region_index];
    region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_train_line_endpoints(&cusp->cusp_region_begin[region_index],
                                              &cusp->cusp_region_end[region_index],
                                              p_region, node, start_endpoint, START);
    split_cusp_region_train_line_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                          node, start_endpoint, index, split_type);
    index++;

    // interior edges
    for (node = node->next; node->next != path_end; node = node->next) {
        p_region = g->regions[node->cusp_region_index];
        region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin[region_index],
                                           &cusp->cusp_region_end[region_index],
                                           p_region, node);
        split_cusp_region_path_interior(&cusp->cusp_region_end[region_index], p_region, node, index);
        index++;
    }

    // update last region
    p_region = g->regions[node->cusp_region_index];
    region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   p_region, finish_endpoint, node, FINISH);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                    node, finish_endpoint, index, FINISH);
    index++;

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
}

void update_cusp_triangle_train_line_endpoints(CuspRegion *cusp_region_start, CuspRegion *cusp_region_end,
                                               CuspRegion *region, PathNode *node, PathEndPoint *path_endpoint, int pos) {
    VertexIndex vertex1, vertex2;
    CuspRegion *current_region;

    vertex1 = remaining_face[region->tet_vertex][node->next_face];
    vertex2 = remaining_face[node->next_face][region->tet_vertex];

    for (current_region = cusp_region_start->next; current_region != cusp_region_end; current_region = current_region->next) {
        if (current_region == NULL || current_region->tet_index == -1)
            continue;

        // which triangle are we in?
        if (current_region->tet_index != region->tet_index || current_region->tet_vertex != region->tet_vertex)
            continue;

        if (!current_region->adj_cusp_triangle[node->next_face])
            continue;

        // Curve goes around the vertex or passes through the face that carries it
        if (current_region->curve[node->next_face][vertex1] > region->curve[node->next_face][vertex1]) {
            current_region->curve[node->next_face][vertex1]++;
            current_region->dive[node->next_face][vertex1] = FALSE;

        } else if (current_region->curve[node->next_face][vertex1] < region->curve[node->next_face][vertex1]) {
            current_region->curve[node->next_face][vertex2]++;
            current_region->dive[node->next_face][vertex2] = FALSE;
        }
    }
}

void split_cusp_region_train_line_endpoint(CuspRegion *region_end, CuspRegion *region, PathNode *node,
                                           PathEndPoint *path_endpoint, int index, int split_type) {
    VertexIndex vertex1, vertex2, other_vertex;
    CuspRegion *new_region = NEW_STRUCT(CuspRegion);

    copy_region(region, new_region);
    new_region->index = index;

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

    /*
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */

    if (split_type == 0) {
        if (node->next_face == path_endpoint->vertex) {
            // curve dives through the face opposite the next face
            other_vertex = (VertexIndex) (path_endpoint->face == vertex1 ? vertex2 : vertex1);

            new_region->curve[node->next_face][path_endpoint->face]++;
            new_region->dive[path_endpoint->face][other_vertex]   = region->dive[path_endpoint->face][other_vertex];
            new_region->dive[path_endpoint->vertex][other_vertex] = region->dive[path_endpoint->vertex][other_vertex];
            new_region->adj_cusp_triangle[other_vertex]           = FALSE;

            region->curve[node->next_face][other_vertex]++;
            region->dive[path_endpoint->face][other_vertex]       = FALSE;
            region->dive[path_endpoint->vertex][other_vertex]     = FALSE;
            region->adj_cusp_triangle[path_endpoint->face]        = FALSE;
        } else {
            new_region->curve[node->next_face][path_endpoint->face]++;
            new_region->dive[vertex1][path_endpoint->vertex]        = region->dive[vertex1][path_endpoint->vertex];
            new_region->dive[vertex2][path_endpoint->vertex]        = region->dive[vertex2][path_endpoint->vertex];
            new_region->adj_cusp_triangle[path_endpoint->face]      = FALSE;
            new_region->adj_cusp_triangle[path_endpoint->vertex]    = FALSE;

            region->curve[node->next_face][path_endpoint->vertex]++;
            region->dive[vertex1][path_endpoint->vertex]        = FALSE;
            region->dive[vertex2][path_endpoint->vertex]        = FALSE;
        }
    } else if (split_type == 1 || split_type == 2) {
        other_vertex = (VertexIndex) (path_endpoint->face == vertex1 ? vertex2 : vertex1);
        new_region->curve[path_endpoint->face][other_vertex]++;
        new_region->dive[path_endpoint->face][path_endpoint->vertex] = region->dive[path_endpoint->face][path_endpoint->vertex];
        new_region->adj_cusp_triangle[path_endpoint->vertex] = FALSE;
        new_region->adj_cusp_triangle[other_vertex] = FALSE;

        region->curve[path_endpoint->face][path_endpoint->vertex]++;
        region->dive[path_endpoint->face][path_endpoint->vertex] = FALSE;
    } else
        uFatalError("split_cusp_region_train_line_endpoint", "symplectic_basis");

    INSERT_BEFORE(new_region, region_end);
}

void cycle_path(Graph *g, EdgeNode *node_begin, EdgeNode *node_end,
                Boolean *processed, Boolean *discovered, int *parent,
                int start, int prev, int finish, int cycle_start, int cycle_end) {
    EdgeNode *node, *temp_node, temp_begin, temp_end;

    temp_begin.next = &temp_end;
    temp_begin.prev = NULL;
    temp_end.next = NULL;
    temp_end.prev = &temp_begin;

    // find a path from start -> cycle_end
    find_path(start, cycle_end, parent, node_begin, node_end);

    // duplicate the path start -> cycle_start, and reverse it
    find_path(start, cycle_start, parent, &temp_begin, &temp_end);
    for (node = temp_end.prev; node != &temp_begin; node = node->prev) {
        temp_node = NEW_STRUCT( EdgeNode );
        temp_node->y = node->y;
        INSERT_BEFORE(temp_node, node_end);
    }

    // free temp path
    while (temp_begin.next != &temp_end) {
        node = temp_begin.next;
        REMOVE_NODE(node);
        my_free(node);
    }

    // find a path from visited -> target
    init_search(g, processed, discovered, parent);
    bfs(g, prev, processed, discovered, parent);
    find_path(prev, finish, parent, node_end->prev, node_end);
}

// ------------------------------------

void do_oscillating_curves(CuspStructure **cusps, OscillatingCurves *curves, EndMultiGraph *multi_graph) {
    int i;

    for (i = 0; i < curves->num_curves; i++) {
        do_one_oscillating_curve(cusps, curves, multi_graph, curves->edge_class[i], i);

        if (debug) {
            printf("\n");
            printf("Oscillating Curve %d\n", i);
            printf("\n");
            printf("-------------------------------\n");

            log_structs(cusps[0]->manifold, cusps, curves, "dual_curves");
            log_structs(cusps[0]->manifold, cusps, curves, "endpoints");
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
                              int edge_class, int curve_index) {
    int path_length, orientation = START;
    CuspEndPoint *cusp_end_point, *endpoint;
    CurveComponent *path,
                   *curve_begin = &curves->curve_begin[curve_index],
                   *curve_end = &curves->curve_end[curve_index];

    cusp_end_point = find_multi_graph_path(multi_graph, cusps[0]->manifold, edge_class, &path_length);
    curve_begin->edge_class[FINISH] = edge_class;
    curve_end->edge_class[START]    = edge_class;

    path = setup_first_curve_component(cusps[cusp_end_point->cusp_index], multi_graph, cusp_end_point, curve_begin, curve_end);
    do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    update_path_holonomy(path, edge_class);

    if (debug) {
        printf("\n");
        printf("First Curve Segment\n");
        printf("\n");
        printf("-------------------------------\n");

        log_structs(cusps[0]->manifold, cusps, NULL, "cusp_regions");
        log_structs(cusps[0]->manifold, cusps, NULL, "graph");
    }

    // interior curve components, coming from train lines
    for (endpoint = cusp_end_point->next; endpoint->next != NULL; endpoint = endpoint->next) {
        orientation = (orientation == START ? FINISH : START);

        path = setup_train_line_component(cusps[endpoint->cusp_index], multi_graph, curve_begin, curve_end, endpoint, orientation);
        do_curve_component_on_train_line(cusps[path->cusp_index], path);
        update_path_holonomy(path, edge_class);
    }

    path = setup_last_curve_component(cusps[endpoint->cusp_index], multi_graph, endpoint, curve_begin, curve_end);
    do_curve_component_to_new_edge_class(cusps[path->cusp_index], path);
    update_path_holonomy(path, edge_class);

    update_adj_curve_along_path(cusps, curves, curve_index);

    if (debug) {
        printf("\n");
        printf("Last Curve Segment\n");
        printf("\n");
        printf("-------------------------------\n");

        log_structs(cusps[0]->manifold, cusps, NULL, "cusp_regions");
        log_structs(cusps[0]->manifold, cusps, NULL, "graph");
    }


    endpoint = cusp_end_point;
    while (endpoint != NULL) {
        cusp_end_point = endpoint;
        endpoint = endpoint->next;
        my_free(cusp_end_point);
    }
}

void copy_path_endpoint(PathEndPoint *endpoint1, PathEndPoint *endpoint2) {
    endpoint1->vertex = endpoint2->vertex;
    endpoint1->face = endpoint2->face;
    endpoint1->tri = endpoint2->tri;
    endpoint1->region_index = endpoint2->region_index;
    endpoint1->region = endpoint2->region;
    endpoint1->node = endpoint2->node;
    endpoint1->num_adj_curves = endpoint2->num_adj_curves;
}

PathNode *copy_path_node(PathNode *node) {
    PathNode *new_node = NEW_STRUCT( PathNode );

    new_node->next = NULL;
    new_node->prev = NULL;
    new_node->next_face = node->next_face;
    new_node->prev_face = node->prev_face;
    new_node->inside_vertex = node->inside_vertex;
    new_node->cusp_region_index = node->cusp_region_index;
    new_node->tri = node->tri;

    return new_node;
}

CurveComponent *setup_train_line_component(CuspStructure *cusp, EndMultiGraph *multi_graph, 
                                CurveComponent *curve_begin, CurveComponent *curve_end, 
                                CuspEndPoint *endpoint, int orientation) {
    CurveComponent *path; 

    path = init_curve_component(endpoint->edge_class[orientation],
                                endpoint->edge_class[orientation == START ? FINISH : START],
                                endpoint->cusp_index);

    INSERT_BEFORE(path, curve_end);

    if (path->edge_class[START] == multi_graph->e0 && orientation == START && cusp->extra_endpoint_e0.tri != NULL) {
        copy_path_endpoint(&path->endpoints[START], &cusp->extra_endpoint_e0);
        copy_path_endpoint(&path->endpoints[FINISH], &cusp->train_line_endpoint[path->edge_class[FINISH]]);
    } else if (path->edge_class[FINISH] == multi_graph->e0 && orientation == FINISH 
        && cusp->extra_endpoint_e0.tri != NULL) {
        copy_path_endpoint(&path->endpoints[START], &cusp->train_line_endpoint[path->edge_class[START]]);
        copy_path_endpoint(&path->endpoints[FINISH], &cusp->extra_endpoint_e0);
    } else {
        copy_path_endpoint(&path->endpoints[START], &cusp->train_line_endpoint[path->edge_class[START]]);
        copy_path_endpoint(&path->endpoints[FINISH], &cusp->train_line_endpoint[path->edge_class[FINISH]]);
    }

    return path;
}

/*
 * Find a curve along the train line and copy it to 'curve'
 */

void do_curve_component_on_train_line(CuspStructure *cusp, CurveComponent *curve) {
    int orientation = 0;
    PathNode *node, *new_node, *start_node, *finish_node;
    FaceIndex temp_face;

    start_node = curve->endpoints[START].node;
    finish_node = curve->endpoints[FINISH].node;

    for (node = cusp->train_line_path_begin.next; node != &cusp->train_line_path_end; node = node->next) {
        if (node == start_node) {
            orientation = 1;
            break;
        } else if (node == finish_node) {
            orientation = -1;
            break;
        }
    }

    if (orientation == 1) {
        // copy the train line into the curve
        for (node = start_node; node != finish_node->next; node = node->next) {
            new_node = copy_path_node(node);
            INSERT_BEFORE(new_node, &curve->curves_end);
        }
    } else if (orientation == -1) {
        // copy the train line into the curve
        for (node = start_node; node != finish_node->prev; node = node->prev) {
            new_node = copy_path_node(node);
            INSERT_BEFORE(new_node, &curve->curves_end);

            // reverse direction of faces
            temp_face = new_node->next_face;
            new_node->next_face = new_node->prev_face;
            new_node->prev_face = temp_face;
        }
    } else
        uFatalError("do_curve_component_on_train_line", "symplectic_basis");

    // correct endpoint inside vertices
    curve->curves_begin.next->prev_face = curve->endpoints[START].face;
    curve->curves_end.prev->next_face = curve->endpoints[FINISH].face;
}

CurveComponent *setup_first_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                            CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START],
                                endpoint->edge_class[FINISH],
                                endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end);

    construct_cusp_region_dual_graph(cusp);
    find_single_endpoint(cusp->dual_graph, &path->endpoints[START], path->edge_class[START], START);
    find_train_line_endpoint(cusp, &path->endpoints[FINISH], path->edge_class[FINISH],
                             START, multi_graph->e0, (Boolean) (endpoint->next->next != NULL));
    return path;
}

CurveComponent *setup_last_curve_component(CuspStructure *cusp, EndMultiGraph *multi_graph, CuspEndPoint *endpoint,
                                           CurveComponent *curves_begin, CurveComponent *curves_end) {
    CurveComponent *path;
    path = init_curve_component(endpoint->edge_class[START],
                                endpoint->edge_class[FINISH],
                                endpoint->cusp_index);
    INSERT_BEFORE(path, curves_end);

    construct_cusp_region_dual_graph(cusp);
    find_single_matching_endpoint(cusp->dual_graph, &curves_begin->next->endpoints[START],
                                  &path->endpoints[START], path->edge_class[START], FINISH);

    find_single_matching_endpoint(cusp->dual_graph, &path->prev->endpoints[FINISH],
                                  &path->endpoints[FINISH], path->edge_class[FINISH], FINISH);

    return path;
}

/*
 * Construct oscillating curves on the boundary components
 */

void do_curve_component_to_new_edge_class(CuspStructure *cusp, CurveComponent *curve) {
    int *parent;
    Boolean *processed, *discovered;
    EdgeNode *node_begin, *node_end, *edge_node;

    processed   = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    discovered  = NEW_ARRAY(cusp->dual_graph->num_vertices, Boolean);
    parent      = NEW_ARRAY(cusp->dual_graph->num_vertices, int);

    node_begin  = NEW_STRUCT( EdgeNode );
    node_end    = NEW_STRUCT( EdgeNode );

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    // Find curve using bfs
    init_search(cusp->dual_graph, processed, discovered, parent);
    bfs(cusp->dual_graph, curve->endpoints[START].region_index, processed, discovered, parent);

    find_path(curve->endpoints[START].region_index, curve->endpoints[FINISH].region_index,
              parent, node_begin, node_end);
    graph_path_to_dual_curve(cusp->dual_graph, node_begin, node_end,
                             &curve->curves_begin,&curve->curves_end,
                             &curve->endpoints[START], &curve->endpoints[FINISH]);

    // Reallocate memory
    my_free(processed);
    my_free(discovered);
    my_free(parent);

    // Split the regions along the curve
    split_cusp_regions_along_path(cusp, &curve->curves_begin, &curve->curves_end,
                                  &curve->endpoints[START], &curve->endpoints[FINISH]);

    while (node_begin->next != node_end) {
        edge_node = node_begin->next;
        REMOVE_NODE(edge_node);
        my_free(edge_node);
    }

    my_free(node_begin);
    my_free(node_end);
}


/*
 * Find a cusp region which can dive along a face into a vertex of
 * the cusp triangle which corresponds to 'edge_class' and 'edge_index',
 * and store the result in path_endpoint.
 */

void find_single_endpoint(Graph *g, PathEndPoint *path_endpoint, int edge_class, int edge_index) {
    int i;
    VertexIndex vertex;
    FaceIndex face1, face2, face;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL) {
            continue;
        }

        region = g->regions[i];
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
            path_endpoint->region_index     = i;
            path_endpoint->num_adj_curves   = region->num_adj_curves[path_endpoint->face][path_endpoint->vertex];

            return ;
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

void find_single_matching_endpoint(Graph *g, PathEndPoint *path_endpoint1, PathEndPoint *path_endpoint2, int edge_class, int edge_index) {
    int i;
    Boolean region_index, region_vertex, region_dive, region_curve;
    CuspRegion *region;

    // which cusp region
    for (i = 0; i < g->num_vertices; i++) {
        if (g->regions[i] == NULL)
            continue;

        region = g->regions[i];

        // are we in the matching endpoint
        region_index    = (Boolean) (region->tet_index != path_endpoint1->tri->tet_index);
        region_vertex   = (Boolean) (region->tet_vertex != path_endpoint1->vertex);
        region_dive     = (Boolean) !region->dive[path_endpoint1->face][path_endpoint1->tri->tet_vertex];
        region_curve    = (Boolean) (region->num_adj_curves[path_endpoint1->face][path_endpoint1->tri->tet_vertex] != path_endpoint1->num_adj_curves);

        if (region_index || region_vertex || region_dive || region_curve)
            continue;

        path_endpoint2->region          = region;
        path_endpoint2->tri             = region->tri;
        path_endpoint2->vertex          = path_endpoint1->tri->tet_vertex;
        path_endpoint2->face            = path_endpoint1->face;
        path_endpoint2->region_index    = i;
        path_endpoint2->num_adj_curves  = region->num_adj_curves[path_endpoint2->face][path_endpoint2->vertex];

        return ;
    }

    // didn't find valid path endpoints
    uFatalError("find_single_matching_endpoints", "symplectic_basis");
}

/*
 * find a path endpoint which matches the train line endpoint found during
 * do_manifold_train_lines().
 */

void find_train_line_endpoint(CuspStructure *cusp, PathEndPoint *endpoint, int edge_class, int edge_index,
                                       int e0, Boolean is_train_line) {
    int i;
    Boolean region_index, region_vertex, region_dive, region_curve;
    CuspRegion *region;
    PathEndPoint *train_line_endpoint;

    if (edge_class == e0 && edge_index == FINISH && cusp->extra_endpoint_e0.tri != NULL) {
        train_line_endpoint = &cusp->extra_endpoint_e0;
    } else {
        train_line_endpoint = &cusp->train_line_endpoint[edge_class];
    }

    // which cusp region
    for (i = 0; i < cusp->dual_graph->num_vertices; i++) {
        if (cusp->dual_graph->regions[i] == NULL)
            continue;

        region = cusp->dual_graph->regions[i];
        region_index    = (Boolean) (region->tet_index != train_line_endpoint->tri->tet_index);
        region_vertex   = (Boolean) (region->tet_vertex != train_line_endpoint->tri->tet_vertex);
        region_dive     = (Boolean) !region->dive[train_line_endpoint->face][train_line_endpoint->vertex];

        if (is_train_line) {
            region_curve    = (Boolean) (region->num_adj_curves[train_line_endpoint->face][train_line_endpoint->vertex] !=
                                         train_line_endpoint->num_adj_curves);
        } else {
            region_curve    = (Boolean) (region->num_adj_curves[train_line_endpoint->face][train_line_endpoint->vertex] != 0);
        }

        if (region_index || region_vertex || region_dive || region_curve)
            continue;

        endpoint->region          = region;
        endpoint->tri             = region->tri;
        endpoint->vertex          = train_line_endpoint->vertex;
        endpoint->face            = train_line_endpoint->face;
        endpoint->region_index    = region->index;
        endpoint->num_adj_curves  = region->num_adj_curves[endpoint->face][endpoint->vertex];

        return ;
    }

    uFatalError("find_train_line_endpoint", "symplectic_basis");
}

/*
 * After finding a path, each node contains the index of the region it lies in. 
 * Update path info calculates the face the path crosses to get to the next node 
 * and the vertex it cuts off to simplify combinatorial holonomy calculation.
 */

void graph_path_to_dual_curve(Graph *g, EdgeNode *node_begin, EdgeNode *node_end, PathNode *path_begin,
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
            if (g->regions[edge_node->y]->tet_vertex != face &&
                start_endpoint->vertex != face &&
                finish_endpoint->vertex != face)
                break;

        region = g->regions[edge_node->y];

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
    endpoint_edge_node_to_path_node(g->regions[edge_node->y], path_end, edge_node, start_endpoint, START);

    for (edge_node = node_begin->next->next; edge_node->next != node_end; edge_node = edge_node->next)
        interior_edge_node_to_path_node(g->regions[edge_node->y], path_end, edge_node);

    // Set Tail node
    endpoint_edge_node_to_path_node(g->regions[edge_node->y], path_end, edge_node, finish_endpoint, FINISH);
}

void endpoint_edge_node_to_path_node(CuspRegion *region, PathNode *path_end, EdgeNode *edge_node, PathEndPoint *path_endpoint, int pos) {
    FaceIndex face;
    VertexIndex vertex1, vertex2;
    PathNode *path_node = NEW_STRUCT( PathNode );
    path_node->cusp_region_index = edge_node->y;
    path_node->tri = region->tri;

    if (pos == START) {
        path_node->next_face = -1;
        for (face = 0; face < 4; face++) {
            if (face == region->tet_vertex || !region->adj_cusp_triangle[face] || path_node->next_face != -1)
                continue;

        if (region->adj_cusp_regions[face]->index == edge_node->next->y)
            path_node->next_face = face;
        }

//         next node isn't in an adjacent region
        if (path_node->next_face == -1)
            uFatalError("endpoint_edge_node_to_path_node", "symplectic_basis");

        path_node->prev_face = path_endpoint->face;
    } else {
        path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face], path_end->prev->next_face);
        path_node->next_face = path_endpoint->face;
    }

    vertex1 = remaining_face[region->tet_vertex][path_endpoint->vertex];
    vertex2 = remaining_face[path_endpoint->vertex][region->tet_vertex];

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

    path_node->prev_face = EVALUATE(path_end->prev->tri->tet->gluing[path_end->prev->next_face], path_end->prev->next_face);

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
    FaceIndex face;
    PathNode *node = path_begin->next;
    CuspRegion *region, *p_region, *current_region;
    Graph *g = cusp->dual_graph;

    // empty path
    if (node == path_end)
        return ;

    // path of len 1
    if (node->next == path_end) {
        region = NEW_STRUCT(CuspRegion);
        p_region = g->regions[node->cusp_region_index];
        region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
        INSERT_BEFORE(region, &cusp->cusp_region_end[region_index])
        copy_region(p_region, region);

        face = node->inside_vertex;

        region->index = index;
        region->adj_cusp_triangle[start_endpoint->vertex]  = FALSE;
        region->adj_cusp_triangle[finish_endpoint->vertex] = FALSE;
        region->dive[face][start_endpoint->vertex]  = TRUE;
        region->dive[face][finish_endpoint->vertex] = TRUE;
        region->dive[start_endpoint->vertex][finish_endpoint->vertex] = (Boolean) (face != finish_endpoint->face);
        region->dive[finish_endpoint->vertex][start_endpoint->vertex] = (Boolean) (face != start_endpoint->face);
        region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
        region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;

        p_region->adj_cusp_triangle[face]             = FALSE;
        p_region->dive[face][start_endpoint->vertex]  = (Boolean) (face == start_endpoint->face);
        p_region->dive[face][finish_endpoint->vertex] = (Boolean) (face == finish_endpoint->face);
        p_region->temp_adj_curves[face][start_endpoint->vertex]++;
        p_region->temp_adj_curves[face][finish_endpoint->vertex]++;

        // update other cusp regions
        for (current_region = cusp->cusp_region_begin[region_index].next;
            current_region != &cusp->cusp_region_end[region_index];
            current_region = current_region->next) {

           if (region->tet_index != current_region->tet_index || region->tet_vertex != current_region->tet_vertex)
               continue;

           if (current_region == region || current_region == p_region)
               continue;

           if (current_region->adj_cusp_triangle[start_endpoint->vertex] || current_region->adj_cusp_triangle[finish_endpoint->vertex]) {
               current_region->temp_adj_curves[face][finish_endpoint->vertex]++;
               current_region->temp_adj_curves[face][start_endpoint->vertex]++;

           } else {
               current_region->temp_adj_curves[start_endpoint->vertex][finish_endpoint->vertex]++;
               current_region->temp_adj_curves[finish_endpoint->vertex][start_endpoint->vertex]++;
           }
        }

        update_adj_region_data(cusp);
        cusp->num_cusp_regions++;
        return;
    }

    /*
     * Update first region
     *
     * Standing at the vertex where the curve dives through, and looking
     * at the opposite face, region becomes the cusp region to the right
     * of the curve and region to the left of the curve.
     */
    p_region = g->regions[node->cusp_region_index];
    region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   p_region, start_endpoint, node, START);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                    node, start_endpoint, index, START);
    index++;

    // interior edges
    while ((node = node->next)->next->next != NULL) {
        p_region = g->regions[node->cusp_region_index];
        region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
        update_cusp_triangle_path_interior(&cusp->cusp_region_begin[region_index],
                                           &cusp->cusp_region_end[region_index], p_region, node);
        split_cusp_region_path_interior(&cusp->cusp_region_end[region_index], p_region, node, index);
        index++;
    }

    // update last region
    p_region = g->regions[node->cusp_region_index];
    region_index = tri_to_index(p_region->tet_index, p_region->tet_vertex);
    update_cusp_triangle_endpoints(&cusp->cusp_region_begin[region_index],
                                   &cusp->cusp_region_end[region_index],
                                   p_region, finish_endpoint, node, FINISH);
    split_cusp_region_path_endpoint(&cusp->cusp_region_end[region_index], p_region,
                                    node, finish_endpoint, index, FINISH);
    index++;

    update_adj_region_data(cusp);
    cusp->num_cusp_regions = index;
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
     * Region becomes the cusp region closest to the inside vertex and
     * new_region becomes the cusp region on the other side of the oscillating curve
     */
    copy_region(region, new_region);
    new_region->index = index;

    // Update new region
    new_region->curve[v1][node->inside_vertex]++;
    new_region->curve[v2][node->inside_vertex]++;
    new_region->dive[node->inside_vertex][v1]       = region->dive[node->inside_vertex][v1];
    new_region->dive[v2][v1]                        = region->dive[v2][v1];
    new_region->dive[node->inside_vertex][v2]       = region->dive[node->inside_vertex][v2];
    new_region->dive[v1][v2]                        = region->dive[v1][v2];

    // Update region
    region->curve[v2][v1]++;
    region->curve[v1][v2]++;
    region->dive[v2][v1]                            = FALSE;
    region->dive[node->inside_vertex][v1]           = FALSE;
    region->dive[v1][v2]                            = FALSE;
    region->dive[node->inside_vertex][v2]           = FALSE;
    region->adj_cusp_triangle[node->inside_vertex]  = FALSE;

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
        new_region->dive[path_endpoint->face][path_endpoint->vertex]                         = region->dive[path_endpoint->face][path_endpoint->vertex];
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

    for (current_region = cusp_region_start->next; current_region != cusp_region_end; current_region = current_region->next) {
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

    for (current_region = cusp_region_start->next; current_region != cusp_region_end; current_region = current_region->next) {
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

                continue;
            }

            if (current_region->curve[path_endpoint->vertex][face1] > region->curve[path_endpoint->vertex][face1]) {
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

void update_adj_curve_along_path(CuspStructure **cusps, OscillatingCurves *curves, int curve_index) {
    int i, j;
    CurveComponent *curve,
               *dual_curve_begin = &curves->curve_begin[curve_index],
               *dual_curve_end   = &curves->curve_end[curve_index];
    CuspStructure *cusp;
    Triangulation *manifold = cusps[0]->manifold;

    // Update regions curve data
    for (curve = dual_curve_begin->next; curve != dual_curve_end; curve = curve->next)
        update_adj_curve_on_cusp(cusps[curve->cusp_index]);

    // update train line endpoints
    for (i = 0; i < manifold->num_cusps; i++) {
        cusp = cusps[i];

        if (cusp->extra_endpoint_e0.tri != NULL) {
            update_adj_curve_at_endpoint(&cusp->extra_endpoint_e0, dual_curve_begin->next, 1);
            update_adj_curve_at_endpoint(&cusp->extra_endpoint_e0, dual_curve_end->prev, 1);
        }

        for (j = 0; j < manifold->num_tetrahedra; j++) {
            if (cusp->train_line_endpoint[j].tri == NULL)
                continue;

            update_adj_curve_at_endpoint(&cusp->train_line_endpoint[j], dual_curve_begin->next, 1);
            update_adj_curve_at_endpoint(&cusp->train_line_endpoint[j], dual_curve_end->prev, 1);
        }
    }
}

/*
 * curve_begin and curve_end are header and tailer nodes of a doubly linked list of path
 * components for a new path. Update the path_endpoint->num_adj_curves attribute to account for this
 * new curve.
 */

void update_adj_curve_at_endpoint(PathEndPoint *path_endpoint, CurveComponent *path, int pos) {
    int i;
    PathEndPoint *curve_end_point;

    for (i = 0; i < 2; i++) {
        curve_end_point = &path->endpoints[i];

        if (curve_end_point->tri->tet_index != path_endpoint->tri->tet_index ||
            curve_end_point->tri->tet_vertex != path_endpoint->tri->tet_vertex ||
            curve_end_point->face != path_endpoint->face ||
            curve_end_point->vertex != path_endpoint->vertex)
            continue;

        path_endpoint->num_adj_curves++;
    }
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

    for (path_node = path->curves_begin.next; path_node != &path->curves_end; path_node = path_node->next) {
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->next_face]++;
        path_node->tri->tet->extra[edge_class].curve[path_node->tri->tet_vertex][path_node->prev_face]--;
    }
}

// ----------------------------------

/*
 * Construct the symplectic equations from the dual curves
 */

void calculate_holonomy(Triangulation *manifold, int **symp_eqns, int num_curves) {
    int curve, v, f, ff;
    Tetrahedron *tet;

    // which tet
    for (tet = manifold->tet_list_begin.next; tet != &manifold->tet_list_end; tet = tet->next) {

        // which curve
        for (curve = 0; curve < num_curves; curve++) {
            // which tet v
            for (v = 0; v < 4; v++) {

                for (f = 0; f < 4; f++) {
                    if (f == v)
                        continue;

                    ff = (int) remaining_face[v][f];

                    symp_eqns[curve][3 * tet->index + edge3_between_faces[f][ff]]
                                += FLOW(tet->extra[curve].curve[v][f], tet->extra[curve].curve[v][ff]);
                }
            }
        }
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

CuspEndPoint *find_multi_graph_path(EndMultiGraph *multi_graph, Triangulation *manifold, int edge_class, int *path_length) {
    Graph *g = multi_graph->multi_graph;
    Boolean *processed     = NEW_ARRAY(g->num_vertices, Boolean);
    Boolean *discovered    = NEW_ARRAY(g->num_vertices, Boolean);
    int *parent         = NEW_ARRAY(g->num_vertices, int);
    int start, end, startE0, endE0;
    EdgeNode *node_begin = NEW_STRUCT( EdgeNode ), *node_end = NEW_STRUCT( EdgeNode), *tempNode;
    CuspEndPoint *cusp_end_point;

    node_begin->next = node_end;
    node_begin->prev = NULL;
    node_end->next   = NULL;
    node_end->prev   = node_begin;

    find_edge_ends(g, manifold, edge_class, &start, &end);
    find_edge_ends(g, manifold, multi_graph->e0, &startE0, &endE0);

    init_search(g, processed, discovered, parent);
    bfs(g, start, processed, discovered, parent);

    *path_length = 0;
    *path_length = find_path_len(start, end, parent, *path_length);

    /*
     * Construct a path of even legnth through the end multi graph. 
     * As a convention the path stored in node_begin -> node_end is 
     * the cusps the path visits, the edge class of the edge is 
     * determined by the edge classes in the end multi graph and e0.
     */
    if (*path_length % 2 == 1) {
        find_path(start, end, parent, node_begin, node_end);
    } else {
        init_search(g, processed, discovered, parent);
        bfs(g, start, processed, discovered, parent);

        find_path(start, startE0, parent, node_begin, node_end);

        init_search(g, processed, discovered, parent);
        bfs(g, endE0, processed, discovered, parent);

        find_path(endE0, end, parent, node_end->prev, node_end);
    }

    cusp_end_point = graph_path_to_cusp_path(multi_graph, node_begin, node_end, edge_class);

    while (node_begin->next != node_end) {
        tempNode = node_begin->next;
        REMOVE_NODE(tempNode);
        my_free(tempNode);
    }

    my_free(node_begin);
    my_free(node_end);

    my_free(parent);
    my_free(discovered);
    my_free(processed);

    return cusp_end_point;
}

/*
 * Converts the EdgeNode path through the cusps in the end multigraph to 
 * a CuspEndPoint path which contains the edge classes on each cusp. 
 * A CuspEndPoint corresponding to one section of an oscillating curve, and
 * constructing such a section for all CuspEndPoints gives the whole curve
 *
 * node_begin is a doubly linked list with NULL header and tailer nodes and
 * the return CuspEndPoint is a pointer to a singly linked list with a NULL 
 * next value for the last node. 
 */

CuspEndPoint *graph_path_to_cusp_path(EndMultiGraph *multi_graph, EdgeNode *node_begin, EdgeNode *node_end, int edge_class) {
    int cusp, prev_edge_class;
    EdgeNode *node;
    CuspEndPoint *start_endpoint, *endpoint;

    start_endpoint = NEW_STRUCT( CuspEndPoint );
    start_endpoint->next = NULL;
    endpoint = start_endpoint;

    prev_edge_class = edge_class;
    for (node = node_begin->next; node->next != node_end; node = node->next) {
        cusp = node->y;

        endpoint->next = NEW_STRUCT( CuspEndPoint );
        endpoint->next->next = NULL;
        endpoint->cusp_index = cusp;
        endpoint->edge_class[START] = prev_edge_class;
        endpoint->edge_class[FINISH] = multi_graph->edges[node->y][node->next->y];

        if (endpoint->edge_class[FINISH] == -1)
            endpoint->edge_class[FINISH] = multi_graph->e0;
        
        prev_edge_class = endpoint->edge_class[FINISH];
        endpoint = endpoint->next;
    }

    endpoint->cusp_index = node->y;
    endpoint->edge_class[START] = prev_edge_class;
    endpoint->edge_class[FINISH] = edge_class;
    return start_endpoint;
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

Boolean *edge_classes_on_cusp(EndMultiGraph *multi_graph, int cusp_index) {
    int i;
    EdgeNode *edge_node;
    Boolean *edge_classes = NEW_ARRAY(multi_graph->num_edge_classes, Boolean);

    for (i = 0; i < multi_graph->num_edge_classes; i++)
        edge_classes[i] = FALSE;

    for (edge_node = multi_graph->multi_graph->edge_list_begin[cusp_index].next;
         edge_node != &multi_graph->multi_graph->edge_list_end[cusp_index];
         edge_node = edge_node->next) {
        edge_classes[multi_graph->edges[cusp_index][edge_node->y]] = TRUE;
    }

    return edge_classes;
}
