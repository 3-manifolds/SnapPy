/**
 * diagram.h
 *
 * Data structures (prefixed by Diagram) to encode a knotted graph as
 * planar diagram with crossings (generalizes a knot diagram).
 *
 * Functions to convert a Diagram to a Graph and triangulate the complement.
 *
 * A Diagram has an embedding into the plane and its edges can cross.
 * The conversion to Graph turns each crossing into a four-valent vertex
 * results in a fat graph.
 *
 */

#ifndef _diagram_
#define _diagram_

#include "SnapPea.h"

typedef struct DiagramEndData DiagramEndData;
typedef struct DiagramEdge DiagramEdge;
typedef struct DiagramVertex DiagramVertex;
typedef struct DiagramCrossing DiagramCrossing;
typedef struct Diagram Diagram;
typedef enum DiagramEdgeType DiagramEdgeType;
typedef enum DiagramEndType DiagramEndType;

typedef struct Graph Graph;

/* Corresponds to DiagramCanvas::DiagramCanvas in gui/diagram_canvas.cpp */
void initialize_diagram(Diagram *);
/* Corresponds to DiagramCanvas::clearDiagram in gui/diagram_canvas.cpp */
void free_diagram(Diagram *);

void initialize_diagram_vertex(DiagramVertex *vertex);
void add_end_data_to_diagram_vertex(DiagramEndData * data, DiagramVertex * vertex);

void initialize_diagram_edge(DiagramEdge * edge);
void add_crossing_to_diagram_edge(DiagramCrossing * crossing, DiagramEdge * edge);

char * dump_diagram(Diagram * diagram);

/* Corresponds to DiagramCanvas::assign_arcs in gui/interface.cpp */
void assign_diagram_arcs(Diagram *);
/* Corresponds to DiagramCanvas::assign_links in gui/interface.cpp */
void assign_diagram_links(Diagram *);
/* Corresponds to DiagramCanvas::getCrossingSigns in gui/interface.cpp */
void assign_diagram_crossing_signs(Diagram * diagram);
/* Corresponds to DiagramCanvas::ed_angles in gui/interface.cpp */
void assign_diagram_end_data_angles(Diagram * diagram);
/* Corresponds to DiagramCanvas::assign_crossings_to_edges in gui/interface.cpp */
void assign_crossings_to_diagram_edges(Diagram * diagram);
/* Corresponds to DiagramCanvas::prepare_components_for_output in gui/interface.cpp */
void prepare_diagram_components_for_output(Diagram * diagram);

/* Corresponds to get_strand in gui/misc_functions.cpp */
int get_diagram_strand(DiagramEdge * e, DiagramVertex * v);

/* Corresponds to get_next_crossing in gui/misc_functions.cpp */
DiagramCrossing * get_next_diagram_crossing(DiagramEdge *e, DiagramCrossing *c);
/* Corresponds to get_prev_crossing in gui/misc_functions.cpp */
DiagramCrossing * get_prev_diagram_crossing(DiagramEdge *e, DiagramCrossing *c);

Graph * diagram_to_graph(Diagram *);
/* Corresponds to DiagramCanvas::outputTriangulation in gui/interface.cpp */
Triangulation * triangulate_diagram_complement(Diagram *, Boolean remove_vertices);

/* Corresponds to EndType in gui/diagram_canvas.h */
enum DiagramEndType
{
    diagramBegin = 0,
    diagramEnd
};

/* Corresponds to EdgeType in gui/diagram_canvas.h */
enum DiagramEdgeType
{
    diagramSingular = 0,
    diagramDrilled
};

/* Corresponds to EndData in gui/diagram_canvas.h */
struct DiagramEndData
{
    DiagramEdge *edge;
    DiagramEndType type;
    Boolean singular;
    double angle;
};

/* Corresponds to Vertex in gui/diagram_canvas.h */
struct DiagramVertex
{
    int x, y;
    int connected_component;
    int vertex_id;
    int link_id;
    int num_incident_end_data;
    DiagramEndData **incident_end_data;
};

/* Corresponds to Edge in gui/diagram_canvas.h */
struct DiagramEdge
{
    DiagramVertex *vertex[2];
    int num_crossings;
    DiagramCrossing **crossings;
    int arc_id;
    int link_id;
    int edge_id;
    DiagramEdgeType edge_type;
};

/* Corresponds to Crossing in gui/diagram_canvas.h */
struct DiagramCrossing
{
    int x, y;
    int crossing_id;
    int crossing_sign;
    DiagramEdge *over, *under;
    double position_on_overstrand, position_on_understrand;
};

/* Corresponds to DiagramCanvas in gui/diagram_canvas.h */
struct Diagram
{
    int num_arcs;
    int num_links;
    int num_vertices;
    DiagramVertex **vertices;
    int num_edges;
    DiagramEdge **edges;
    int num_crossings;
    DiagramCrossing **crossings;
};

#endif
