/**
 * diagram.h
 *
 * Data structures (prefixed by OrbDiagram) to encode a knotted graph as
 * planar diagram with crossings (generalizes a knot diagram).
 *
 * Functions to convert a OrbDiagram to a Graph and triangulate the complement.
 *
 * A OrbDiagram has an embedding into the plane and its edges can cross.
 * The conversion to Graph turns each crossing into a four-valent vertex
 * results in a fat graph.
 *
 */

#ifndef _orb_diagram_
#define _orb_diagram_

#include "SnapPea.h"

typedef struct OrbDiagramEndData OrbDiagramEndData;
typedef struct OrbDiagramEdge OrbDiagramEdge;
typedef struct OrbDiagramVertex OrbDiagramVertex;
typedef struct OrbDiagramCrossing OrbDiagramCrossing;
typedef struct OrbDiagram OrbDiagram;
typedef enum OrbDiagramEdgeType OrbDiagramEdgeType;
typedef enum OrbDiagramEndType OrbDiagramEndType;

typedef struct Graph Graph;

/* Corresponds to DiagramCanvas::DiagramCanvas in gui/diagram_canvas.cpp */
void orb_initialize_diagram(OrbDiagram *);
/* Corresponds to DiagramCanvas::clearDiagram in gui/diagram_canvas.cpp */
void orb_free_diagram(OrbDiagram *);

void orb_initialize_diagram_vertex(OrbDiagramVertex *vertex);
void orb_add_end_data_to_diagram_vertex(OrbDiagramEndData * data, OrbDiagramVertex * vertex);

void orb_initialize_diagram_edge(OrbDiagramEdge * edge);
void orb_add_crossing_to_diagram_edge(OrbDiagramCrossing * crossing, OrbDiagramEdge * edge);

char * orb_dump_diagram(OrbDiagram * diagram);

/* Corresponds to DiagramCanvas::assign_arcs in gui/interface.cpp */
void orb_assign_diagram_arcs(OrbDiagram *);
/* Corresponds to DiagramCanvas::assign_links in gui/interface.cpp */
void orb_assign_diagram_links(OrbDiagram *);
/* Corresponds to DiagramCanvas::getCrossingSigns in gui/interface.cpp */
void orb_assign_diagram_crossing_signs(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::ed_angles in gui/interface.cpp */
void orb_assign_diagram_end_data_angles(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::assign_crossings_to_edges in gui/interface.cpp */
void orb_assign_crossings_to_diagram_edges(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::prepare_components_for_output in gui/interface.cpp */
void orb_prepare_diagram_components_for_output(OrbDiagram * diagram);

/* Corresponds to get_strand in gui/misc_functions.cpp */
int orb_get_diagram_strand(OrbDiagramEdge * e, OrbDiagramVertex * v);

/* Corresponds to get_next_crossing in gui/misc_functions.cpp */
OrbDiagramCrossing * orb_get_next_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c);
/* Corresponds to get_prev_crossing in gui/misc_functions.cpp */
OrbDiagramCrossing * orb_get_prev_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c);

Graph * orb_diagram_to_graph(OrbDiagram *);
/* Corresponds to DiagramCanvas::outputTriangulation in gui/interface.cpp */
Triangulation * orb_triangulate_diagram_complement(OrbDiagram *, Boolean remove_vertices);

/* Corresponds to EndType in gui/diagram_canvas.h */
enum OrbDiagramEndType
{
    diagramBegin = 0,
    diagramEnd
};

/* Corresponds to EdgeType in gui/diagram_canvas.h */
enum OrbDiagramEdgeType
{
    diagramSingular = 0,
    diagramDrilled
};

/* Corresponds to EndData in gui/diagram_canvas.h */
struct OrbDiagramEndData
{
    OrbDiagramEdge *edge;
    OrbDiagramEndType type;
    Boolean singular;
    double angle;
};

/* Corresponds to Vertex in gui/diagram_canvas.h */
struct OrbDiagramVertex
{
    int x, y;
    int connected_component;
    int vertex_id;
    int link_id;
    int num_incident_end_data;
    OrbDiagramEndData **incident_end_data;
};

/* Corresponds to Edge in gui/diagram_canvas.h */
struct OrbDiagramEdge
{
    OrbDiagramVertex *vertex[2];
    int num_crossings;
    OrbDiagramCrossing **crossings;
    int arc_id;
    int link_id;
    int edge_id;
    OrbDiagramEdgeType edge_type;
};

/* Corresponds to Crossing in gui/diagram_canvas.h */
struct OrbDiagramCrossing
{
    int x, y;
    int crossing_id;
    int crossing_sign;
    OrbDiagramEdge *over, *under;
    double position_on_overstrand, position_on_understrand;
};

/* Corresponds to DiagramCanvas in gui/diagram_canvas.h */
struct OrbDiagram
{
    int num_arcs;
    int num_links;
    int num_vertices;
    OrbDiagramVertex **vertices;
    int num_edges;
    OrbDiagramEdge **edges;
    int num_crossings;
    OrbDiagramCrossing **crossings;
};

#endif
