/**
 *  @file orb_diagram.h
 *
 *  This file defines data structures (OrbDiagram...) to encode a knotted
 *  graph as planar diagram with crossings - generalizing a knot diagram.
 *
 *  That is, it has edges that can cross each other. These crossings are
 *  recorded in OrbDiagramEdge.
 *
 *  orb_diagram_to_graph turns an OrbDiagram into an OrbGraph by
 *  turning each crossing into a four-valent vertex of a fat graph.
 *
 *  Ported from gui/diagram_canvas.h in Orb:
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h
 *
 */

#ifndef _orb_diagram_
#define _orb_diagram_

#include "SnapPea.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

typedef struct OrbDiagramEndData OrbDiagramEndData;
typedef struct OrbDiagramEdge OrbDiagramEdge;
typedef struct OrbDiagramVertex OrbDiagramVertex;
typedef struct OrbDiagramCrossing OrbDiagramCrossing;
typedef struct OrbDiagram OrbDiagram;

typedef struct OrbGraph OrbGraph;

/*
 * Corresponds to EndType in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L27
 */
enum OrbDiagramEndType
{
    diagramBegin = 0,
    diagramEnd
};

typedef enum OrbDiagramEndType OrbDiagramEndType;

/* Corresponds to EdgeType in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L28
 */
enum OrbDiagramEdgeType
{
    diagramSingular = 0,
    diagramDrilled
};

typedef enum OrbDiagramEdgeType OrbDiagramEdgeType;

/* Corresponds to EndData in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L34-L42
 */
struct OrbDiagramEndData
{
    OrbDiagramEdge    *edge;
    OrbDiagramEndType type;
    Boolean           singular;
    double            angle;
};

/* Corresponds to Vertex in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L44-L53
 */
struct OrbDiagramVertex
{
    int               x, y;
    int               connected_component;
    int               vertex_id;
    int               link_id;
    int               num_incident_end_data;
    OrbDiagramEndData **incident_end_data;
};

/* Corresponds to Edge in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L66-L85
 */
struct OrbDiagramEdge
{
    OrbDiagramVertex   *vertex[2];
    int                num_crossings;
    OrbDiagramCrossing **crossings;
    int                arc_id;
    int                link_id;
    int                edge_id;
    OrbDiagramEdgeType edge_type;
};

/* Corresponds to Crossing in gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L55-L63
 */
struct OrbDiagramCrossing
{
    int            x, y;
    int            crossing_id;
    int            crossing_sign;
    OrbDiagramEdge *over, *under;
    double         position_on_overstrand, position_on_understrand;
};

/* Corresponds to the diagram data stored by DiagramCanvas in
 * gui/diagram_canvas.h:
 * https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.h#L87-L127
 */
struct OrbDiagram
{
    int                num_arcs;
    int                num_links;
    int                num_vertices;
    OrbDiagramVertex   **vertices;
    int                num_edges;
    OrbDiagramEdge     **edges;
    int                num_crossings;
    OrbDiagramCrossing **crossings;
};

SNAPPEA_NAMESPACE_END_SCOPE

#endif
