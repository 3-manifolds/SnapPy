/**
 *  @file orb_diagram.c
 *
 *  Ported from:
 *  gui/diagram_canvas.cpp
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/diagram_canvas.cpp
 *  gui/interface.cpp
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/interface.cpp
 *  gui/misc_functions.cpp
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/gui/misc_functions.cpp
 */

#include "kernel.h"

#include <stdio.h>
#include <stdlib.h>

SNAPPEA_NAMESPACE_BEGIN_SCOPE


static void orb_add_crossing_to_diagram_edge(OrbDiagramCrossing * crossing, OrbDiagramEdge * edge);

static char * orb_dump_diagram(OrbDiagram * diagram);

/* Corresponds to DiagramCanvas::getCrossingSigns in gui/interface.cpp */
static void orb_assign_diagram_crossing_signs(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::ed_angles in gui/interface.cpp */
static void orb_assign_diagram_end_data_angles(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::assign_crossings_to_edges in gui/interface.cpp */
static void orb_assign_crossings_to_diagram_edges(OrbDiagram * diagram);
/* Corresponds to DiagramCanvas::prepare_components_for_output in gui/interface.cpp */
static void orb_prepare_diagram_components_for_output(OrbDiagram * diagram);

/* Corresponds to get_strand in gui/misc_functions.cpp */
static int orb_get_diagram_strand(OrbDiagramEdge * e, OrbDiagramVertex * v);

/* Corresponds to get_next_crossing in gui/misc_functions.cpp */
static OrbDiagramCrossing * orb_get_next_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c);
/* Corresponds to get_prev_crossing in gui/misc_functions.cpp */
static OrbDiagramCrossing * orb_get_prev_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c);


/* Corresponds to DiagramCanvas::DiagramCanvas in gui/diagram_canvas.cpp */
void orb_initialize_diagram(
    OrbDiagram * diagram)
{
    diagram->num_arcs = 0;
    diagram->num_links = 0;
    diagram->num_vertices = 0;
    diagram->vertices = NULL;
    diagram->num_edges = 0;
    diagram->edges = NULL;
    diagram->num_crossings = 0;
    diagram->crossings = NULL;
}

void orb_initialize_diagram_vertex(
    OrbDiagramVertex *vertex)
{
    vertex->vertex_id = 0;
    vertex->link_id = -1;
    vertex->num_incident_end_data = 0;
    vertex->incident_end_data = NULL;
}

static void orb_free_diagram_vertex(
    OrbDiagramVertex *vertex)
{
    int i;

    if (!vertex) {
	return;
    }

    if (vertex->incident_end_data) {
	for (i = 0; i < vertex->num_incident_end_data; ++i) {
	    my_free(vertex->incident_end_data[i]);
	}
	my_free(vertex->incident_end_data);
    }

    my_free(vertex);
}

/* Corresponds to DiagramCanvas::clearDiagram in gui/diagram_canvas.cpp */
void orb_free_diagram(
    OrbDiagram *diagram)
{
    int i;

    if (!diagram) {
	return;
    }

    if (diagram->vertices) {
	for (i = 0; i < diagram->num_vertices; ++i) {
	    orb_free_diagram_vertex(diagram->vertices[i]);
	}
	my_free(diagram->vertices);
    }

    if (diagram->edges) {
	for (i = 0; i < diagram->num_edges; ++i) {
	    my_free(diagram->edges[i]);
	}
	my_free(diagram->edges);
    }

    if (diagram->crossings) {
	for (i = 0; i < diagram->num_crossings; ++i) {
	    my_free(diagram->crossings[i]);
	}
	my_free(diagram->crossings);
    }

    my_free(diagram);
}

void orb_initialize_diagram_edge(
    OrbDiagramEdge * edge)
{
    edge->vertex[0] = NULL;
    edge->vertex[1] = NULL;
    edge->num_crossings = 0;
    edge->crossings = NULL;
    edge->arc_id = -1;
    edge->link_id = -1;
    edge->edge_id = -1;
    edge->edge_type = diagramSingular;
}

void orb_add_end_data_to_diagram_vertex(
    OrbDiagramEndData * data,
    OrbDiagramVertex * vertex)
{
    OrbDiagramEndData ** new_incident_end_data =
	NEW_ARRAY(vertex->num_incident_end_data + 1, OrbDiagramEndData*);

    for (int i = 0; i < vertex->num_incident_end_data; i++) {
	new_incident_end_data[i] = vertex->incident_end_data[i];
    }
    new_incident_end_data[vertex->num_incident_end_data++] = data;
    my_free(vertex->incident_end_data);
    vertex->incident_end_data = new_incident_end_data;
}

void orb_add_crossing_to_diagram_edge(
    OrbDiagramCrossing * crossing,
    OrbDiagramEdge * edge)
{
    OrbDiagramCrossing ** new_crossings =
	NEW_ARRAY(edge->num_crossings + 1, OrbDiagramCrossing*);

    for (int i = 0; i < edge->num_crossings; i++) {
        new_crossings[i] = edge->crossings[i];
    }
    new_crossings[edge->num_crossings++] = crossing;
    my_free(edge->crossings);
    edge->crossings = new_crossings;
}

/* Corresponds to DiagramCanvas::assign_arcs in gui/interface.cpp */
void orb_assign_diagram_arcs(
    OrbDiagram * diagram)
{
    diagram->num_arcs = 0;

    Boolean drilled_arc;

    int queue_begin = 0;
    int queue_end = 0;

    /* Each edge appears on the queue at most ones */
    OrbDiagramEdge ** queue = NEW_ARRAY(diagram->num_edges, OrbDiagramEdge *);

    Boolean * visited = NEW_ARRAY(diagram->num_edges, Boolean);

    for (int i = 0; i < diagram->num_edges; i++) {
	OrbDiagramEdge * edge = diagram->edges[i];
	edge->edge_id = i;
	edge->arc_id = -1;
	visited[i] = FALSE;
    }

    {
	OrbDiagramEdge * edge = diagram->edges[0];
	queue[queue_end++] = edge;
	visited[0] = TRUE;
	drilled_arc = edge->edge_type == diagramDrilled;
	if (!drilled_arc)
	{
	    edge->arc_id = diagram->num_arcs;
	}
    }

    while (queue_begin < queue_end)
    {
	while (queue_begin < queue_end)
	{
	    OrbDiagramEdge * e = queue[queue_begin++];
	    for (int i = 0; i < 2; i++)
	    {
		OrbDiagramVertex * v = e->vertex[i];

		if (e->vertex[i]->num_incident_end_data == 2)
		{
		    int j = (v->incident_end_data[0]->edge == e) ? 1 : 0;
		    OrbDiagramEdge *e1 = v->incident_end_data[j]->edge;
		    Boolean needsSwitch = v->incident_end_data[j]->type == i;
		    if (!visited[e1->edge_id])
		    {
			visited[e1->edge_id] = TRUE;

			if (!drilled_arc) {
			    e1->arc_id = diagram->num_arcs;
			}

			if (needsSwitch)
			{
			    OrbDiagramVertex *temp;
			    temp = e1->vertex[0];
			    e1->vertex[0] = e1->vertex[1];
			    e1->vertex[1] = temp;

			    for (int j = 0; j < 2; j++)
			    {
				for (int k = 0; k < e1->vertex[j]->num_incident_end_data; k++)
				{
				    if (e1->vertex[j]->incident_end_data[k]->edge == e1) {
					e1->vertex[j]->incident_end_data[k]->type =
                                            (OrbDiagramEndType)(                                           
                                                1 - e1->vertex[j]->incident_end_data[k]->type);
				    }
				}
			    }
			}

			queue[queue_end++] = e1;
		    }
		}
	    }
	}

	if (!drilled_arc) {
	    diagram->num_arcs++;
	}

	for (int i = 0; i < diagram->num_edges; i++)
	{
	    if (visited[i]) {
		continue;
	    }

	    OrbDiagramEdge * edge = diagram->edges[i];

	    queue[queue_end++] = edge;
	    visited[edge->edge_id] = TRUE;

	    drilled_arc = edge->edge_type == diagramDrilled;
	    if (!drilled_arc) {
		edge->arc_id = diagram->num_arcs;
	    }

	    break;
	}
    }

    my_free(visited);
    my_free(queue);
}

/* Corresponds to DiagramCanvas::assign_links in gui/interface.cpp */
void orb_assign_diagram_links(
    OrbDiagram * diagram)
{
    diagram->num_links = 0;

    int queue_begin = 0;
    int queue_end = 0;

    /* Each edge appears on the queue at most ones */
    OrbDiagramEdge ** queue = NEW_ARRAY(diagram->num_edges, OrbDiagramEdge *);

    Boolean * visited = NEW_ARRAY(diagram->num_edges, Boolean);

    for (int i = 0; i < diagram->num_edges; i++)
    {
	OrbDiagramEdge * edge = diagram->edges[i];
	edge->edge_id = i;
	edge->link_id = -2;
	visited[i] = edge->edge_type != diagramDrilled;

	if (queue_begin == queue_end && !visited[i])
	{
	    queue[queue_end++] = edge;
	    visited[i] = TRUE;
	    edge->link_id = diagram->num_links;
	}
    }

    while (queue_begin < queue_end)
    {
	while (queue_begin < queue_end)
	{
	    OrbDiagramEdge * e = queue[queue_begin++];
	    for (int i = 0; i < 2; i++)
	    {
		OrbDiagramVertex * v = e->vertex[i];

		for (int j = 0; j < v->num_incident_end_data; j++)
		{
		    OrbDiagramEdge *e1 = v->incident_end_data[j]->edge;
		    if (!visited[e1->edge_id])
		    {
			visited[e1->edge_id] = TRUE;
			e1->link_id = diagram->num_links;
			queue[queue_end++] = e1;
		    }
		}
	    }
	}

	diagram->num_links++;

	for (int i = 0; i < diagram->num_edges; i++)
	{
	    if (visited[i]) {
		continue;
	    }
	    OrbDiagramEdge * e = diagram->edges[i];
	    queue[queue_end++] = e;
	    e->link_id = diagram->num_links;
	}
    }

    for (int i = 0; i < diagram->num_vertices; i++)
    {
	OrbDiagramVertex * v = diagram->vertices[i];
	v->link_id = -1;
	for (int j = 0; j < v->num_incident_end_data; j++)
	{
	    if (v->incident_end_data[j]->edge->edge_type == diagramDrilled)
	    {
		v->link_id = v->incident_end_data[j]->edge->link_id;
		break;
	    }
	}
	if (v->link_id >= 0) {
	    continue;
	}
	if (v->num_incident_end_data > 2)
	{
	    v->link_id = diagram->num_links++;
	    continue;
	}
	if (v->num_incident_end_data == 2)
	{
	    Boolean singular_loop = TRUE;
	    int arc = v->incident_end_data[0]->edge->arc_id;

	    for (int j = 0; j < diagram->num_edges; j++)
	    {
		OrbDiagramEdge * e = diagram->edges[j];
		if (e->arc_id != arc) {
		    continue;
		}
		for (int k = 0; k < 2; k++) {
		    OrbDiagramVertex * v1 = e->vertex[k];
		    if (v1->num_incident_end_data != 2 ||
			v1->link_id > -1 ||
			v1->incident_end_data[0]->edge->edge_type == diagramDrilled ||
			v1->incident_end_data[1]->edge->edge_type == diagramDrilled)
		    {
			singular_loop = FALSE;
			break;
		    }
		}
		if (!singular_loop) {
		    break;
		}
	    }
	    if (singular_loop) {
		v->link_id = diagram->num_links++;
	    }
	}
    }

    my_free(visited);
    my_free(queue);
}

char * orb_dump_diagram(
    OrbDiagram * diagram)
{
    size_t size = 10000000;

    char * buffer = NEW_ARRAY(size, char);
    char * p = buffer;
    char * end = buffer + size - 1;

    p += snprintf(p, end - p, "num_arcs = %d\n", diagram->num_arcs);
    p += snprintf(p, end - p, "num_links = %d\n", diagram->num_links);

    for (int i = 0; i < diagram->num_vertices; i++)
    {
	OrbDiagramVertex * v = diagram->vertices[i];
	p += snprintf(p, end - p, "Vertex:\n");
	p += snprintf(p, end - p, "    %d %d\n", v->x, v->y);
	p += snprintf(p, end - p, "    %d %d %d\n", v->connected_component, v->vertex_id, v->link_id);
	for (int j = 0; j < v->num_incident_end_data; j++) {
	    OrbDiagramEndData * e = v->incident_end_data[j];
	    p += snprintf(p, end - p, "    End data:\n");
	    p += snprintf(p, end - p, "        Edge: %d\n", e->edge->edge_id);
	    p += snprintf(p, end - p, "        %d %d %lf\n",
			  e->type, e->singular, e->angle);
	}
    }

    for (int i = 0; i < diagram->num_edges; i++)
    {
	OrbDiagramEdge * e = diagram->edges[i];
	p += snprintf(p, end - p, "Edge:\n");
	p += snprintf(p, end - p, "    %d %d\n", e->vertex[0]->vertex_id, e->vertex[1]->vertex_id);
	p += snprintf(p, end - p, "    %d %d %d %d\n", e->edge_id, e->arc_id, e->link_id, e->edge_type);
	for (int j = 0; j < e->num_crossings; j++) {
	    OrbDiagramCrossing * c = e->crossings[j];
	    p += snprintf(p, end - p, "    Crossing %d\n", c->crossing_id);
	}
    }

    for (int i = 0; i < diagram->num_crossings; i++)
    {
	OrbDiagramCrossing * c = diagram->crossings[i];
	p += snprintf(p, end - p, "Crossing:\n");
	p += snprintf(p, end - p, "    %d\n", c->crossing_id);
	p += snprintf(p, end - p, "    %d %d  %d\n", c->x, c->y, c->crossing_sign);
	p += snprintf(p, end - p, "    %d %d\n", c->over->edge_id, c->under->edge_id);
	p += snprintf(p, end - p, "    %lf %lf\n", c->position_on_overstrand, c->position_on_understrand);
    }

    return buffer;
}

static
Complex
orb_point_position_to_complex(
    OrbDiagramVertex * v)
{
    Complex r;
    r.real = v->x;
    r.imag = v->y;
    return r;
}

void
orb_assign_diagram_crossing_signs(
    OrbDiagram * diagram)
{
    for (int i = 0; i < diagram->num_crossings; i++)
    {
        OrbDiagramCrossing *c = diagram->crossings[i];

        Complex z1 = orb_point_position_to_complex(c->over->vertex[diagramBegin]);
        Complex z2 = orb_point_position_to_complex(c->over->vertex[diagramEnd]);
        Complex z3 = orb_point_position_to_complex(c->under->vertex[diagramBegin]);

        /* w = ( z3 - z1 ) / ( z2 - z1 ); */
	Complex w = complex_div(complex_minus(z3, z1), complex_minus(z2, z1));

        c->crossing_sign = ( w.imag > 0 ) ? 1 : -1;
    }
}

void
orb_assign_diagram_end_data_angles(
    OrbDiagram * diagram)
{
    for (int i = 0; i < diagram->num_vertices; i++)
    {
        OrbDiagramVertex * v = diagram->vertices[i];
	for (int j = 0; j < v->num_incident_end_data; j++)
	{
  	    OrbDiagramEndData * e = v->incident_end_data[j];
	    Complex z =
  	        complex_minus(
		    orb_point_position_to_complex(e->edge->vertex[diagramEnd]),
		    orb_point_position_to_complex(e->edge->vertex[diagramBegin]));
	    if (e->type == diagramEnd) {
    	        z = complex_negate(z);
	    }
	    e->angle = -atan2(z.imag, z.real);/* minus???? */
	}
    }
}

void
orb_assign_crossings_to_diagram_edges(
    OrbDiagram * diagram)
{
    for (int i = 0; i < diagram->num_crossings; i++)
    {
        OrbDiagramCrossing * c = diagram->crossings[i];
	orb_add_crossing_to_diagram_edge(c, c->over);
	orb_add_crossing_to_diagram_edge(c, c->under);
    }
}

void
orb_prepare_diagram_components_for_output(
    OrbDiagram * diagram)
{
    {
        int     link_id = -5;

	for (int i = 0; i < diagram->num_arcs; i++)
	{
    	    for (int j = 0; j < diagram->num_edges; j++)
	    {
	        OrbDiagramEdge *e = diagram->edges[j];
		if (e->arc_id == i)
		{
  		    for(int k = 0; k < 2; k++)
		    {
		        if (e->vertex[k]->link_id != -1 )
		        {
			    e->vertex[k]->incident_end_data[orb_get_diagram_strand(e,e->vertex[k])]->singular=(k==1);
			    if (k == 0) {
			        link_id = e->vertex[k]->link_id;
			    }
			}
		    }
		}
	    }
    	    for (int j = 0; j < diagram->num_edges; j++)
	    {
	        OrbDiagramEdge *e = diagram->edges[j];
		if (e->arc_id == i) {
  		    e->link_id = link_id;
		}
	    }
	}
    }

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        OrbDiagramVertex *v = diagram->vertices[i];
        if (v->link_id == -1 )
	{
  	    for( int j = 0; j < v->num_incident_end_data; j++)
	    {
	        if (v->incident_end_data[j]->edge->edge_type == diagramDrilled)
		{
 		    v->link_id = v->incident_end_data[j]->edge->link_id;
		}
	    }
	}
    }
}

int orb_get_diagram_strand(OrbDiagramEdge * e, OrbDiagramVertex * v)
{
    for (int i = 0; i < v->num_incident_end_data; i++)
    {
        if (e == v->incident_end_data[i]->edge) {
	    return i;
        }
    }

    fprintf(stderr, "ERR: diagram_get_strand\n");
    return -1;
}

OrbDiagramCrossing * orb_get_next_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c)
{
    int i;
    for (i = 0; i < e->num_crossings; i++) {
	if (e->crossings[i] == c) {
	    break;
	}
    }

    if (i + 1 < e->num_crossings) {
	return e->crossings[i + 1];
    }
    return NULL;
}

OrbDiagramCrossing * orb_get_prev_diagram_crossing(OrbDiagramEdge *e, OrbDiagramCrossing *c)
{
    int i;
    for (i = 0; i < e->num_crossings; i++) {
	if (e->crossings[i] == c) {
	    break;
	}
    }

    if (0 < i && i < e->num_crossings) {
	return e->crossings[i - 1];
    }
    return NULL;
}

static int cmp(const double d1, const double d2)
{
    if (d1 < d2) {
	return -1;
    }
    if (d1 > d2) {
	return +1;
    }
    return 0;
}

EXTERN_C_BEGIN_SCOPE
static int orb_ed_more(const void *d1, const void *d2)
{
    /* return ed1->angle > ed2->angle; */

    OrbDiagramEndData * ed1 = *(OrbDiagramEndData **)d1;
    OrbDiagramEndData * ed2 = *(OrbDiagramEndData **)d2;

    return cmp(ed2->angle, ed1->angle);
}

static int orb_crossing_less(const void *d1, const void *d2)
{
    OrbDiagramCrossing *c1 = *(OrbDiagramCrossing **)d1;
    OrbDiagramCrossing *c2 = *(OrbDiagramCrossing **)d2;

    if (c1->over == c2->over)
    {
	return cmp(c1->position_on_overstrand, c2->position_on_overstrand);
    }

    if (c1->over == c2->under)
    {
	return cmp(c1->position_on_overstrand, c2->position_on_understrand);
    }

    if (c1->under == c2->over)
    {
	return cmp(c1->position_on_understrand, c2->position_on_overstrand);
    }

    if (c1->under == c2->under)
    {
	return cmp(c1->position_on_understrand, c2->position_on_understrand);
    }

    return 0;
}

EXTERN_C_END_SCOPE

void orb_remove_meeting( OrbGraph *graph, int index )
{
    if ( index>= graph->num_meetings) {
	return;
    }

    OrbGraphMeeting *old_array = graph->meeting;
    OrbGraphMeeting *new_array = NEW_ARRAY(graph->num_meetings-1, OrbGraphMeeting );

    for (int i = 0, j=0; i< graph->num_meetings-1; i++ )
    {
        if (i == index) {
	    j++;
	}

        new_array[i].type = old_array[i+j].type;
	new_array[i].handedness = old_array[i+j].handedness;
        new_array[i].num_strands = old_array[i+j].num_strands;
        new_array[i].label = old_array[i+j].label;
        new_array[i].strand = old_array[i+j].strand;
        new_array[i].component = old_array[i+j].component;
        new_array[i].neighbor = old_array[i+j].neighbor;
    }

    my_free(old_array);
    graph->meeting = new_array;
    graph->num_meetings--;

    for(int i = 0; i < graph->num_meetings; i++) {
	for(int j=0; j < graph->meeting[i].num_strands; j++) {
	    if (graph->meeting[i].neighbor[j] >= index) {
		graph->meeting[i].neighbor[j]--;
	    }
	}
    }
}

/* Ported from OrbDiagramCanvas::outputTriangulation in interface.cpp */
OrbGraph * orb_diagram_to_graph(
    OrbDiagram * diagram)
{
    int num_meetings = 0;

    orb_assign_diagram_crossing_signs(diagram);
    orb_assign_diagram_end_data_angles(diagram);
    orb_assign_crossings_to_diagram_edges(diagram);
    orb_prepare_diagram_components_for_output(diagram);

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        diagram->vertices[i]->vertex_id = num_meetings++;
	qsort(diagram->vertices[i]->incident_end_data,
	      diagram->vertices[i]->num_incident_end_data,
	      sizeof(OrbDiagramEndData*),
	      &orb_ed_more);
    }

    for(int i = 0; i < diagram->num_crossings; i++) {
        diagram->crossings[i]->crossing_id = num_meetings++;
    }

    for(int i = 0; i < diagram->num_edges; i++) {
        qsort(diagram->edges[i]->crossings,
	      diagram->edges[i]->num_crossings,
	      sizeof(OrbDiagramCrossing*),
	      &orb_crossing_less);
    }

    OrbGraph *graph = NEW_STRUCT(OrbGraph);
    graph->num_meetings = num_meetings;
    graph->num_components = diagram->num_links;
    graph->num_free_loops = 0;
    graph->meeting = NEW_ARRAY( num_meetings, OrbGraphMeeting );

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        OrbGraphMeeting *meeting = &graph->meeting[i];
	OrbDiagramVertex *vertex = diagram->vertices[i];

	meeting->type = Inter;
	meeting->num_strands = vertex->num_incident_end_data;
	meeting->strand = NEW_ARRAY(meeting->num_strands,int);
	meeting->neighbor = NEW_ARRAY(meeting->num_strands,int);
	meeting->component = NEW_ARRAY(meeting->num_strands,int);
	meeting->label = NEW_ARRAY(meeting->num_strands,int);
	meeting->tet = NULL;
	meeting->handedness = 0;

	for(int j = 0; j < meeting->num_strands; j++)
        {
	    OrbDiagramEndData *ed = vertex->incident_end_data[j];
	    OrbDiagramEdge *e = ed->edge;
	    meeting->label[j] = (ed->singular) ? e->arc_id : -1;

	    meeting->component[j] = diagram->vertices[i]->link_id;

	    if (ed->type == diagramBegin)
	    {
	        if (e->num_crossings == 0)
		{
		    meeting->strand[j] = orb_get_diagram_strand(e, e->vertex[diagramEnd]);
		    meeting->neighbor[j] = e->vertex[diagramEnd]->vertex_id;
		}
		else
		{
  		    OrbDiagramCrossing *c = e->crossings[0];
		    meeting->neighbor[j] = c->crossing_id;

		    if (c->over == e) {
		        meeting->strand[j] = 0;
		    } else {
		        meeting->strand[j] = 1;
		    }
                }
	    }
	    else
	    {
	        if (e->num_crossings == 0)
		{
		    meeting->strand[j] = orb_get_diagram_strand(e, e->vertex[diagramBegin]);
		    meeting->neighbor[j] = e->vertex[diagramBegin]->vertex_id;
		}
		else
		{
		    OrbDiagramCrossing *c = e->crossings[e->num_crossings-1];
		    meeting->neighbor[j] = c->crossing_id;

		    if (c->over == e) {
		        meeting->strand[j] = 2;
		    } else {
		        meeting->strand[j] = 3;
		    }

		}
	    }
	}
    }

    for(int i = 0; i < diagram->num_crossings; i++)
    {
	OrbGraphMeeting *meeting = &graph->meeting[diagram->num_vertices+i];
	OrbDiagramCrossing *c = diagram->crossings[i], *other;

	meeting->type = Cross;
	meeting->num_strands = 4;
	meeting->strand = NEW_ARRAY( 4, int );
	meeting->neighbor = NEW_ARRAY( 4, int );
	meeting->component = NEW_ARRAY( 4, int );
	meeting->label = NEW_ARRAY( 4, int );
	meeting->tet = NULL;
	meeting->handedness = c->crossing_sign;

	for (int j = 0; j < 4; j++) {
	    meeting->label[j] = -1;
	}

	OrbDiagramEdge * e = c->over;

	meeting->component[0] = e->link_id;
	meeting->component[2] = e->link_id;

	if ( (other = orb_get_prev_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[0] = orb_get_diagram_strand( e, e->vertex[diagramBegin] );
	    meeting->neighbor[0] = e->vertex[diagramBegin]->vertex_id;
	}
	else
	{
	    meeting->neighbor[0] = other->crossing_id;
	    meeting->strand[0] = (other->over==e) ? 2 : 3;
	}

	if ( (other = orb_get_next_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[2] = orb_get_diagram_strand( e, e->vertex[diagramEnd] );
	    meeting->neighbor[2] = e->vertex[diagramEnd]->vertex_id;
	}
	else
	{
	    meeting->neighbor[2] = other->crossing_id;
	    meeting->strand[2] = (other->over==e) ? 0 : 1;
	}

	e = c->under;

	meeting->component[1] = e->link_id;
	meeting->component[3] = e->link_id;

	if ( (other = orb_get_prev_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[1] = orb_get_diagram_strand( e, e->vertex[diagramBegin] );
	    meeting->neighbor[1] = e->vertex[diagramBegin]->vertex_id;
	}
	else
	{
	    meeting->neighbor[1] = other->crossing_id;
	    meeting->strand[1] = (other->over==e) ? 2 : 3;
	}

	if ( (other = orb_get_next_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[3] = orb_get_diagram_strand( e, e->vertex[diagramEnd] );
	    meeting->neighbor[3] = e->vertex[diagramEnd]->vertex_id;
	}
	else
	{
	    meeting->neighbor[3] = other->crossing_id;
	    meeting->strand[3] = (other->over==e) ? 0 : 1;
	}

    }

    for(int i = 0; i < diagram->num_crossings; i++)
    {
	if (diagram->crossings[i]->crossing_sign == 1)
	{
	    int temp_strand , temp_neighbor;
	    OrbGraphMeeting *meeting, *nbr1, *nbr3;

	    meeting = &graph->meeting[diagram->crossings[i]->crossing_id];

	    nbr1 = &graph->meeting[meeting->neighbor[1]];
	    nbr3 = &graph->meeting[meeting->neighbor[3]];

	    nbr1->strand[ meeting->strand[1] ] = 3;
	    nbr3->strand[ meeting->strand[3] ] = 1;

	    temp_strand = meeting->strand[3];
	    temp_neighbor = meeting->neighbor[3];

	    meeting->strand[3] = meeting->strand[1];
	    meeting->neighbor[3] = meeting->neighbor[1];

	    meeting->strand[1] = temp_strand;
	    meeting->neighbor[1] = temp_neighbor;
	}
    }


    int i = 0;
    while (i < graph->num_meetings)
    {
	OrbGraphMeeting *meeting = &graph->meeting[i];
	if (meeting->num_strands == 2
	    && meeting->label[0] == -1 && meeting->label[1] == -1 )
	{
	    OrbGraphMeeting *nbr0, *nbr1;

	    nbr0 = &graph->meeting[meeting->neighbor[0]];
	    nbr1 = &graph->meeting[meeting->neighbor[1]];

	    if (nbr0 == meeting && nbr1 == meeting)
	    {
		/* let's put a twist in this component */

		int component = meeting->component[0];

		my_free( meeting->component );
		my_free( meeting->label );
		my_free( meeting->strand );
		my_free( meeting->neighbor );

		meeting->type		= Cross;
		meeting->label 		= NEW_ARRAY( 4, int );
		meeting->component	= NEW_ARRAY( 4, int );
		meeting->neighbor	= NEW_ARRAY( 4, int );
		meeting->strand		= NEW_ARRAY( 4, int );

		meeting->num_strands = 4;

		for ( int j = 0; j< 4; j++)
		{
		    meeting->label[j] = -1;
		    meeting->component[j] = component;
		    meeting->neighbor[j] = i;
		}

		meeting->strand[0] = 3;
		meeting->strand[1] = 2;
		meeting->strand[2] = 1;
		meeting->strand[3] = 0;
		meeting->handedness = 1;
		i++;
	    }
	    else
	    {
		nbr0->strand[meeting->strand[0]] = meeting->strand[1];
		nbr1->strand[meeting->strand[1]] = meeting->strand[0];

		nbr0->neighbor[meeting->strand[0]] = meeting->neighbor[1];
		nbr1->neighbor[meeting->strand[1]] = meeting->neighbor[0];
		my_free(meeting->strand);
		my_free(meeting->neighbor);
		my_free(meeting->component);
		my_free(meeting->label);
		orb_remove_meeting( graph, i);
	    }
	}
	else i++;
    }

    for (int i = 0; i < diagram->num_edges; i++) {
	diagram->edges[i]->num_crossings = 0;
	my_free(diagram->edges[i]->crossings);
	diagram->edges[i]->crossings = NULL;
    }

    return graph;
}

/* Corresponds to DiagramCanvas::outputTriangulation in gui/interface.cpp */
/* Ported from DiagramCanvas::outputTriangulation in interface.cpp */
Triangulation *
orb_triangulate_diagram_complement(
    OrbDiagram *diagram,
    Boolean remove_finite_vertices)
{
    Triangulation * t = NULL;

    OrbGraph * g = orb_diagram_to_graph(diagram);
    if (g != NULL) {
	t = orb_triangulate_graph_complement(g, remove_finite_vertices);
    }
    orb_free_graph(g);
    return t;
}

SNAPPEA_NAMESPACE_END_SCOPE
