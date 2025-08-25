#include "diagram.h"
#include "graph.h"

/* for my_free */
#include "kernel.h"
#include "SnapPea.h"

#include <stdio.h>
#include <stdlib.h>

void initialize_diagram(
    Diagram * diagram)
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

void initialize_diagram_vertex(
    DiagramVertex *vertex)
{
    vertex->vertex_id = 0;
    vertex->link_id = -1;
    vertex->num_incident_end_data = 0;
    vertex->incident_end_data = NULL;
}

static void free_diagram_vertex(
    DiagramVertex *vertex)
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

void free_diagram(
    Diagram *diagram)
{
    int i;

    if (!diagram) {
	return;
    }

    if (diagram->vertices) {
	for (i = 0; i < diagram->num_vertices; ++i) {
	    free_diagram_vertex(diagram->vertices[i]);
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

void initialize_diagram_edge(
    DiagramEdge * edge)
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

void add_end_data_to_diagram_vertex(
    DiagramEndData * data,
    DiagramVertex * vertex)
{
    DiagramEndData ** new_incident_end_data =
	NEW_ARRAY(vertex->num_incident_end_data + 1, DiagramEndData*);

    for (int i = 0; i < vertex->num_incident_end_data; i++) {
	new_incident_end_data[i] = vertex->incident_end_data[i];
    }
    new_incident_end_data[vertex->num_incident_end_data++] = data;
    my_free(vertex->incident_end_data);
    vertex->incident_end_data = new_incident_end_data;
}

void add_crossing_to_diagram_edge(
    DiagramCrossing * crossing,
    DiagramEdge * edge)
{
    DiagramCrossing ** new_crossings =
	NEW_ARRAY(edge->num_crossings + 1, DiagramCrossing*);

    for (int i = 0; i < edge->num_crossings; i++) {
        new_crossings[i] = edge->crossings[i];
    }
    new_crossings[edge->num_crossings++] = crossing;
    my_free(edge->crossings);
    edge->crossings = new_crossings;
}

void assign_diagram_arcs(
    Diagram * diagram)
{
    diagram->num_arcs = 0;

    Boolean drilled_arc;

    int queue_begin = 0;
    int queue_end = 0;

    /* Each edge appears on the queue at most ones */
    DiagramEdge ** queue = NEW_ARRAY(diagram->num_edges, DiagramEdge *);

    Boolean * visited = NEW_ARRAY(diagram->num_edges, Boolean);

    for (int i = 0; i < diagram->num_edges; i++) {
	DiagramEdge * edge = diagram->edges[i];
	edge->edge_id = i;
	edge->arc_id = -1;
	visited[i] = FALSE;
    }

    {
	DiagramEdge * edge = diagram->edges[0];
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
	    DiagramEdge * e = queue[queue_begin++];
	    for (int i = 0; i < 2; i++)
	    {
		DiagramVertex * v = e->vertex[i];

		if (e->vertex[i]->num_incident_end_data == 2)
		{
		    int j = (v->incident_end_data[0]->edge == e) ? 1 : 0;
		    DiagramEdge *e1 = v->incident_end_data[j]->edge;
		    Boolean needsSwitch = v->incident_end_data[j]->type == i;
		    if (!visited[e1->edge_id])
		    {
			visited[e1->edge_id] = TRUE;

			if (!drilled_arc) {
			    e1->arc_id = diagram->num_arcs;
			}

			if (needsSwitch)
			{
			    DiagramVertex *temp;
			    temp = e1->vertex[0];
			    e1->vertex[0] = e1->vertex[1];
			    e1->vertex[1] = temp;

			    for (int j = 0; j < 2; j++)
			    {
				for (int k = 0; k < e1->vertex[j]->num_incident_end_data; k++)
				{
				    if (e1->vertex[j]->incident_end_data[k]->edge == e1) {
					e1->vertex[j]->incident_end_data[k]->type =
					    1 - e1->vertex[j]->incident_end_data[k]->type;
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

	    DiagramEdge * edge = diagram->edges[i];

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

void assign_diagram_links(
    Diagram * diagram)
{
    diagram->num_links = 0;

    int queue_begin = 0;
    int queue_end = 0;

    /* Each edge appears on the queue at most ones */
    DiagramEdge ** queue = NEW_ARRAY(diagram->num_edges, DiagramEdge *);

    Boolean * visited = NEW_ARRAY(diagram->num_edges, Boolean);

    for (int i = 0; i < diagram->num_edges; i++)
    {
	DiagramEdge * edge = diagram->edges[i];
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
	    DiagramEdge * e = queue[queue_begin++];
	    for (int i = 0; i < 2; i++)
	    {
		DiagramVertex * v = e->vertex[i];

		for (int j = 0; j < v->num_incident_end_data; j++)
		{
		    DiagramEdge *e1 = v->incident_end_data[j]->edge;
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
	    DiagramEdge * e = diagram->edges[i];
	    queue[queue_end++] = e;
	    e->link_id = diagram->num_links;
	}
    }

    for (int i = 0; i < diagram->num_vertices; i++)
    {
	DiagramVertex * v = diagram->vertices[i];
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
		DiagramEdge * e = diagram->edges[j];
		if (e->arc_id != arc) {
		    continue;
		}
		for (int k = 0; k < 2; k++) {
		    DiagramVertex * v1 = e->vertex[k];
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

char * dump_diagram(
    Diagram * diagram)
{
    size_t size = 10000000;

    char * buffer = my_malloc(size);
    char * p = buffer;
    char * end = buffer + size - 1;

    p += snprintf(p, end - p, "num_arcs = %d\n", diagram->num_arcs);
    p += snprintf(p, end - p, "num_links = %d\n", diagram->num_links);

    for (int i = 0; i < diagram->num_vertices; i++)
    {
	DiagramVertex * v = diagram->vertices[i];
	p += snprintf(p, end - p, "Vertex:\n");
	p += snprintf(p, end - p, "    %d %d\n", v->x, v->y);
	p += snprintf(p, end - p, "    %d %d %d\n", v->connected_component, v->vertex_id, v->link_id);
	for (int j = 0; j < v->num_incident_end_data; j++) {
	    DiagramEndData * e = v->incident_end_data[j];
	    p += snprintf(p, end - p, "    End data:\n");
	    p += snprintf(p, end - p, "        Edge: %d\n", e->edge->edge_id);
	    p += snprintf(p, end - p, "        %d %d %lf\n",
			  e->type, e->singular, e->angle);
	}
    }

    for (int i = 0; i < diagram->num_edges; i++)
    {
	DiagramEdge * e = diagram->edges[i];
	p += snprintf(p, end - p, "Edge:\n");
	p += snprintf(p, end - p, "    %d %d\n", e->vertex[0]->vertex_id, e->vertex[1]->vertex_id);
	p += snprintf(p, end - p, "    %d %d %d %d\n", e->edge_id, e->arc_id, e->link_id, e->edge_type);
	for (int j = 0; j < e->num_crossings; j++) {
	    DiagramCrossing * c = e->crossings[j];
	    p += snprintf(p, end - p, "    Crossing %d\n", c->crossing_id);
	}
    }

    for (int i = 0; i < diagram->num_crossings; i++)
    {
	DiagramCrossing * c = diagram->crossings[i];
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
point_position_to_complex(
    DiagramVertex * v)
{
    Complex r;
    r.real = v->x;
    r.imag = v->y;
    return r;
}

void
assign_diagram_crossing_signs(
    Diagram * diagram)
{
    for (int i = 0; i < diagram->num_crossings; i++)
    {
        DiagramCrossing *c = diagram->crossings[i];

        Complex z1 = point_position_to_complex(c->over->vertex[diagramBegin]);
        Complex z2 = point_position_to_complex(c->over->vertex[diagramEnd]);
        Complex z3 = point_position_to_complex(c->under->vertex[diagramBegin]);

        /* w = ( z3 - z1 ) / ( z2 - z1 ); */
	Complex w = complex_div(complex_minus(z3, z1), complex_minus(z2, z1));

        c->crossing_sign = ( w.imag > 0 ) ? 1 : -1;
    }
}

void
assign_diagram_end_data_angles(
    Diagram * diagram)
{
    for (int i = 0; i < diagram->num_vertices; i++)
    {
        DiagramVertex * v = diagram->vertices[i];
	for (int j = 0; j < v->num_incident_end_data; j++)
	{
  	    DiagramEndData * e = v->incident_end_data[j];
	    Complex z =
  	        complex_minus(
		    point_position_to_complex(e->edge->vertex[diagramEnd]),
		    point_position_to_complex(e->edge->vertex[diagramBegin]));
	    if (e->type == diagramEnd) {
    	        z = complex_negate(z);
	    }
	    e->angle = -atan2(z.imag, z.real);/* minus???? */
	}
    }
}

void
assign_crossings_to_diagram_edges(
    Diagram * diagram)
{
    for (int i = 0; i < diagram->num_crossings; i++)
    {
        DiagramCrossing * c = diagram->crossings[i];
	add_crossing_to_diagram_edge(c, c->over);
	add_crossing_to_diagram_edge(c, c->under);
    }
}

void
prepare_diagram_components_for_output(
    Diagram * diagram)
{
    {
        int     link_id = -5;

	for (int i = 0; i < diagram->num_arcs; i++)
	{
    	    for (int j = 0; j < diagram->num_edges; j++)
	    {
	        DiagramEdge *e = diagram->edges[j];
		if (e->arc_id == i)
		{
  		    for(int k = 0; k < 2; k++)
		    {
		        if (e->vertex[k]->link_id != -1 )
		        {
			    e->vertex[k]->incident_end_data[get_diagram_strand(e,e->vertex[k])]->singular=(k==1);
			    if (k == 0) {
			        link_id = e->vertex[k]->link_id;
			    }
			}
		    }
		}
	    }
    	    for (int j = 0; j < diagram->num_edges; j++)
	    {
	        DiagramEdge *e = diagram->edges[j];
		if (e->arc_id == i) {
  		    e->link_id = link_id;
		}
	    }
	}
    }

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        DiagramVertex *v = diagram->vertices[i];
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

int get_diagram_strand(DiagramEdge * e, DiagramVertex * v)
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

DiagramCrossing * get_next_diagram_crossing(DiagramEdge *e, DiagramCrossing *c)
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

DiagramCrossing * get_prev_diagram_crossing(DiagramEdge *e, DiagramCrossing *c)
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

static int cmp(double d1, double d2)
{
    if (d1 < d2) {
	return -1;
    }
    if (d1 > d2) {
	return +1;
    }
    return 0;
}

static int ed_more(const void *d1, const void *d2)
{
    /* return ed1->angle > ed2->angle; */

    DiagramEndData * ed1 = *(DiagramEndData **)d1;
    DiagramEndData * ed2 = *(DiagramEndData **)d2;

    return cmp(ed2->angle, ed1->angle);
}

static int crossing_less(const void *d1, const void *d2)
{
    DiagramCrossing *c1 = *(DiagramCrossing **)d1;
    DiagramCrossing *c2 = *(DiagramCrossing **)d2;

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

void remove_meeting( Graph *graph, int index )
{
    if ( index>= graph->num_meetings) {
	return;
    }

    GraphMeeting *old_array = graph->meeting;
    GraphMeeting *new_array = NEW_ARRAY(graph->num_meetings-1, GraphMeeting );

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

/* Ported from DiagramCanvas::outputTriangulation in interface.cpp */
Graph * diagram_to_graph(
    Diagram * diagram)
{
    int num_meetings = 0;

    assign_diagram_crossing_signs(diagram);
    assign_diagram_end_data_angles(diagram);
    assign_crossings_to_diagram_edges(diagram);
    prepare_diagram_components_for_output(diagram);

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        diagram->vertices[i]->vertex_id = num_meetings++;
	qsort(diagram->vertices[i]->incident_end_data,
	      diagram->vertices[i]->num_incident_end_data,
	      sizeof(DiagramEndData*),
	      &ed_more);
    }

    for(int i = 0; i < diagram->num_crossings; i++) {
        diagram->crossings[i]->crossing_id = num_meetings++;
    }

    for(int i = 0; i < diagram->num_edges; i++) {
        qsort(diagram->edges[i]->crossings,
	      diagram->edges[i]->num_crossings,
	      sizeof(DiagramCrossing*),
	      &crossing_less);
    }

    Graph *graph = NEW_STRUCT(Graph);
    graph->num_meetings = num_meetings;
    graph->num_components = diagram->num_links;
    graph->num_free_loops = 0;
    graph->meeting = NEW_ARRAY( num_meetings, GraphMeeting );

    for(int i = 0; i < diagram->num_vertices; i++)
    {
        GraphMeeting *meeting = &graph->meeting[i];
	DiagramVertex *vertex = diagram->vertices[i];

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
	    DiagramEndData *ed = vertex->incident_end_data[j];
	    DiagramEdge *e = ed->edge;
	    meeting->label[j] = (ed->singular) ? e->arc_id : -1;

	    meeting->component[j] = diagram->vertices[i]->link_id;

	    if (ed->type == diagramBegin)
	    {
	        if (e->num_crossings == 0)
		{
		    meeting->strand[j] = get_diagram_strand(e, e->vertex[diagramEnd]);
		    meeting->neighbor[j] = e->vertex[diagramEnd]->vertex_id;
		}
		else
		{
  		    DiagramCrossing *c = e->crossings[0];
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
		    meeting->strand[j] = get_diagram_strand(e, e->vertex[diagramBegin]);
		    meeting->neighbor[j] = e->vertex[diagramBegin]->vertex_id;
		}
		else
		{
		    DiagramCrossing *c = e->crossings[e->num_crossings-1];
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
	GraphMeeting *meeting = &graph->meeting[diagram->num_vertices+i];
	DiagramCrossing *c = diagram->crossings[i], *other;

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

	DiagramEdge * e = c->over;

	meeting->component[0] = e->link_id;
	meeting->component[2] = e->link_id;

	if ( (other = get_prev_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[0] = get_diagram_strand( e, e->vertex[diagramBegin] );
	    meeting->neighbor[0] = e->vertex[diagramBegin]->vertex_id;
	}
	else
	{
	    meeting->neighbor[0] = other->crossing_id;
	    meeting->strand[0] = (other->over==e) ? 2 : 3;
	}

	if ( (other = get_next_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[2] = get_diagram_strand( e, e->vertex[diagramEnd] );
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

	if ( (other = get_prev_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[1] = get_diagram_strand( e, e->vertex[diagramBegin] );
	    meeting->neighbor[1] = e->vertex[diagramBegin]->vertex_id;
	}
	else
	{
	    meeting->neighbor[1] = other->crossing_id;
	    meeting->strand[1] = (other->over==e) ? 2 : 3;
	}

	if ( (other = get_next_diagram_crossing( e, c)) == NULL)
	{
	    meeting->strand[3] = get_diagram_strand( e, e->vertex[diagramEnd] );
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
	    GraphMeeting *meeting, *nbr1, *nbr3;

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
	GraphMeeting *meeting = &graph->meeting[i];
	if (meeting->num_strands == 2
	    && meeting->label[0] == -1 && meeting->label[1] == -1 )
	{
	    GraphMeeting *nbr0, *nbr1;

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
		remove_meeting( graph, i);
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

/* Ported from DiagramCanvas::outputTriangulation in interface.cpp */
Triangulation *
triangulate_diagram_complement(
    Diagram *diagram,
    Boolean remove_vertices)
{
    Triangulation * t = NULL;

    Graph * g = diagram_to_graph(diagram);
    if (g != NULL) {
	t = triangulate_graph_complement(g, remove_vertices);
    }
    free_graph(g);
    return t;
}
