#include "diagram_io.h"

#include "diagram.h"
#include "ostream.h"
#include "kernel.h"

#include <stdio.h>

static Boolean fill_diagram(
    OrbDiagram * diagram, char *file_data)
{
    int chars_consumed;

    int num_vertices;
    if (sscanf(file_data, "%d%n", &num_vertices, &chars_consumed) != 1) {
	return FALSE;
    }
    file_data += chars_consumed;

    diagram->vertices = NEW_ARRAY( num_vertices, OrbDiagramVertex* );

    for (diagram->num_vertices = 0;
	 diagram->num_vertices < num_vertices;
	 diagram->num_vertices++) {
	OrbDiagramVertex * v = NEW_STRUCT( OrbDiagramVertex );
	diagram->vertices[diagram->num_vertices] = v;
        orb_initialize_diagram_vertex(v);
	int index;
	if (sscanf(
		file_data,
		"%d%d%d%n",
		&index, &v->x, &v->y,
		&chars_consumed) != 3) {
	    return FALSE;
	}
	file_data += chars_consumed;
    }

    int num_edges;
    if (sscanf(file_data, "%d%n", &num_edges, &chars_consumed) != 1) {
	return FALSE;
    }
    file_data += chars_consumed;

    diagram->edges = NEW_ARRAY( num_edges, OrbDiagramEdge* );

    for (diagram->num_edges = 0;
	 diagram->num_edges < num_edges;
	 diagram->num_edges++) {
	OrbDiagramEdge * e = NEW_STRUCT( OrbDiagramEdge );
	diagram->edges[diagram->num_edges] = e;
	orb_initialize_diagram_edge(e);

	int index, vertex_id0, vertex_id1, edge_type;
	if (sscanf(file_data,
		   "%d%d%d%d%n",
		   &index,
		   &vertex_id0,
		   &vertex_id1,
		   &edge_type,
		   &chars_consumed) != 4) {
	    return FALSE;
	}

	e->vertex[diagramBegin] = diagram->vertices[vertex_id0];
	e->vertex[diagramEnd]   = diagram->vertices[vertex_id1];
	e->edge_type = edge_type;

	OrbDiagramEndData * begin_data = NEW_STRUCT(OrbDiagramEndData);
	begin_data->edge = e;
	begin_data->type = diagramBegin;
	begin_data->singular = FALSE;
	begin_data->angle = 0.0;
	orb_add_end_data_to_diagram_vertex(begin_data, diagram->vertices[vertex_id0]);

	OrbDiagramEndData * end_data = NEW_STRUCT(OrbDiagramEndData);
	end_data->edge = e;
	end_data->type = diagramEnd;
	end_data->singular = FALSE;
	end_data->angle = 0.0;
	orb_add_end_data_to_diagram_vertex(end_data, diagram->vertices[vertex_id1]);

	file_data += chars_consumed;
    }

    int num_crossings;
    if (sscanf(file_data, "%d%n", &num_crossings, &chars_consumed) != 1) {
	return FALSE;
    }
    file_data += chars_consumed;

    diagram->crossings = NEW_ARRAY( num_crossings, OrbDiagramCrossing* );

    for (diagram->num_crossings = 0;
	 diagram->num_crossings < num_crossings;
	 diagram->num_crossings++) {
	OrbDiagramCrossing * c = NEW_STRUCT( OrbDiagramCrossing );
	diagram->crossings[diagram->num_crossings] = c;

	int index, edge_id0, edge_id1;
	if (sscanf(file_data,
		   "%d%d%d%d%d%lf%lf%n",
		   &index,
		   &c->x, &c->y,
		   &edge_id0, &edge_id1,
		   &c->position_on_overstrand, &c->position_on_understrand,
		   &chars_consumed) != 7) {
	    return FALSE;
	}

	c->over  = diagram->edges[edge_id0];
	c->under = diagram->edges[edge_id1];

	file_data += chars_consumed;
    }

    return TRUE;
}

OrbDiagram *orb_read_diagram_from_string(
    char * file_data)
{
    OrbDiagram *diagram = NEW_STRUCT(OrbDiagram);
    orb_initialize_diagram(diagram);
    if (fill_diagram(diagram, file_data)) {
	orb_assign_diagram_arcs(diagram);
	orb_assign_diagram_links(diagram);
	return diagram;
    }

    orb_free_diagram(diagram);
    return NULL;
}

/* Ported from DiagramCanvas::saveDiagram in interface.cpp */
void
orb_write_diagram_to_stream(
    OStream * stream,
    OrbDiagram * diagram)
{
    ostream_printf(stream, "%d\n", diagram->num_vertices);
    for (int i = 0; i < diagram->num_vertices; i++)
    {
	ostream_printf(
	    stream,
	    "%d\t%d\t%d\n",
	    i, diagram->vertices[i]->x, diagram->vertices[i]->y);
	diagram->vertices[i]->vertex_id = i;
    }

    ostream_printf(stream, "\n%d\n", diagram->num_edges);
    for (int i = 0; i < diagram->num_edges; i++)
    {
	ostream_printf(
	    stream,
	    "%d\t%d\t%d\t%d\n",
	    i,
	    diagram->edges[i]->vertex[diagramBegin]->vertex_id,
	    diagram->edges[i]->vertex[diagramEnd]->vertex_id,
	    (int)diagram->edges[i]->edge_type);
	diagram->edges[i]->edge_id = i;
    }

    ostream_printf(stream, "\n%d\n", diagram->num_crossings);
    for (int i = 0; i < diagram->num_crossings; i++)
    {
	ostream_printf(
	    stream,
	    "%d\t%d\t%d,\t%d\t%d\t%lf\t%lf\n",
	    i,
	    diagram->crossings[i]->x,
	    diagram->crossings[i]->y,
	    diagram->crossings[i]->over->edge_id,
	    diagram->crossings[i]->under->edge_id,
	    diagram->crossings[i]->position_on_overstrand,
	    diagram->crossings[i]->position_on_understrand);
    }
}

char *
orb_write_diagram_to_string(
    OrbDiagram * diagram)
{
    OStream stream;
    string_stream_init(&stream);

    orb_write_diagram_to_stream(&stream, diagram);

    return stream.buffer;
}
