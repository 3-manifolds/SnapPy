#include "diagram_canvas.h"




//#include <qpoint.h>
//#include <qstring.h>
//#include <qfiledialog.h>
//#include <qmessagebox.h>
#include <complex>
#include <cstdio>
#include <set>
#include <algorithm>
#include <queue>

extern std::complex<double> point_to_complex( QPoint p );
extern bool ed_more( EndData *ed1, EndData *ed2 );
extern bool crossing_less( Crossing *c1, Crossing *c2 );
extern int  get_strand( Edge *e, Vertex *v );
extern Crossing *get_prev_crossing( Edge *e, Crossing *c);
extern Crossing *get_next_crossing( Edge *e, Crossing *c);
extern void remove_meeting( Graph *graph, int i );
extern "C" Triangulation *triangulate_graph_complement( Graph *gamma, int *num_singular_edges );
extern "C" void dump_triangulation( Triangulation *manifold, Boolean is_cover );


void DiagramCanvas::highlightEdge( int e )
{

    std::vector<Edge *>::iterator edge_it;

    edge_it = edgeList.begin();

    while ( edge_it < edgeList.end() )
    {
    	if ( (*edge_it)->arc_id == e)
	{
		(*edge_it)->selected = TRUE;
		(*edge_it)->thickness = THICK;
	}
    	else
	{
		(*edge_it)->selected = FALSE;
		(*edge_it)->thickness = THIN;
	}

	edge_it++;
    }

//     buffer.fill( CANVAS );
//     show();
//     update();
}

void DiagramCanvas::exportSlot()
{
	Triangulation *manifold = NULL;

	if (untouched)
	{
            //	QMessageBox::information(this,"New Diagram",
//				"The canvas is empty.");	
		return;
	}

	if (edge_creation_in_progress)
	{
            //	QMessageBox::information(this,"New Diagram",
            //		"Finish the edge first.");
		return;
	}

	if (vertex_drag_in_progress)
	{
//		QMessageBox::information(this,"New Diagram",
            //			"Finish vertex placement first.");
		return;
	}

	if (!isReadyForTriangulation())
	{
//		QMessageBox::information(this,"New Diagram",
            //			"The current diagram is not ready for\n"
            //		"the triangulation process.");
		return;
	}

//	assign_cuffs();
	manifold = outputTriangulation();

//	emit exporting(manifold);
}

void DiagramCanvas::clearDiagram()
{
    int i;

    vertex_drag_in_progress = FALSE;

    if (edge_creation_in_progress)
    {
        //       paint( currentEdge, TRUE );
         currentEdge->vertex[end]->~Vertex();
         currentEdge->~Edge();
         edge_creation_in_progress = FALSE;
    }
  
    //  for ( i=0; i< (int) edgeList.size(); ++i ) paint( edgeList[i], TRUE );
//    update();

    for(i=0; i< (int) vertexList.size(); i++)
		vertexList[i]->~Vertex();
    vertexList.clear();

    for(i=0; i< (int) edgeList.size(); i++)                                      
               edgeList[i]->~Edge(); 
    edgeList.clear();

    for(i=0; i< (int) crossingList.size(); i++)                                      
                crossingList[i]->~Crossing(); 
    crossingList.clear();

    untouched = TRUE;
    num_arcs = 0;
    num_links = 0;

//    buffer.resize( size() );
//    buffer.fill( CANVAS );
//    update();
}

void DiagramCanvas::saveDiagram( QTextStream &stream )
{
    int i, num_edges, num_crossings, num_vertices;

    num_vertices = vertexList.size();
    num_edges = edgeList.size();
    num_crossings = crossingList.size();

    stream << "\n" << num_vertices << "\n";
    for ( i=0; i<num_vertices; ++i )
    {
          vertexList[i]->vertex_id = i;
          stream << i << "\t";
	  stream <<  vertexList[i]->position.x() << "\t";
	  stream << vertexList[i]->position.y() << "\n";
    }

    stream << "\n" << num_edges << "\n";

    for ( i=0; i<num_edges; ++i )
    {
          edgeList[i]->edge_id = i;
	  stream << i << "\t";
	  stream << edgeList[i]->vertex[begin]->vertex_id << "\t";
	  stream << edgeList[i]->vertex[end]->vertex_id << "\t";
	  stream << ((edgeList[i]->edge_type == drilled) ? 1 : 0) << "\n";
    }

    stream << "\n" << num_crossings << "\n";

    for (i=0; i<num_crossings; ++i )
    {
        Crossing *c = crossingList[i];

        c->crossing_id = i;
	stream << i << "\t";
	stream << c->position.x() << "\t";
	stream << c->position.y() << "\t";
	stream << c->over->edge_id << "\t";
	stream << c->under->edge_id << "\t";
	stream << c->position_on_overstrand << "\t";
	stream << c->position_on_understrand << "\n";
    }
}

void DiagramCanvas::readDiagram( QTextStream &stream )
{
    int i, num_vertices, num_edges, num_crossings;

    clearDiagram();

    vertexList.clear();
  
    stream >> num_vertices;

    for (i=0;i<num_vertices;++i)
    {
            int index, x, y;

            stream >> index >> x >> y;
            Vertex *v = new Vertex( QPoint(x, y) );
            vertexList.insert( vertexList.end(), v );
    }

    edgeList.clear();
 
    stream >> num_edges;
    for( i = 0; i<num_edges;++i)
    {
           int index, begin_id, end_id, type;

           stream >> index >> begin_id >> end_id >> type;
           Edge *e = new Edge( vertexList[begin_id], vertexList[end_id] );
           e->edge_type = (type==1) ? drilled : singular;
           edgeList.insert( edgeList.end(), e );
    }  

    crossingList.clear(); 
   
    stream >> num_crossings;
   
    for ( i=0; i<num_crossings; ++i )
    {
          int index, x, y, over, under;
          double t1, t2;

          stream >> index >> x >> y >> over >> under >> t1 >> t2;
          Crossing *c = new Crossing;
          c->position = QPoint( x, y); 
          c->over = edgeList[over];
          c->under = edgeList[under];
          c->position_on_overstrand  = t1;
          c->position_on_understrand = t2;

          crossingList.insert( crossingList.end(), c );
    }
   
    assign_arcs();
    assign_links();
//    show();
//    update();

    untouched =  (num_arcs + num_links > 0) ? FALSE :TRUE;

}

void DiagramCanvas::deleteSelectedEdges()
{
    std::vector<Edge *>::iterator edge_it;
    std::vector<EndData *>::iterator ed_it;
    std::vector<Crossing *>::iterator crossing_it;
    Vertex *vertex1, *vertex2;

    edge_it = edgeList.begin();

    while ( edge_it < edgeList.end() )
    if ((*edge_it)->selected)
    {
        vertex1 = (*edge_it)->vertex[begin];
        vertex2 = (*edge_it)->vertex[end];

        // paint( *edge_it, TRUE );

        crossing_it = crossingList.begin();

        while ( crossing_it < crossingList.end() )
        {
               Crossing *c = *crossing_it;

               if ( c->over == *edge_it || c->under == *edge_it ) 
               {
                      crossingList.erase( crossing_it );
                      //     update();
               }
               else ++crossing_it;
        }

        if (vertex1->incidentEndData.size() == 1)
        {
            std::vector< Vertex *>::iterator vertex_it;
 
            vertex_it = vertexList.begin();
            while (vertex_it < vertexList.end())
            if (*vertex_it == vertex1)
                vertexList.erase( vertex_it );
            else
                vertex_it++;

            vertex1->~Vertex();
        }
        else
        {

            ed_it = vertex1->incidentEndData.begin();

            while( ed_it < vertex1->incidentEndData.end())
            if ( (*ed_it)->edge == (*edge_it) )
            {
                   (*ed_it)->~EndData();
                   vertex1->incidentEndData.erase( ed_it );
            }
            else ++ed_it;
        }

        if (vertex2->incidentEndData.size() == 1)
        {
            std::vector< Vertex *>::iterator vertex_it;
         
            vertex_it = vertexList.begin();
            while (vertex_it < vertexList.end())
            if (*vertex_it == vertex2)
                vertexList.erase( vertex_it );
            else
                vertex_it++;

            vertex2->~Vertex();
        }
        else 
        { 

            ed_it = vertex2->incidentEndData.begin();

            while( ed_it < vertex2->incidentEndData.end())
            if ( (*ed_it)->edge == (*edge_it) )
            {
                   (*ed_it)->~EndData();
                   vertex2->incidentEndData.erase( ed_it );
            }
            else ++ed_it;
        }


        (*edge_it)->~Edge();
        edgeList.erase(edge_it);
    }
    else ++edge_it;

    if (edgeList.size()==0)
       untouched = TRUE;
    else
    {
       assign_arcs();
       assign_links();
    }

//     buffer.fill( CANVAS );
//     show();
//     update();
}

void DiagramCanvas::deleteArc( int a )
{
    std::vector<Edge *>::iterator edge_it;
    std::vector<EndData *>::iterator ed_it;
    std::vector<Crossing *>::iterator crossing_it;
    Vertex *vertex1, *vertex2;


    edge_it = edgeList.begin();

    while ( edge_it < edgeList.end() )
    if ((*edge_it)->arc_id == a && (*edge_it)->edge_type == singular )
    {
        vertex1 = (*edge_it)->vertex[begin];
        vertex2 = (*edge_it)->vertex[end];

//        paint( *edge_it, TRUE );

        crossing_it = crossingList.begin();

        while ( crossing_it < crossingList.end() )
        {
               Crossing *c = *crossing_it;

               if ( c->over == *edge_it || c->under == *edge_it ) 
               {
                      crossingList.erase( crossing_it );
                      //                    update();
               }
               else ++crossing_it;
        }

        if (vertex1->incidentEndData.size() == 1)
        {
            std::vector< Vertex *>::iterator vertex_it;
 
            vertex_it = vertexList.begin();
            while (vertex_it < vertexList.end())
            if (*vertex_it == vertex1)
                vertexList.erase( vertex_it );
            else
                vertex_it++;

            vertex1->~Vertex();
        }
        else
        {

            ed_it = vertex1->incidentEndData.begin();

            while( ed_it < vertex1->incidentEndData.end())
            if ( (*ed_it)->edge == (*edge_it) )
            {
                   (*ed_it)->~EndData();
                   vertex1->incidentEndData.erase( ed_it );
            }
            else ++ed_it;
        }

        if (vertex2->incidentEndData.size() == 1)
        {
            std::vector< Vertex *>::iterator vertex_it;
         
            vertex_it = vertexList.begin();
            while (vertex_it < vertexList.end())
            if (*vertex_it == vertex2)
                vertexList.erase( vertex_it );
            else
                vertex_it++;

            vertex2->~Vertex();
        }
        else 
        { 
            ed_it = vertex2->incidentEndData.begin();

            while( ed_it < vertex2->incidentEndData.end())
            if ( (*ed_it)->edge == (*edge_it) )
            {
                   (*ed_it)->~EndData();
                   vertex2->incidentEndData.erase( ed_it );
            }
            else ++ed_it;
        }


        (*edge_it)->~Edge();
        edgeList.erase(edge_it);
    }
    else ++edge_it;

    num_arcs--;

    if (edgeList.size()==0)
       untouched = TRUE;
    else
    {
       reassign_arcs( a );
       reassign_links( FALSE );
    }

//     buffer.fill( CANVAS );
//     show();
//     update();
}

void DiagramCanvas::drillArc( int a )
{
    int i;

    for( i=0; i<edgeList.size(); i++)
    if (edgeList[i]->edge_type == singular && edgeList[i]->arc_id == a)
    {
        edgeList[i]->edge_type = drilled;
	edgeList[i]->link_id   = -255;
    }

    num_arcs--;
//    buffer.fill( CANVAS );
    reassign_arcs( a );
    reassign_links( TRUE );
//    update();
}

void DiagramCanvas::reassign_arcs( int a )
{
    int i,j, m, k;
    bool go_again = TRUE;

    for( i=0; i<edgeList.size(); i++)
    if (edgeList[i]->edge_type == singular && edgeList[i]->arc_id > a)
        edgeList[i]->arc_id--;

    while (go_again)
    {
	go_again = FALSE;
    	for(i=0;i<vertexList.size();i++)
    	if (vertexList[i]->incidentEndData.size()==2)
    	{
		if (	vertexList[i]->incidentEndData[0]->edge->edge_type == singular &&
			vertexList[i]->incidentEndData[1]->edge->edge_type == singular &&	
			vertexList[i]->incidentEndData[0]->edge->arc_id    !=
			vertexList[i]->incidentEndData[1]->edge->arc_id )
		{
			int reset_id = MIN( vertexList[i]->incidentEndData[0]->edge->arc_id,
						vertexList[i]->incidentEndData[1]->edge->arc_id);
			int dead_id =  MAX( vertexList[i]->incidentEndData[0]->edge->arc_id,
						vertexList[i]->incidentEndData[1]->edge->arc_id);

			bool needs_flip = FALSE;
			if (vertexList[i]->incidentEndData[0]->type == vertexList[i]->incidentEndData[1]->type)
				/* we need to flip the arrows on arc labelled dead_id */
				needs_flip = TRUE;

			for(j=0;j<edgeList.size();j++)
			if (edgeList[j]->edge_type == singular && edgeList[j]->arc_id == dead_id )
			{
				if (needs_flip)
				{
                                	Vertex *temp = edgeList[j]->vertex[0];
                                	edgeList[j]->vertex[0] = edgeList[j]->vertex[1];
                                	edgeList[j]->vertex[1] = temp;

                                	for(m=0;m<2;m++)
                                		for(k=0;k<edgeList[j]->vertex[m]->incidentEndData.size(); k++ )
                                		if (edgeList[j]->vertex[m]->incidentEndData[k]->edge == edgeList[j] )
                                			edgeList[j]->vertex[m]->incidentEndData[k]->type =
                                				(edgeList[j]->vertex[m]->incidentEndData[k]->type==begin) ?
										end : begin;

					std::vector<Crossing *>::iterator c_it;
					c_it = crossingList.begin();
					while ( c_it < crossingList.end() )
					{
						Crossing *c = *c_it;
							
						if ( c->over == edgeList[j] || c->under == edgeList[j] )
						{
							refresh_crossing( c );
                                                        //				update();
						}

						++c_it;
					}
				}

				edgeList[j]->arc_id = reset_id;
			}

    			for( j=0; j<edgeList.size(); j++)
    			if (edgeList[j]->edge_type == singular && edgeList[j]->arc_id > dead_id )
       				edgeList[j]->arc_id--;

			num_arcs--;
			go_again = TRUE;
		}
   	 }
    }

}

void DiagramCanvas::reassign_links( bool arc_was_drilled )
{
	int i, j, k, id[2], link_id, dead_id;

	if (arc_was_drilled)
	{
		k = 0;
		for(i=0;i<vertexList.size();i++)
//		if (vertexList[i]->incidentEndData.size()>2)
			for(j=0;j<vertexList[i]->incidentEndData.size();j++)
			if (	vertexList[i]->incidentEndData[j]->edge->edge_type == drilled &&
				vertexList[i]->incidentEndData[j]->edge->link_id == -255 &&
				vertexList[i]->link_id > -1 )
				id[k++] = vertexList[i]->link_id;

		if (k!=2)
			uFatalError("reassign_links","interface");
	
		link_id = MIN(id[0], id[1]);
		dead_id = MAX(id[0], id[1] );	

		for(i=0;i<edgeList.size();i++)
		if (edgeList[i]->edge_type == drilled && (edgeList[i]->link_id == -255 || edgeList[i]->link_id == dead_id) )
		{
			edgeList[i]->link_id = link_id;
			edgeList[i]->vertex[0]->link_id = link_id;
			edgeList[i]->vertex[1]->link_id = link_id;
		}

		if (link_id != dead_id )
		{
			for(i=0;i<edgeList.size();i++)
			if (edgeList[i]->edge_type == drilled && edgeList[i]->link_id > dead_id )
				edgeList[i]->link_id--;

			for(i=0;i<vertexList.size();i++)
			if (vertexList[i]->link_id > dead_id )
				vertexList[i]->link_id--;
		}
	}
	else
	{
		/* need to remove an bi-valent cusps if they don't lie on singular loops.
			otherwise we need to select the lowest link_id on the loop to keep */

		for(i = 0; i<num_arcs;i++)
		{
			bool arc_is_loop = TRUE;

			for(j = 0; j < edgeList.size() && arc_is_loop; j++)
			if (	edgeList[j]->edge_type == singular &&
				edgeList[j]->arc_id == i )
			{
				if (	edgeList[j]->vertex[begin]->incidentEndData.size() != 2 ||
					edgeList[j]->vertex[end]->incidentEndData.size() != 2 )
					arc_is_loop = FALSE;
			}

			if (!arc_is_loop)
			for(j = 0; j < edgeList.size(); j++)
			if (	edgeList[j]->edge_type == singular &&
				edgeList[j]->arc_id == i )
			{
				if (	edgeList[j]->vertex[begin]->incidentEndData.size() == 2 &&
					edgeList[j]->vertex[begin]->link_id > -1 )
				{
					dead_id = edgeList[j]->vertex[begin]->link_id;
					edgeList[j]->vertex[begin]->link_id = -1;

					for(k = 0; k < vertexList.size();k++)
					if (	vertexList[k]->link_id > dead_id )
						vertexList[k]->link_id--;

					for(k = 0; k < edgeList.size();k++)
					if (	edgeList[k]->edge_type == drilled &&
						edgeList[k]->link_id > dead_id )
						edgeList[k]->link_id--;
				}

				if (	edgeList[j]->vertex[end]->incidentEndData.size() == 2 &&
					edgeList[j]->vertex[end]->link_id > -1 )
				{
                                        dead_id = edgeList[j]->vertex[end]->link_id;
                                        edgeList[j]->vertex[end]->link_id = -1;

                                        for(k = 0; k < vertexList.size();k++)
                                        if (    vertexList[k]->link_id > dead_id )
                                                vertexList[k]->link_id--;

                                        for(k = 0; k < edgeList.size();k++)
                                        if (    edgeList[k]->edge_type == drilled &&
                                                edgeList[k]->link_id > dead_id )
                                                edgeList[k]->link_id--;
				}
			}
		}

		/* the remain bi-valent  will lie on singular loops. */
		/* keep the lowest link_id on each loop */

		for( i = 0; i < num_arcs; i++)
		{
			Vertex *curV = NULL;

			for( j = 0; j < vertexList.size(); j++)
			if (	vertexList[i]->incidentEndData.size() == 2 &&
				vertexList[i]->link_id > -1 &&
				vertexList[i]->incidentEndData[0]->edge->arc_id == i )
			{
				if (curV==NULL)
					curV = vertexList[i];
				else
				{
					if (vertexList[i]->link_id < curV->link_id )
					{
						dead_id = curV->link_id;
						curV->link_id = -1;
						curV = vertexList[i];
					}
					else
					{
						dead_id = vertexList[i]->link_id;
						vertexList[i]->link_id = -1;
					}

                                        for(k = 0; k < vertexList.size();k++)
                                        if (    vertexList[k]->link_id > dead_id )
                                                vertexList[k]->link_id--;

                                        for(k = 0; k < edgeList.size();k++)
                                        if (    edgeList[k]->edge_type == drilled &&
                                                edgeList[k]->link_id > dead_id )
                                                edgeList[k]->link_id--;
				}
			}
		}
	}
}

void DiagramCanvas::getCrossingSigns()
{
    int i, num_crossings;
    std::complex<double> z1, z2, z3, z4, w;
   
    num_crossings = crossingList.size();
   
    for ( i=0; i<num_crossings; ++i )
    {
        Crossing *c = crossingList[i];

        z1 = point_to_complex( c->over->vertex[begin]->position );
        z2 = point_to_complex( c->over->vertex[end]->position );
        z3 = point_to_complex( c->under->vertex[begin]->position );
        z4 = point_to_complex( c->under->vertex[end]->position );

        w = ( z3 - z1 ) / ( z2 - z1 );

        c->crossing_sign = ( w.imag() > 0 ) ? 1 : -1;
    }
}

void DiagramCanvas::ed_angles()
{
    int i;
    std::vector<EndData *>::iterator ed_it;
    std::complex<double> z;

    for(i=0; i<vertexList.size(); i++)
    for( ed_it = vertexList[i]->incidentEndData.begin();
         ed_it != vertexList[i]->incidentEndData.end();
         ed_it++ )
    { 
        if ((*ed_it)->type == begin)
             z = point_to_complex( (*ed_it)->edge->vertex[end]->position
                          - (*ed_it)->edge->vertex[begin]->position );
        else
             z = point_to_complex( (*ed_it)->edge->vertex[begin]->position
                          - (*ed_it)->edge->vertex[end]->position );

        (*ed_it)->angle = -arg(z);/* minus???? */
    }
}

void DiagramCanvas::assign_crossings_to_edges()
{
    int i;

    for(i = 0; i < crossingList.size(); i++)
    {
           crossingList[i]->over->crossings.insert(
                       crossingList[i]->over->crossings.end(), crossingList[i] );
           crossingList[i]->under->crossings.insert(
                       crossingList[i]->under->crossings.end(), crossingList[i] ); 
    }
}


Triangulation *DiagramCanvas::outputTriangulation()
{
    int num_meetings, i, num_singular_edges;

    if (edge_creation_in_progress || vertex_drag_in_progress )
        fprintf( stderr, "Err: diagram is still being modified.\n");

    getCrossingSigns();
    ed_angles();
    assign_crossings_to_edges();
    prepare_components_for_output();

    num_meetings = 0;

    for(i=0;i<vertexList.size();i++)
    {
          vertexList[i]->vertex_id = num_meetings++;
          sort( vertexList[i]->incidentEndData.begin(),
                vertexList[i]->incidentEndData.end(), ed_more );
    }

    for(i=0;i<crossingList.size();i++)
          crossingList[i]->crossing_id = num_meetings++;

    for(i=0;i<edgeList.size();i++)
          sort( edgeList[i]->crossings.begin(),
                edgeList[i]->crossings.end(), crossing_less );

    Graph *graph = new Graph;
    graph->num_meetings = num_meetings;
    graph->num_components = num_links;
    graph->num_free_loops = 0;
    graph->meeting = NEW_ARRAY( num_meetings, GraphMeeting );
   
    for(i=0;i<vertexList.size();i++)
    {
          GraphMeeting *meeting = &graph->meeting[i];
          Vertex *vertex = vertexList[i];
          int j;

          meeting->type = Inter;
          meeting->num_strands = vertex->incidentEndData.size();
          meeting->strand = NEW_ARRAY(meeting->num_strands,int);
          meeting->neighbor = NEW_ARRAY(meeting->num_strands,int);
          meeting->component = NEW_ARRAY(meeting->num_strands,int);
          meeting->label = NEW_ARRAY(meeting->num_strands,int);
          meeting->tet = NULL;
	  meeting->handedness = 0;

          for(j=0;j<meeting->num_strands;j++)
          {
                EndData *ed = vertex->incidentEndData[j];
                Edge *e = ed->edge;
                meeting->label[j] = (ed->singular) ? e->arc_id : -1;

                meeting->component[j] = vertexList[i]->link_id;

                if (ed->type == begin)
                {
                      if (e->crossings.empty())
                      {
                              meeting->strand[j] = get_strand(e, e->vertex[end]);
                              meeting->neighbor[j] = e->vertex[end]->vertex_id;
                      }
                      else
                      {
                              Crossing *c = *(e->crossings.begin());
                              meeting->neighbor[j] = c->crossing_id;

                              if (c->over == e)
                                      meeting->strand[j] = 0;
                              else    meeting->strand[j] = 1;
                              
                      } 
                }
                else
                {
                     if (e->crossings.empty())
                      {
                              meeting->strand[j] = get_strand(e, e->vertex[begin]);
                              meeting->neighbor[j] = e->vertex[begin]->vertex_id;
                      }
                      else
                      {
                              Crossing *c = *(e->crossings.end()-1);
                              meeting->neighbor[j] = c->crossing_id;

                              if (c->over == e)
                                      meeting->strand[j] = 2;
                              else    meeting->strand[j] = 3;

                      }
                }
          }
    } 

    for( i=0;i<crossingList.size();i++)
    {
          GraphMeeting *meeting = &graph->meeting[vertexList.size()+i];
          Crossing *c = crossingList[i], *other;
          Edge *e;
          int j;

          meeting->type = Cross;
          meeting->num_strands = 4;
          meeting->strand = NEW_ARRAY( 4, int );
          meeting->neighbor = NEW_ARRAY( 4, int ); 
          meeting->component = NEW_ARRAY( 4, int ); 
          meeting->label = NEW_ARRAY( 4, int ); 
          meeting->tet = NULL;
	  meeting->handedness = c->crossing_sign;

          for (j=0;j<4;j++)
                meeting->label[j] = -1;

          e = c->over;

          meeting->component[0] = e->link_id;
          meeting->component[2] = e->link_id;

          if ( (other = get_prev_crossing( e, c)) == NULL)
          {
                meeting->strand[0] = get_strand( e, e->vertex[begin] );
                meeting->neighbor[0] = e->vertex[begin]->vertex_id;
          }
          else
          {
                meeting->neighbor[0] = other->crossing_id;
                meeting->strand[0] = (other->over==e) ? 2 : 3;
          }

          if ( (other = get_next_crossing( e, c)) == NULL)
          {
                meeting->strand[2] = get_strand( e, e->vertex[end] );
                meeting->neighbor[2] = e->vertex[end]->vertex_id;
          }
          else
          {
                meeting->neighbor[2] = other->crossing_id;
                meeting->strand[2] = (other->over==e) ? 0 : 1;
          }

          e = c->under;

          meeting->component[1] = e->link_id;
          meeting->component[3] = e->link_id;

          if ( (other = get_prev_crossing( e, c)) == NULL)
          {
                meeting->strand[1] = get_strand( e, e->vertex[begin] );
                meeting->neighbor[1] = e->vertex[begin]->vertex_id;
          }
          else
          {
                meeting->neighbor[1] = other->crossing_id;
                meeting->strand[1] = (other->over==e) ? 2 : 3;
          }

          if ( (other = get_next_crossing( e, c)) == NULL)
          {
                meeting->strand[3] = get_strand( e, e->vertex[end] );
                meeting->neighbor[3] = e->vertex[end]->vertex_id;
          }
          else
          {
                meeting->neighbor[3] = other->crossing_id;
                meeting->strand[3] = (other->over==e) ? 0 : 1;
          }

    }

    for(i=0;i<crossingList.size();i++)
    if (crossingList[i]->crossing_sign == 1)
    {
          int temp_strand , temp_neighbor;
          GraphMeeting *meeting, *nbr1, *nbr3;
        
          meeting = &graph->meeting[crossingList[i]->crossing_id];

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

    for(i=0;i<edgeList.size();i++)
         edgeList[i]->crossings.clear();

    i =0;
    while (i<graph->num_meetings)
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

        		for ( int j=0; j< 4; j++)
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


   Triangulation *manifold = triangulate_graph_complement( graph, FALSE );

   return manifold;

}


void DiagramCanvas::invert_arc( Edge *e )
{
     int i;

     for( i =0 ; i<edgeList.size(); i++ )
     if (edgeList[i]->edge_type == e->edge_type )
     {
           if ((edgeList[i]->edge_type == singular &&
               edgeList[i]->arc_id != e->arc_id) ||
               (edgeList[i]->edge_type == drilled &&
               edgeList[i]->link_id != e->link_id))
               continue;

           edgeList[i]->selected = !edgeList[i]->selected;
           edgeList[i]->thickness = (edgeList[i]->selected) ? THICK : THIN;
     }

//     buffer.fill( CANVAS );
//     update();
}



void DiagramCanvas::drillToggle()
{
    int i;

    drill_on = (drill_on) ? FALSE : TRUE;

    for( i=0; i<edgeList.size(); i++)
    if (edgeList[i]->selected)
        edgeList[i]->edge_type = (drill_on) ? drilled : singular;

    if (edge_creation_in_progress)
        currentEdge->edge_type = (drill_on) ? drilled : singular;

//    buffer.fill( CANVAS );
    assign_arcs();
    assign_links();
//    update();
}


void DiagramCanvas::assign_arcs()
{
    std::queue<Edge *> queue;
    int i;

    bool *visited = new bool[edgeList.size()];
    bool drilled_arc;
    num_arcs = 0;

    for( i = 0; i < edgeList.size(); i++ )
    {
          edgeList[i]->edge_id = i;
          edgeList[i]->arc_id = -1;
          visited[i] = FALSE;

          if (queue.empty() && !visited[i])
          {
                queue.push(edgeList[i]);
                visited[i] = TRUE;
		drilled_arc = (edgeList[i]->edge_type == drilled);

		if (!drilled_arc)
                	edgeList[i]->arc_id = num_arcs;
          }
    }

    while (!queue.empty())
    {
          while (!queue.empty())
          {
                  Edge *e = queue.front();

                  for ( i = 0; i < 2; i++ )
                  {
                          Vertex *v = e->vertex[i];

                          if ( e->vertex[i]->incidentEndData.size() == 2 )
                          {
                                 int j = (v->incidentEndData[0]->edge == e) ? 1 : 0;
                                 Edge *e1 = v->incidentEndData[j]->edge;

                                 bool needsSwitch = (v->incidentEndData[j]->type==i);
                                 if ((edge_creation_in_progress &&
                                     e1 == currentEdge ) ||
					e->edge_type != e1->edge_type )
                                     continue;

                                 if (!visited[e1->edge_id])
                                 {
                                         visited[e1->edge_id] = TRUE;

					 if (!drilled_arc)
                                         	e1->arc_id = num_arcs;

                                         if (needsSwitch)
                                         {
                                                Vertex *temp;
                                                temp = e1->vertex[0];
                                                e1->vertex[0] = e1->vertex[1];
                                                e1->vertex[1] = temp;

                                                for(j=0;j<2;j++)
                                                for(int k=0;k<e1->vertex[j]->incidentEndData.size(); k++ )
                                                if (e1->vertex[j]->incidentEndData[k]->edge == e1 )
                                                        e1->vertex[j]->incidentEndData[k]->type =
                                                                (e1->vertex[j]->incidentEndData[k]->type==begin) ? end : begin;

						std::vector<Crossing *>::iterator c_it;
						c_it = crossingList.begin();
						while ( c_it < crossingList.end() )
						{
							Crossing *c = *c_it;
							if ( c->over == e1 || c->under == e1 )
							{
								refresh_crossing( c );
//								update();
							}
							++c_it;
						}
                                         }
                                         queue.push(e1);
                                 }
                          }
                  }

                  queue.pop();
          }

	if (!drilled_arc)
          	num_arcs++;

          for ( i = 0; i< edgeList.size(); i++ )
          if (!visited[i])
          {
                  queue.push(edgeList[i]);
                  visited[edgeList[i]->edge_id] = TRUE;

		  drilled_arc = (edgeList[i]->edge_type == drilled);

		if (!drilled_arc)
			edgeList[i]->arc_id = num_arcs;

                  break;
          }
    }


}


void DiagramCanvas::assign_links()
{
    std::queue<Edge *> queue;
    int i;
    bool *visited = new bool[edgeList.size()];
    num_links = 0;

    for( i = 0; i < edgeList.size(); i++ )
    {
          edgeList[i]->edge_id = i;
          edgeList[i]->link_id = -2;
          visited[i] = (edgeList[i]->edge_type == drilled) ? FALSE : TRUE;

          if (queue.empty() && !visited[i] )
          {
                queue.push(edgeList[i]);
                visited[i] = TRUE;
                edgeList[i]->link_id = num_links;
          }
    }

    while (!queue.empty())
    {
          while (!queue.empty())
          {
                  Edge *e = queue.front();

                  for ( i = 0; i < 2; i++ )
                  {
                          Vertex *v = e->vertex[i];
                          int j;

                          for ( j = 0; j <  v->incidentEndData.size(); j++ )
                          {
                                 Edge *e1 = v->incidentEndData[j]->edge;

                                 if ( edge_creation_in_progress && e1 == currentEdge )
                                      continue;
                                 else if (!visited[e1->edge_id])
                                 {
                                      visited[e1->edge_id] = TRUE;
                                      e1->link_id = num_links;
                                      queue.push(e1);
                                 }
                          }
                  }

                  queue.pop();
          }

          num_links++;

          for ( i = 0; i< edgeList.size(); i++ )
          if (!visited[i])
          {
                  queue.push(edgeList[i]);
                  visited[edgeList[i]->edge_id] = TRUE;
                  edgeList[i]->link_id = num_links;
                  break;
          }
    }

	/* now for the orbifold links */ 

	if (!edge_creation_in_progress)
    for( i=0; i<vertexList.size();i++)
    {
            int j;
            Vertex *v = vertexList[i];
            v->link_id = -1;

	for( j = 0; j < v->incidentEndData.size() && v->link_id == -1; j++ )
	if ( v->incidentEndData[j]->edge->edge_type == drilled )
		v->link_id = v->incidentEndData[j]->edge->link_id;

	if (v->link_id < 0 && v->incidentEndData.size() > 2 )
		v->link_id = num_links++;

	if (v->link_id < 0 && v->incidentEndData.size() == 2 )
	{
		/* check if v lies on a singular loop */
		bool singular_loop = TRUE;
		Edge *e = v->incidentEndData[0]->edge;
		int arc = e->arc_id;
		int k;

		for(j=0;j<edgeList.size() && singular_loop;j++)
		if ( edgeList[j]->arc_id == arc )
			for(k=0;k<2 && singular_loop;k++)
			{
				if (edgeList[j]->vertex[k]->incidentEndData.size() != 2 ||
					edgeList[j]->vertex[k]->link_id > -1 )
					singular_loop = FALSE;
				else if ( edgeList[j]->vertex[k]->incidentEndData[0]->edge->edge_type == drilled ||
					edgeList[j]->vertex[k]->incidentEndData[1]->edge->edge_type == drilled )
					singular_loop = FALSE;
			}


		if (singular_loop)
			v->link_id = num_links++;
	}
    }

}

void DiagramCanvas::assign_cuffs()
{

	num_cuffs = 0;

	for( int i = 0; i < edgeList.size(); i++)
		edgeList[i]->cuff_id = -1;

	for( int i = 0; i < num_links; i++ )
	{
		bool singular = FALSE;
		bool cuff_marked = FALSE;

		for( int j = 0; j < vertexList.size() && !singular; j ++ )
		if (vertexList[j]->link_id ==  i )
			for( int k = 0;
				 k < vertexList[j]->incidentEndData.size() && !singular;
				 k++ )
			if ( vertexList[j]->incidentEndData[k]->edge->edge_type == singular)
				singular = TRUE;

		if  (singular) continue;

		for( int j = 0; j < vertexList.size(); j ++ )
		if (vertexList[j]->link_id ==  i
			&& vertexList[j]->incidentEndData.size() ==3 )
		for( int k = 0; k < 3; k++ )
		{
		/* mark cuffs */	
			cuff_marked = TRUE;
			Edge *e = vertexList[j]->incidentEndData[k]->edge;
			EndType type = vertexList[j]->incidentEndData[k]->type;

			if (e->cuff_id != -1 )
				continue;
	
			e->cuff_id = num_cuffs;

			while ( TRUE )
			{
				type = (type == begin) ? end : begin;

				Vertex *v = e->vertex[type];

				if ( v->incidentEndData.size() > 2 )
					break;
				int m;

				for( m = 0; m < 2; m ++ )
				if ( v->incidentEndData[m]->edge != e ) 
					break;

				e = v->incidentEndData[m]->edge;
				e->cuff_id = num_cuffs;
				type = v->incidentEndData[m]->type;
			};

			num_cuffs++;
		}

		if (!cuff_marked)
		{
			/* must have been a torus */
			for (int j = 0; j < edgeList.size(); j++ )
			if( edgeList[j]->link_id == i )
				edgeList[j]->cuff_id = num_cuffs;

			num_cuffs++;
		}
	}

}

void DiagramCanvas::prepare_components_for_output()
{
    int     link_id = -5, i, j, k;

    for(i=0;i<num_arcs;i++)
    {
	for(j=0;j<edgeList.size();j++)
	{
		Edge *e = edgeList[j];
		if (e->arc_id==i)
			for(k=0;k<2;k++)
			if (e->vertex[k]->link_id != -1 )
			{
				e->vertex[k]->incidentEndData[get_strand(e,e->vertex[k])]->singular=(k==1);
				if (k==0) link_id = e->vertex[k]->link_id;
			}	
	}
	for( j = 0; j<edgeList.size();j++)
	if (edgeList[j]->arc_id == i )
		edgeList[j]->link_id = link_id;
    }

    for(i=0;i<vertexList.size();i++)
	if (vertexList[i]->link_id == -1 )
	{
		for( int j = 0; j < vertexList[i]->incidentEndData.size(); j++)
		if ( vertexList[i]->incidentEndData[j]->edge->edge_type == drilled )
			vertexList[i]->link_id =
				vertexList[i]->incidentEndData[j]->edge->link_id;

	}
}

bool DiagramCanvas::isReadyForTriangulation()
{
       int i;

       for(i=0;i<vertexList.size(); i++)
       if (vertexList[i]->incidentEndData.size() < 2)
            return FALSE;

       return TRUE;
}
