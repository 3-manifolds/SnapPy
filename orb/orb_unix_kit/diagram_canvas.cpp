#include "diagram_canvas.h"

#include "SnapPea.h"


/*
#include <qapplication.h>
#include <qlayout.h>
#include <qevent.h>
#include <qpainter.h>
#include <qpen.h>
#include <qrect.h>
#include <qregion.h>
#include <qpoint.h>
#include <qwmatrix.h>
*/
#include <complex>
//#include <qpointarray.h>
//#include <qfiledialog.h>
#include <cstdio>
#include <string>

extern QPoint time_to_point( Edge *edge, double t );
extern QPoint complex_to_point( std::complex<double> z );
extern std::complex<double> point_to_complex( QPoint p );

DiagramCanvas::DiagramCanvas( /*QWidget *parent, */ const char *name, bool r )
//: /*QWidget( parent, name ), buffer( width(), height() )*/
{
	readOnly = r;

	initialize_colors();
//	setBackgroundColor( CANVAS );
//	buffer.resize( 600, 600 );
//	buffer.fill( CANVAS );
//	setBackgroundMode( QWidget::NoBackground );
//	setMouseTracking( TRUE );
	edge_creation_in_progress = FALSE;
	vertex_drag_in_progress = FALSE;
	proximity_tolerance = 3;
	num_arcs = 0;
	num_links = 0;
	num_cuffs = 0;
	untouched = TRUE;
	drill_on = FALSE;

}

DiagramCanvas::~DiagramCanvas()
{
	clearDiagram();
}

void DiagramCanvas::setReadOnly()
{
	readOnly = TRUE;
}

void DiagramCanvas::deselectEdges()
{
	for( int i = 0; i < (int) edgeList.size(); i++)
	{
		edgeList[i]->selected = FALSE;
		edgeList[i]->thickness = THIN;
	}

//     buffer.fill( CANVAS );
        //   show();
//     update();
}

#if 0
void DiagramCanvas::mousePressEvent(QMouseEvent *e)
{
    if (readOnly)
    {
    	for( int i = 0; i<edgeList.size();i++)
    	if (proximity(e->pos(),edgeList[i]) && edgeList[i]->edge_type == singular )
                	emit selectionChanged( edgeList[i]->arc_id );
	return;
    }


    untouched = FALSE;


    if (e->button() == LeftButton )
    {
          if (edge_creation_in_progress)
          {
                    int vertex_index;

                    if (completed_edge_general_position(currentEdge, &vertex_index ))
                    {
                              if (joinState == join)
                              {
                                      Vertex  *vertex, *new_vertex;

                                      vertex = vertexList[vertex_index];

                                      paint( currentEdge, TRUE );
                                      currentEdge->vertex[end]->~Vertex();
                                      currentEdge->vertex[end] = vertex;

                                      EndData *ed = new EndData( currentEdge, end );

                                      vertex->incidentEndData.insert(
                                             vertex->incidentEndData.end(), ed );

                                      edgeList.insert( edgeList.end(), currentEdge );

				      new_vertex = NULL;
					Edge *edge = currentEdge;
					currentEdge = NULL;
					edge_creation_in_progress = FALSE;
                                      assign_arcs();
                                      assign_links();
						assign_new_crossings( edge);
						paint( edge, FALSE );
                              }
                              else
                              {
                                      Vertex *vertex, *new_vertex;

                                      vertex = currentEdge->vertex[end];
                                      vertexList.insert( vertexList.end(), vertex );

                                      edgeList.insert( edgeList.end(), currentEdge );

                                      new_vertex = new Vertex( vertex->position );
					Edge *edge = currentEdge;
                                      currentEdge = new Edge( vertex, new_vertex );
                                      currentEdge->edge_type = (drill_on) ? drilled : singular;

                                      edge_creation_in_progress = TRUE;
                                      assign_arcs();
                                      assign_links();
					assign_new_crossings( edge );
                              }
                    }
          }
          else
          {
                   int vertex_index;

                   if (new_point_general_position(e->pos(), &vertex_index ))
                   {
                              if (joinState == join)
                              {
                                      Vertex *vertex, *new_vertex;

                                      vertex = vertexList[vertex_index];
                                      new_vertex = new Vertex( vertex->position );
                                            
                                      currentEdge = new Edge( vertex, new_vertex);
                                      currentEdge->edge_type = (drill_on) ? drilled : singular;
                                     
                                      edge_creation_in_progress = TRUE;
                                      assign_arcs();
                                      assign_links();
                              }
                              else
                              {
                                      Vertex *new_vertex1, *new_vertex2;

                                      new_vertex1 = new Vertex( e->pos() );
                                      vertexList.insert( vertexList.end(), new_vertex1 );
                                      new_vertex2 = new Vertex( e->pos() );

                                      currentEdge = new Edge( new_vertex1, new_vertex2);
                                      currentEdge->edge_type = (drill_on) ? drilled : singular;

                                      edge_creation_in_progress = TRUE;
                                      assign_arcs();
                                      assign_links();
                              }
                   }
                   else
                   {
                             bool close_to_crossing = FALSE;
                             int i;

                             for( i = 0; i<crossingList.size();i++)
                             if (proximity(e->pos(),crossingList[i]->position))
                             {
                                        crossingList[i]->switchCrossing();
                                        update();
                                    	 close_to_crossing = TRUE;
                             }


                             if (!close_to_crossing)
                                     for( i = 0; i<edgeList.size();i++)
                                     if (proximity(e->pos(),edgeList[i]))
                                            invert_arc(edgeList[i]);

                   }
          }
    }

    if (e->button() == RightButton )
    {
                 if (edge_creation_in_progress)
                 {
                             Vertex *vertex;

                             paint( currentEdge, TRUE );

                             vertex = currentEdge->vertex[begin];
                             currentEdge->vertex[end]->~Vertex();

                             currentEdge->~Edge();

                             if (vertex->incidentEndData.size() == 1)
                             {
                                 vector<Vertex *>::iterator v_it;

                                 v_it = vertexList.begin();
                                 while( v_it != vertexList.end() )
                                 if (*v_it == vertex)
                                     break;
                                 else
                                     v_it++;

                                 (*v_it)->~Vertex();
                                 vertexList.erase(v_it);       

                             }
                             else
                             {
                                 vector<EndData *>::iterator ed_it;

                                 ed_it = vertex->incidentEndData.end() - 1;
                                 (*ed_it)->~EndData();
                                 vertex->incidentEndData.erase( ed_it);
                             }

                             edge_creation_in_progress = FALSE;
                             assign_arcs();
                             assign_links();
                 }
                 else
                 {
                             int i;
                             bool action_taken = FALSE;

                             for( i = 0; i<crossingList.size();i++)
                             if (proximity(e->pos(),crossingList[i]->position))
                             {
                                        crossingList[i]->switchCrossing();
                                        update();
                                        action_taken = TRUE;
                             }

                             if (!action_taken) for( i = 0; i<vertexList.size();i++)
                             if (proximity(e->pos(),vertexList[i]->position))
                             {
                                        draggedVertex = vertexList[i];
                                        saved_vertex_position = draggedVertex->position;
                                        delete_dragged_crossings();
                                        vertex_drag_in_progress = TRUE;
                                        action_taken = TRUE;
                             }

                             if (!action_taken) for( i = 0; i<edgeList.size();i++)
                             if (proximity(e->pos(),edgeList[i]))
                             {
                                       if (edgeList[i]->selected)
                                       {
                                                paint( edgeList[i], TRUE );
                                                edgeList[i]->selected = FALSE;
                                                edgeList[i]->thickness = THIN;
                                                paint( edgeList[i], FALSE);
                                       }
                                       else
                                       {
                                                paint( edgeList[i], TRUE );
                                                edgeList[i]->selected = TRUE;
                                                edgeList[i]->thickness = THICK; 
                                                paint( edgeList[i], FALSE);
                                       }
                                       update();
                             }
                 }
    }

    if (e->button() == MidButton )
          emit drillToggled();
}


void DiagramCanvas::mouseReleaseEvent(QMouseEvent *e)
{
    vector< EndData *>::iterator ed_it;

    if ( vertex_drag_in_progress)
    {
           if (!dragged_vertex_general_position())
           {
                     ed_it = draggedVertex->incidentEndData.begin();

                     while( ed_it != draggedVertex->incidentEndData.end() )
                     {
                              paint( (*ed_it)->edge, TRUE );
                              ed_it++;
                     }

                     draggedVertex->position = saved_vertex_position;

                     ed_it = draggedVertex->incidentEndData.begin();

                     while( ed_it != draggedVertex->incidentEndData.end() )
                     {
                              paint( (*ed_it)->edge, FALSE );
                              ed_it++;
                     }
           }

           ed_it = draggedVertex->incidentEndData.begin();

           while( ed_it != draggedVertex->incidentEndData.end() )
           {
                      assign_new_crossings( (*ed_it)->edge );
                      ed_it++;
           }
           vertex_drag_in_progress = FALSE;
           update();
    }
}

void DiagramCanvas::mouseMoveEvent(QMouseEvent *e)
{
    int i;
    bool status_bar_changed = FALSE, near_vertex = FALSE;
    QString status;

    setCursor(Qt::arrowCursor);

    if (edge_creation_in_progress)
         if (!completed_edge_general_position(currentEdge, &i ))
         {
              status.sprintf(" Right click to kill cord."); 
              status_bar_changed = TRUE;
         }

    for(i=0;i<vertexList.size();i++)
    if (proximity(e->pos(),vertexList[i]->position))
    {
          near_vertex = TRUE;
          setCursor(Qt::pointingHandCursor);

          if (!status_bar_changed && !edge_creation_in_progress)
          {
               if (vertexList[i]->link_id == -1)
		{
			if (!readOnly)
		       		status.sprintf("Right click and drag to move.");
		}
	       else    status.sprintf(" Vertex:\t%d", vertexList[i]->link_id+1 ); 

               status_bar_changed = TRUE;
          }
    }

    if (!edge_creation_in_progress && !status_bar_changed && !near_vertex )
         for(i=0;i<edgeList.size();i++)
         if (proximity(e->pos(),edgeList[i]))
         {
              setCursor(Qt::pointingHandCursor);

              if (edgeList[i]->edge_type == singular)
                      status.sprintf(" Coloured Edge:\t%d", edgeList[i]->arc_id+1);
              else if (readOnly && edgeList[i]->cuff_id > -1 )
                      status.sprintf(" Vertex:\t%d\tCuff:\t%d",
				edgeList[i]->link_id+1, edgeList[i]->cuff_id+1);
     	      else  status.sprintf(" Vertex:\t%d", edgeList[i]->link_id+1);

              status_bar_changed = TRUE;
         }

    if (edge_creation_in_progress)
        redrawMovedEdge( e );
    else for(i=0;i<crossingList.size();i++)
    if (!readOnly && proximity(e->pos(),crossingList[i]->position))
    {
          setCursor(Qt::pointingHandCursor);
	  status.sprintf(" Click to switch crossing.");
	  status_bar_changed = TRUE;
    }

    if (vertex_drag_in_progress)
        redrawDraggedVertex( e );

    if (!status_bar_changed)
         status.sprintf(" ");

    emit statusBarChanged(status);
}

void DiagramCanvas::resizeEvent(QResizeEvent *e)
{
    buffer.resize( size() );
    buffer.fill( CANVAS );
    update();
}

void DiagramCanvas::paintEvent(QPaintEvent *e)
{
    QRect r = e->rect();
    paint_all( r);
}


void DiagramCanvas::paint_all( QRect& r )
{
    int i;
    QPainter p( &buffer );
    QPen pen;
    QBrush brush;
    QPointArray a;
    QColor color;
    QRegion mask, crossing_disk;
    set<double>::iterator dbl_it;

    for ( i=0; i<crossingList.size(); ++i )
    {
        Crossing *c = crossingList[i];

        c->under->underpasses.insert( c->under->underpasses.end(), c->position_on_understrand );
    }

    for ( i=0; i<edgeList.size(); ++i )
    {
        Edge *e = edgeList[i];

        pen.setWidth( e->thickness );
        pen.setColor( CANVAS );
        p.setPen( pen );

        QPoint v, v1, v2, p1, p2 ;

        p1 = e->vertex[begin]->position;
        p2 = e->vertex[end]->position;
        p.drawLine( p1, p2 );
        */
/*
      if (e->edge_type==singular)
        {

                v = p1 - p2;
                v1.setX( v.x() - v.y() );
                v1.setY( v.x() + v.y() );
                v2.setX( v.x() + v.y() );
                v2.setY( -v.x() + v.y() );

                v1 = 5 * v1 / sqrt( double (v1.x() * v1.x() + v1.y() * v1.y()) );
                v2 = 5 * v2 / sqrt( double (v2.x() * v2.x() + v2.y() * v2.y()) );
                v  = v / 2;
                pen.setWidth( ARROW );
                p.setPen( pen );
                p.drawLine( p2+v, p2+v1+v );
                p.drawLine( p2+v, p2+v2+v );
        }
*/
/*
    }

    for ( i=0; i<edgeList.size(); ++i )
    {
        Edge *e = edgeList[i];

        color = (e->edge_type == drilled) ? DRILLED : colorList[e->arc_id % 18 ];

        pen.setWidth( e->thickness );
        pen.setColor( color );
        p.setPen( pen );

        brush.setColor( color );
        brush.setStyle( SolidPattern );
        p.setBrush( brush );

        if ( e->underpasses.size() > 0 )
        {
            p.setClipping( TRUE );
            mask = QRegion( 0, 0, width(), height(), QRegion::Rectangle );

            for ( dbl_it=e->underpasses.begin(); dbl_it!=e->underpasses.end(); ++dbl_it )
            {
                QPoint center = time_to_point( e, *dbl_it );
                crossing_disk = QRegion( center.x()-7, center.y()-7, 14, 14 , QRegion::Ellipse );
                mask -= crossing_disk;
            }

            p.setClipRegion( mask );
            QPoint v, v1, v2, p1, p2 ;

            p1 = e->vertex[begin]->position;
            p2 = e->vertex[end]->position;
            p.drawLine( p1, p2 );
*/
/*
            if (e->edge_type==singular)
            {
                v = p1 - p2;
                v1.setX( v.x() - v.y() );
                v1.setY( v.x() + v.y() );
                v2.setX( v.x() + v.y() );
                v2.setY( -v.x() + v.y() );

                v1 = 5 * v1 / sqrt( double (v1.x() * v1.x() + v1.y() * v1.y()) );
                v2 = 5 * v2 / sqrt( double (v2.x() * v2.x() + v2.y() * v2.y()) );
                v  = v / 2;
                pen.setWidth( ARROW );
                p.setPen( pen );
                p.drawLine( p2+v, p2+v1+v );
                p.drawLine( p2+v, p2+v2+v );
            }
*/
            p.setClipping( FALSE );
        }

        else
        {
            QPoint v, v1, v2, p1, p2 ;

            p1 = e->vertex[begin]->position;
            p2 = e->vertex[end]->position;
            p.drawLine( p1, p2 );
/*
            if (e->edge_type==singular)
            {
                v = p1 - p2;
                v1.setX( v.x() - v.y() );
                v1.setY( v.x() + v.y() );
                v2.setX( v.x() + v.y() );
                v2.setY( -v.x() + v.y() );

                v1 = 5 * v1 / sqrt( double (v1.x() * v1.x() + v1.y() * v1.y()) );
                v2 = 5 * v2 / sqrt( double (v2.x() * v2.x() + v2.y() * v2.y()) );
                v  = v / 2;
                pen.setWidth( ARROW );
                p.setPen( pen );
                p.drawLine( p2+v, p2+v1+v );
                p.drawLine( p2+v, p2+v2+v );
            }
*/
        }

        e->underpasses.clear();
    }

    p.end();

    bitBlt( this, r.x(), r.y(), &buffer, r.x(), r.y(), r.width(), r.height() );

}




void DiagramCanvas::paint( Edge *edge, bool paint_over )
{
    QPainter p( &buffer );
    QPen pen;
    QBrush brush;
    QColor paint_edge_color, edge_color;

    if (edge_creation_in_progress && edge == currentEdge)
        edge_color = NEW_EDGE; 
    else edge_color = (edge->edge_type == drilled) ? DRILLED : colorList[edge->arc_id % 18 ];

    paint_edge_color = (paint_over) ? CANVAS : edge_color;

    pen.setWidth( edge->thickness );
    pen.setColor( paint_edge_color );

    p.setPen( pen );

    QPoint v, v1, v2, p1, p2 ;

    p1 = edge->vertex[begin]->position;
    p2 = edge->vertex[end]->position;
    p.drawLine( p1, p2 );
/*
    if (edge->edge_type==singular && edge_color != NEW_EDGE)
    {
        v = p1 - p2;
        v1.setX( v.x() - v.y() );
        v1.setY( v.x() + v.y() );
        v2.setX( v.x() + v.y() );
        v2.setY( -v.x() + v.y() );

        v1 = 5 * v1 / sqrt( double (v1.x() * v1.x() + v1.y() * v1.y()) );
        v2 = 5 * v2 / sqrt( double (v2.x() * v2.x() + v2.y() * v2.y()) );
        v  = v / 2;
        pen.setWidth( ARROW );
        p.setPen( pen );
        p.drawLine( p2+v, p2+v1+v );
        p.drawLine( p2+v, p2+v2+v );
    }
*/
    p.end();

    update();
}

void DiagramCanvas::redrawMovedEdge( QMouseEvent *e )
{
    QPoint newEndPoint = e->pos();

    paint( currentEdge, TRUE);

    currentEdge->vertex[end]->position = newEndPoint;

    paint( currentEdge, FALSE );
}

void DiagramCanvas::redrawDraggedVertex( QMouseEvent *e )
{
    vector<EndData *>::iterator ed_it;

    QPoint newPosition = e->pos();
    ed_it = draggedVertex->incidentEndData.begin();

    while( ed_it != draggedVertex->incidentEndData.end() )
    {
            paint( (*ed_it)->edge, TRUE );
            ed_it++;
    }

    draggedVertex->position = newPosition;

    ed_it = draggedVertex->incidentEndData.begin();

    while( ed_it != draggedVertex->incidentEndData.end() )
    {
            paint( (*ed_it)->edge, FALSE );
            ed_it++;
    }
    
}

#endif

void DiagramCanvas::delete_dragged_crossings()
{
    std::vector<Crossing *>::iterator i;
   
    i = crossingList.begin();
   
    while ( i < crossingList.end() )
    {
        Crossing *c = *i;

        if ( c->over->vertex[begin] == draggedVertex || c->over->vertex[end] == draggedVertex
                || c->under->vertex[begin] == draggedVertex || c->under->vertex[end] == draggedVertex )
        {
            crossingList.erase( i );
//            update();
        }
        else ++i;
    }
}

bool DiagramCanvas::isUntouched()
{
	return untouched;
}

Triangulation * diagram_data_to_triangulation(const char *d)
{
    std::stringstream s;

    s << d;
    
    DiagramCanvas c;
    c.readDiagram(s);

    Triangulation * t = c.outputTriangulation();

    basic_simplification(t);

    return t;
}
