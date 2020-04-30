#include "diagram_canvas.h"


//#include <qpoint.h>
#include <complex>
#include <cstdio>


extern std::complex<double> point_to_complex( QPoint p );

bool DiagramCanvas::new_point_general_position( QPoint point, int *vertex_index )
{
    int i;

    if ( edgeList.size() == 0 )
    {
        joinState = no_join;
        return TRUE;
    }

    joinState = no_join;

    /*
     *  First check whether point is near to some vertex
     */

    for ( i=0; i<vertexList.size(); ++i ) if ( proximity( point, vertexList[i]->position ) )
    {
            joinState = join;
            *vertex_index = i;
            return TRUE;
    }


    /*
     *  Now check whether point is near some illegal part of the diagram
     */

    for ( i=0; i<edgeList.size(); ++i ) if ( proximity( point, edgeList[i] ) )
        return FALSE;

    return TRUE;
}

bool DiagramCanvas::completed_edge_general_position( Edge *edge, int *vertex_index )
{
    int i, j;

    if ( edgeList.size() == 0 ) return TRUE;

    if (proximity(edge->vertex[begin]->position, edge->vertex[end]->position))
          return FALSE;

    joinState = no_join;

    /*
     *  First check whether edge->end is near to some vertex
     */

    for ( i=0; i<vertexList.size(); ++i ) if ( proximity( edge->vertex[end]->position, vertexList[i]->position ) )
    { 
         joinState = join;
         *vertex_index = i;
    }

    if ( joinState == no_join ) *vertex_index = -1;

    /*
     * make sure there is no edge already running between these two vertices.
     */

    if (joinState == join)
    for(i=0;i<vertexList[*vertex_index]->incidentEndData.size();i++)
    {
          int type = (vertexList[*vertex_index]->incidentEndData[i]->type == begin) ? end : begin;
          
          if (vertexList[*vertex_index]->incidentEndData[i]->edge->vertex[type] == edge->vertex[begin] )
                 return FALSE;
    }

    /*
     *  Check whether edge->end is near some illegal edge
     */


    for ( i=0; i<edgeList.size(); ++i )
    if ( proximity( edge->vertex[end]->position, edgeList[i] ) )
    {
        if (joinState == no_join &&
             (edge->vertex[begin] != edgeList[i]->vertex[begin]
             || edge->vertex[begin] != edgeList[i]->vertex[end] ))
                return FALSE;
        else if ( edgeList[i]->vertex[end] != vertexList[*vertex_index]
                   && edgeList[i]->vertex[begin] != vertexList[*vertex_index])
                return FALSE;
    }


    /*
     *  Check whether some existing vertex is near edge
     */
    for ( i=0; i<vertexList.size(); ++i )
    if ( i != *vertex_index
           && edge->vertex[begin] != vertexList[i]
           && proximity( vertexList[i]->position, edge ))
        return FALSE;
    /*
     *  Check whether some existing crossing is near edge
     */

    for ( i=0; i<crossingList.size(); ++i ) if ( proximity( crossingList[i]->position, edge ) ) return FALSE;

    /*
     *  Check whether any new crossing is near any existing edge or crossing
     */

    std::vector<double> cp;
    double t1, t2;
    QPoint crossing_point;
    double len = abs(point_to_complex( edge->vertex[end]->position - edge->vertex[begin]->position));

    for ( i=0; i<edgeList.size(); ++i ) if ( edge_intersect( edge, edgeList[i], &t1, &t2, &crossing_point ) )
    {
        if ( t1 * len < 2*proximity_tolerance || ( 1 - t1 ) * len < 2*proximity_tolerance )  // too close to endpoints
        {
            cp.clear();
            return FALSE;
        }

        for ( j=0; j<cp.size(); ++j ) if ( abs( cp[j] - t1 ) * len < 2*proximity_tolerance )  // crossings too close together
        {
            cp.clear();
            return FALSE;
        }

        cp.insert( cp.end(), t1 );

        for ( j=0; j<edgeList.size(); ++j ) if ( j != i && proximity( crossing_point, edgeList[j] ) )
        {
            cp.clear();
            return FALSE;
        }

        for ( j=0; j<crossingList.size(); ++j ) if ( proximity( crossing_point, crossingList[j]->position ) )
        {
            cp.clear();
            return FALSE;
        }
    }

    cp.clear();

    return TRUE;
}


bool DiagramCanvas::dragged_vertex_general_position()
{
    int i,j;
    std::vector<EndData *>::iterator ed_it;
    std::vector<EndData *>::iterator ed_it1;

    for(i=0;i<edgeList.size();i++)
    if (proximity(draggedVertex->position,edgeList[i])
        && edgeList[i]->vertex[begin] != draggedVertex
        && edgeList[i]->vertex[end] != draggedVertex )
        return FALSE;

    for(ed_it = draggedVertex->incidentEndData.begin();
        ed_it != draggedVertex->incidentEndData.end();
        ed_it++)
    {
        for(i=0;i<vertexList.size();i++)
        if (proximity(vertexList[i]->position,(*ed_it)->edge)
            && (*ed_it)->edge->vertex[begin] != draggedVertex
            && (*ed_it)->edge->vertex[end] != draggedVertex )
            return FALSE;

        for(i=0;i<crossingList.size();i++)
        if (proximity(crossingList[i]->position,(*ed_it)->edge))
            return FALSE;
    }

    std::vector<double> cp;
    double t1, t2;
    QPoint crossing_point;
    double len;

    for ( ed_it = draggedVertex->incidentEndData.begin();
          ed_it != draggedVertex->incidentEndData.end();
          ed_it++ )
    {
        len = abs( point_to_complex( (*ed_it)->edge->vertex[end]->position
                    - (*ed_it)->edge->vertex[begin]->position));

        for ( i=0; i<edgeList.size(); ++i ) if ( edgeList[i] != (*ed_it)->edge &&
            edge_intersect( (*ed_it)->edge, edgeList[i], &t1, &t2, &crossing_point ) )
        {
            if ( t1 * len < 2*proximity_tolerance
                  || ( 1 - t1 ) * len < 2*proximity_tolerance )  // too close to endpoints
            {
                cp.clear();
                return FALSE;
            }

            for ( j=0; j<cp.size(); ++j )
            if ( abs( cp[j] - t1 ) * len < 2*proximity_tolerance )  // crossings too close together
            {
                cp.clear();
                return FALSE;
            }

            cp.insert( cp.end(), t1 );

            for ( j=0; j<edgeList.size(); ++j )
            if ( j != i &&  proximity( crossing_point, edgeList[j] ) )
            {
                bool not_dragged_edge = TRUE;

                for( ed_it1 = draggedVertex->incidentEndData.begin();
                     ed_it1 != draggedVertex->incidentEndData.end();
                     ed_it1++ )
                if ( (*ed_it)->edge == edgeList[j] ) 
                     not_dragged_edge = FALSE;

                if (not_dragged_edge)
                {
                     cp.clear();
                     return FALSE;
                }
            }

            for ( j=0; j<crossingList.size(); ++j )
            if ( proximity( crossing_point, crossingList[j]->position ) )
            {
                cp.clear();
                return FALSE;
            }

        }

        cp.clear();
     }
    return TRUE;
}



bool DiagramCanvas::proximity( QPoint point, Edge *edge )
{
    double p = abs( point_to_complex( point - edge->vertex[begin]->position ) );
    double q = abs( point_to_complex( point - edge->vertex[end]->position ) );

    double r = edge->distance_from_edge( point );
    double t = edge->projection_to_edge( point );

    bool close_to_interior = r < proximity_tolerance && t >= 0 && t <= 1;
    bool close_to_begin = p < proximity_tolerance;
    bool close_to_end = q < proximity_tolerance;

    if ( close_to_begin || close_to_end || close_to_interior ) return TRUE;
    else return FALSE;
}



bool  DiagramCanvas::proximity( QPoint point1, QPoint point2 )
{
    double d = abs( point_to_complex( point1 - point2 ) );

    if ( d < 2*proximity_tolerance ) return TRUE;
    else return FALSE;
}


bool DiagramCanvas::edge_intersect( Edge *edge1, Edge *edge2, double *t1, double *t2, QPoint *point )
{
    std::complex<double> z1, z2, z3, z4;

    z1 = point_to_complex( edge1->vertex[begin]->position );
    z2 = point_to_complex( edge1->vertex[end]->position );
    z3 = point_to_complex( edge2->vertex[begin]->position );
    z4 = point_to_complex( edge2->vertex[end]->position );

    if ( z1 == z2 || z3 == z4 ) return FALSE;

    std::complex<double> w3 = (z3 - z1) / (z2 - z1);
    std::complex<double> w4 = (z4 - z1) / (z2 - z1);

    if ( w3.imag() >= 0 && w4.imag() >= 0 || w3.imag() <= 0 && w4.imag() <= 0 ) return FALSE;

    *t1 = w3.real() - w3.imag() * ( w4.real() - w3.real() ) / ( w4.imag() - w3.imag() );

    if ( *t1 <= 0 || *t1 >= 1 ) return FALSE;

    std::complex<double> w1 = (z1 - z3) / (z4 - z3);
    std::complex<double> w2 = (z2 - z3) / (z4 - z3);

    *t2 = w1.real() - w1.imag() * ( w2.real() - w1.real() ) / ( w2.imag() - w1.imag() );

    if ( *t2 <= 0 || *t2 >= 1 ) return FALSE;

    std::complex<double> y = ( 1 - *t1 ) * z1 + *t1 * z2;

    *point = QPoint( (int) y.real() , (int) y.imag() );

    return TRUE;
}


void DiagramCanvas::assign_new_crossings( Edge *edge )
{
    int i;
    double t1, t2;
    QPoint crossing_point;

    for ( i=0; i<edgeList.size(); ++i )
        if ( edge != edgeList[i] && edge_intersect( edge, edgeList[i], &t1, &t2, &crossing_point ) )
    {
        Edge *ue = edgeList[i];  // "under_edge"

        currentCrossing = new Crossing;
        crossingList.insert( crossingList.end(), currentCrossing );
        currentCrossing->position = crossing_point;
        currentCrossing->over = edge;
        currentCrossing->under = ue;
        currentCrossing->position_on_overstrand = t1;
        currentCrossing->position_on_understrand = t2;
    }
}

void DiagramCanvas::refresh_crossing( Crossing *c )
{
     edge_intersect( c->over, c->under, &c->position_on_overstrand, &c->position_on_understrand, &c->position );
}
