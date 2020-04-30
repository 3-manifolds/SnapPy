#include "diagram_canvas.h"

#include <iostream>
//#include <qpoint.h>
#include <complex>
#include <cstdio>

std::complex<double> point_to_complex( QPoint p );
QPoint time_to_point( Edge *edge, double t );
bool ed_less( EndData *ed1, EndData *ed2 );
bool crossing_less( Crossing *c1, Crossing *c2 );
Crossing *get_prev_crossing( Edge *e, Crossing *c);
Crossing *get_prev_crossing( Edge *e, Crossing *c);
void remove_meeting( Graph *graph, int i);



void DiagramCanvas::initialize_colors()
{/*
colorList[0] = Qt::red;
colorList[1] = Qt::green;
colorList[2] = Qt::blue;
colorList[3] = Qt::cyan;
colorList[4] = Qt::magenta;
colorList[5] = QColor(255,180,15);

colorList[6] = QColor( 220, 100, 150 );
colorList[7] = QColor( 150, 220, 100 );
colorList[8] = QColor( 150, 100, 220  );
colorList[9] = QColor( 100, 150, 220 );
colorList[10] = QColor( 220, 150, 100  );
colorList[11] = QColor( 100, 220, 150 );

colorList[12] = QColor( 230, 60, 180 );
colorList[13] = QColor( 180, 230, 60 );
colorList[14] = QColor( 180, 60, 230 );
colorList[15] = QColor( 60, 180, 230 );
colorList[16] = QColor( 230, 180, 60 );
colorList[17] = QColor( 60, 230, 180 );

colorList[18] = Qt::lightGray;
 */

}

void Crossing::switchCrossing()
{
    Edge *temp_edge;
    double t1;

    t1 = position_on_overstrand;
    position_on_overstrand = position_on_understrand;
    position_on_understrand = t1;

    temp_edge = over;
    over = under;
    under = temp_edge;

}

std::complex<double> point_to_complex( QPoint p )
{
    return std::complex<double>( (double) p.x(), (double) p.y() );
}

QPoint time_to_point( Edge *edge, double t )
{
    return (1 - t)*edge->vertex[begin]->position + t*edge->vertex[end]->position;
}

bool ed_more( EndData *ed1, EndData *ed2 )
{
      return ed1->angle > ed2->angle;
}

bool crossing_less( Crossing *c1, Crossing *c2 )
{
      Edge *e;

      e = c1->over;

      if (c2->over == e)
           return c1->position_on_overstrand < c2->position_on_overstrand;
      else if (c2->under == e)
           return c1->position_on_overstrand < c2->position_on_understrand;

      e = c1->under;

      if (c2->over == e)
           return c1->position_on_understrand < c2->position_on_overstrand;
      else if (c2->under == e)
           return c1->position_on_understrand < c2->position_on_understrand;

      printf("ERR ERR\n");

      return TRUE;
}

int get_strand( Edge *e, Vertex *v )
{
     int i;

     for(i=0;i<v->incidentEndData.size();i++)
     if (e==v->incidentEndData[i]->edge)
         return i;

     fprintf(stderr,"ERR: in get_strand.\n");

     return -1;
}

Crossing *get_next_crossing( Edge *e, Crossing *c)
{
      int i;

      for(i=0;i<e->crossings.size();i++)
      if ( e->crossings[i] == c)
          break; 

      return (i+1 < e->crossings.size()) ?  e->crossings[i+1] : NULL;
}


Crossing *get_prev_crossing( Edge *e, Crossing *c)
{
      int i;

      for(i=0;i<e->crossings.size();i++)
      if ( e->crossings[i] == c)
          break;

      return (i==0 || i ==  e->crossings.size())  ? NULL : e->crossings[i-1];
}

void remove_meeting( Graph *graph, int index )
{

 if ( index>= graph->num_meetings)
             return;

 GraphMeeting   *old_array,
                *new_array;
 int            i,
                j;


 old_array = graph->meeting;
 new_array = NEW_ARRAY(graph->num_meetings-1, GraphMeeting );

 for ( i = 0, j=0; i< graph->num_meetings-1; i++ )
 {
        if (i == index) j++;

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

 for( i=0; i< graph->num_meetings;i++)
      for( j=0; j<graph->meeting[i].num_strands;j++)
      if (graph->meeting[i].neighbor[j] >= index)
           graph->meeting[i].neighbor[j]--;

}


