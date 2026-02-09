/*
 *  Constructor for class Link
 *
 */

#include "link_diagram.h"
#include <cstdlib>
#include <cstdio>

using namespace std; 

Crossing :: ~Crossing()
{
}



Link::Link( const string& s )
{
    int i, j;
    
    link_record = string( s );
    
    crossing_number = s[0] - 96;
    
    num_components = s[1] - 96;

    component_length.resize( num_components );
    for ( i=0; i<num_components; ++i ) component_length[i] = 2*( s[i+2] - 96 );
    
    dt.resize( crossing_number + 1 );
    alt_sign.resize( 2*crossing_number+1 );
    
    dt[0] = crossing_number;
    for ( i=1; i<=crossing_number; ++i )
    {
        char c = s[num_components+1+i];
        if ( c > 96 )
        {
            dt[i] = 2*( c - 96 );
            alt_sign[2*i-1] = 1;
        }
        else
        {
            dt[i] = 2*( c - 64 );
            alt_sign[2*i-1] = -1;
        }
    }
    
    /*
     *  Expand dt[] to array dowker[]
     */
    
    dowker.resize( 2*crossing_number+1 );
    for (i=1; i<=crossing_number; ++i) dowker[2*i-1] = dt[i];
    
    for (i=1; i<=2*crossing_number-1; i+=2)
    {
        dowker[dowker[i]] = i;
        alt_sign[dowker[i]] = -alt_sign[i];
    }

    component.resize( 2*crossing_number+1 );
    next_index.resize( 2*crossing_number+1 );
    prev_index.resize( 2*crossing_number+1 );
    comp_start.resize( num_components );
    comp_end.resize( num_components);
    
    comp_start[0] = 1;
    comp_end[0] = component_length[0];
    
    for ( i=1; i<num_components; ++i )
    {
        comp_start[i] = comp_start[i-1] + component_length[i-1];
        comp_end[i] = comp_end[i-1] + component_length[i];
    }
    
    for (i=0; i<num_components; ++i) 
      for ( j=comp_start[i]; j<=comp_end[i]; ++j )
	{
	  component[j] = i;
	  next_index[j] = ( j < comp_end[i] ) ? j+1 : comp_start[i];
	  prev_index[j] = ( j > comp_start[i] ) ? j-1 : comp_end[i];
	}
    

    passList.resize(2*crossing_number+1); 
    for (i=0; i<=2*crossing_number; ++i)
    {
        passList[i] = new Pass;
    }
    
    crossing_sign.resize( 2*crossing_number + 1 );
    compute_crossing_signs();

}

bool Link::make_geometry()
{
    int i, j;
    
    /*
     *  Allocate memory for crossings
     */
    
    crossingList.resize( crossing_number+1 );
    
    for (i=0; i<=crossing_number; ++i)
    {
        crossingList[i] = new Crossing;
    }

    for (i=0; i<num_components; ++i) 
      for ( j=comp_start[i]; j<=comp_end[i]; ++j )
    {
        Pass *p = passList[j];
        
        p->index = j;
        p->component = i;
        p->partner = passList[ dowker[j] ];
        
        p->next = passList[ next_index[j] ];
        p->prev = passList[ prev_index[j] ];
    }
        
    for (i=1; i<=crossing_number; ++i)
    {   
        Crossing *c = crossingList[i];
                 
        c->index = i;
        c->sign = crossing_sign[2*i-1];
        c->alt_sign = alt_sign[2*i-1];
        c->orientation = c->sign * c->alt_sign;
        c->odd = passList[2*i-1];
        c->even = passList[dowker[2*i-1]];
        c->over = ( c->alt_sign == 1 ) ? c->odd : c->even;
        c->under = ( c->alt_sign == 1 ) ? c->even : c->odd; 
        c->over->crossing = c;
        c->under->crossing = c;
        c->over->alt_sign = c->alt_sign;
        c->under->alt_sign = c->alt_sign;
        c->over->level = over;
        c->under->level = under;
    }
}

void Link::make_millett()
{
    int compass[8] = { 4, 2, 3, 1, 2, 4, 3, 1 };
    int i, p, q, r;
    
    make_geometry();
     
    for ( i=0; i<crossing_number; ++i )
    {
        Crossing *c = crossingList[i+1];
        millett_component[i][0] = c->over->component;
        millett_component[i][1] = c->under->component;
 
    }
    
    for ( i=1; i<=crossing_number; ++i )
    {
        Crossing *c = crossingList[i];
        Pass *ps;
        
        millett[i-1][0] = i;
        
        millett[i-1][1] = ( c->orientation == 1 ) ? 1 : 0;
        
        /*
         * p, q, r are the following three attributes of the neighbor strand:
         * p = crossing sign,  q = under/over,  r = backwards/forwards.
         * These three determine the compass direction of base relative to neighbor.
         */
        
        ps = c->over->next;
        millett[i-1][2] = ps->crossing->index;
        p = ( ps->crossing->orientation == 1 ) ? 1 : 0;
        q = ( ps->level == over ) ? 1 : 0;
        r = 0;
        millett[i-1][3] = compass[ 4*p + 2*q + r ];
        
        ps = ( c->orientation == 1 ) ? c->under->prev : c->under->next;
        millett[i-1][4] = ps->crossing->index;
        p = ( ps->crossing->orientation == 1 ) ? 1 : 0;
        q = ( ps->level == over ) ? 1 : 0;
        r = ( c->orientation == 1 ) ? 1 : 0;
        millett[i-1][5] = compass[ 4*p + 2*q + r ];
        
        ps = c->over->prev;
        millett[i-1][6] = ps->crossing->index;
        p = ( ps->crossing->orientation == 1 ) ? 1 : 0;
        q = ( ps->level == over ) ? 1 : 0;
        r = 1;
        millett[i-1][7] = compass[ 4*p + 2*q + r ];
        
        ps = ( c->orientation == 1 ) ? c->under->next : c->under->prev;
        millett[i-1][8] = ps->crossing->index;
        p = ( ps->crossing->orientation == 1 ) ? 1 : 0;
        q = ( ps->level == over ) ? 1 : 0;
        r = ( c->orientation == 1 ) ? 0 : 1;
        millett[i-1][9] = compass[ 4*p + 2*q + r ];
    }

    clean_up_diagram();
}

  
void Link::clean_up_diagram()  //  Clears memory allocated in make_geometry().
{
    int i;

    for (i=0; i<crossingList.size(); ++i) delete crossingList[i];
    
    for (i=0; i<=2*crossing_number; ++i)
    {
        delete passList[i];
    }
   
    passList.clear();
    crossingList.clear();
    next_index.clear();
    prev_index.clear();
}


Link::~Link()
{
}
