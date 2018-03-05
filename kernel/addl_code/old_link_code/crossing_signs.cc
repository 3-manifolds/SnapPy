#include "link_diagram.h"
#include <cstdio>
#include "stdlib.h"

using std::vector;

bool orient(int nc, vector<int>& a, vector<int>& f);

Pass *new_Pass(void);

Pass *new_Pass(void){
  Pass *thePass;
  thePass = new Pass;
  thePass->index = 0;
  thePass->new_index = 0;
  thePass->crossing_sign = 0;
  thePass->orientation = 0;
  thePass->component = 0;
  thePass->alt_sign = 0;
  thePass->crossing = NULL;
  thePass->next = NULL;
  thePass->prev = NULL;
  thePass->partner = NULL;
  thePass->level = over;
  return thePass;
}

bool Link::compute_crossing_signs()
{
    int i, j, current_index, current_new_index, ccomp, *add_cr;
    
    if ( num_components == 1 )
    {
    
        orient( crossing_number, dowker, crossing_sign );
        return true;
        
    }

    Pass **orig_pass = new Pass*[ 2*crossing_number + 1 ];
    Pass **new_pass_1 = new Pass*[ num_components - 1 ];
    Pass **new_pass_2 = new Pass*[ num_components - 1 ];
    
    Pass **join_1 = new Pass*[ num_components - 1 ];
    Pass **join_2 = new Pass*[ num_components - 1 ];
    Pass **join_3 = new Pass*[ num_components - 1 ];
    Pass **join_4 = new Pass*[ num_components - 1 ];
    
    vector<int> a;
    vector<int> f;
    
    a.resize( 2*crossing_number + 2*num_components - 1 );
    f.resize( 2*crossing_number + 2*num_components - 1 );
    
    for ( i=1; i<=2*crossing_number; ++i ) orig_pass[i] = new_Pass();
    
    for ( i=0; i<num_components-1; ++i ) new_pass_1[i] = new_Pass();
    for ( i=0; i<num_components-1; ++i ) new_pass_2[i] = new_Pass();
        
    for ( i=1; i<=2*crossing_number; ++i )
    {
        orig_pass[i]->index = i;
        orig_pass[i]->partner = orig_pass[dowker[i]];
    }
 
    for ( i=0; i<num_components; ++i )
    {
        for ( j=comp_start[i]; j<comp_end[i]; ++j ) orig_pass[j]->next = orig_pass[j+1];   
        orig_pass[ comp_end[i] ]->next = orig_pass[ comp_start[i] ];
    }     
 
    /*
     *  Add extra crossings at add_cr locations
     */
     
    add_cr = new int[ num_components - 1 ];
    
    for ( i=0; i<num_components-1; ++i ) add_cr[i] = dowker[ comp_start[i+1]+1 ];
    
    for ( ccomp=0; ccomp < num_components-1; ++ccomp )
    {
        join_1[ ccomp ] = orig_pass[ add_cr[ ccomp ] ];
        join_2[ ccomp ] = orig_pass[ add_cr[ ccomp ] ]->partner;
        join_3[ ccomp ] = orig_pass[ add_cr[ ccomp ] ]->partner->next;
        join_4[ ccomp ] = orig_pass[ add_cr[ ccomp ] ]->next;
        
        new_pass_1[ ccomp ]->partner = new_pass_2[ ccomp ];
        new_pass_2[ ccomp ]->partner = new_pass_1[ ccomp ];
        
        join_1[ ccomp ]->next = new_pass_1[ ccomp ];
        new_pass_1[ ccomp ]->next = join_3[ ccomp ];
        join_2[ ccomp ]->next = new_pass_2[ ccomp ];
        new_pass_2[ ccomp ]->next = join_4[ ccomp ];
    }
 
    /*
     *  Form dt sequence for resulting knot
     */
 
    Pass *c = orig_pass[1];
    current_new_index = 0;
  
    do {
        c->new_index = ++current_new_index;
        c = c->next;
    } while ( c != orig_pass[1] );
 
    c = orig_pass[1];

    for ( i=1; i<=2*( crossing_number + num_components - 1 ); ++i )
    {
        a[ c->new_index ] = c->partner->new_index;
        c = c->next;
    }
    
    if ( !orient( crossing_number+num_components-1, a, f ) )
    {
      throw 123; 
    }
    
    for (i=1; i<=2*crossing_number; ++i)
        orig_pass[i]->crossing_sign = f[ orig_pass[i]->new_index ];
    
    for (i=1; i<=2*crossing_number; ++i) crossing_sign[i] = orig_pass[i]->crossing_sign;  
    
    /*
     *  Clean up
     */

    for ( i=0; i<num_components-1; ++i ) delete new_pass_1[i];
    for ( i=0; i<num_components-1; ++i ) delete new_pass_2[i];
    delete [] new_pass_1;
    delete [] new_pass_2;
    
    for (i=1; i<=2*crossing_number; ++i) delete orig_pass[i];
    delete [] orig_pass;
    
    delete [] join_1;
    delete [] join_2;
    delete [] join_3;
    delete [] join_4;
    
    delete [] add_cr;
    
    return true;
}


bool orient(int nc, vector<int>& a, vector<int>& f)
{
    int e[100], g[100], h[100], np, n;
    int i, s, t;
 
    np = 2*nc;
    n = 2*nc - 1;

    for (i=1; i<=np; ++i)
    {
        e[i]=0;f[i]=0;g[i]=0;h[i]=0;
    }
    
    f[1]=1; f[a[1]]=-1;
    h[1]=1; h[a[1]]=1;
    e[1]=1;
    t=1;                                          // t is beginning of loop
    while (t!=0)                                  // if t is zero, we are done
    {
    
        e[t]=1;                                   // e[] is inside/outside function
        for (i=t+1; i<=np; ++i)                   // go forward from t, determining e[]
            if ((a[i]>=t) && (a[i]<=a[t]))
                e[i]=-e[i-1];
            else
                e[i]=e[i-1];
        for (i=t-1; i>=1; --i)                    // go backwards from t, determining e[]
            if ((a[i+1]>=t) && (a[i+1]<=a[t]))
                e[i]=-e[i+1];
            else
                e[i]=e[i+1];
        for (i=1; i<=t-1; ++i) g[i]=1;            // g[] specifies points outside loop
        for (i=a[t]+1; i<=np; ++i) g[i]=1;        // which have yet to be dealt with
        s=0;i=1;
        while ((s==0) && (i<=np))
        {
            if (g[i]==1) s=i;                     // s is first point with non-zero g[]
            ++i;
        }
        while (s!=0)
        {                                         // an s exists
            g[s] = g[a[s]] = 0;                   // this s is about to be dealt with
            if ((a[s]<t) || (a[s]>a[t]))          // neither s nor its partner is in loop
            {
                if (e[s]*e[a[s]]==-1)             // sequence is not realizable
                {
                    return false;
                }
            }
            else                                  // partner of s is on loop, i.e. [t, a[t]]
            {                                     // is interlaced with [s, a[s]]
                if (f[s]!=0)                      // sign of s has already been computed
                {
                    if (e[s]*e[a[s]]*f[t]!=f[s])
                    {                             // slightly subtle consistency check shows
                        return false;             // that the sequence is not realizable
                    }
                }
                else                              // sign of s has not been computed, but
                {                                 // we can now compute the sign
                    f[s]=e[s]*e[a[s]]*f[t];       // evaluation of sign of s and its partner
                    f[a[s]]=-f[s];
                    if (((s==1) && (abs(a[np]-a[1])==1)) ||
                      ((s!=1) && ((abs(a[s-1]-a[s])==1) || (abs(a[s-1]-a[s])==n))))
                    {                             // if s, s-1 are on a "twist", we'll get no
                    }                             // new information from the loop based at s
                    else
                    {                             // s, s-1 aren't on a "twist", so eventually
                        h[s]=1;                   // we'll need to process the loop
                        h[a[s]]=1;                // based at s; h[] specifies potential t's
                    }
                }
            }
            s=0;
            i=1;
            while ((s==0) && (i<=np))
            {                                     // carry on getting information from
                if (g[i]==1) s=i;                 // current loop based at t by
                ++i;                              // finding next s
            }
        }
        h[t]=0;                                   // t's usefulness has been exhausted
        h[a[t]]=0;
        t=0;
        i=1;
        while ((t==0) && (i<=np))
        {
            if (h[i]==1) t=i;                     // find next t to deal with
            ++i;
        }
    }
    for (i=1; i<=n+1; ++i) if (f[i]==0) return false;
  
    return true;
}

