#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "graph_complement.h"

Triangulation *graph2manifold(FILE *fp, double order, int *n );

Triangulation *triangulate_graph_complement( Graph * );

static int	read_file( FILE *, int *** );
static Graph	*plantri2locus( int , int **, double ); 
static void	get_meeting_info( GraphMeeting *, int , int , int **);
static void	create_labels_and_components( Graph * );
static void	label_torus_components( Graph *, int );
static void     resize_meeting_array( Graph *, int );
static void	create_components( Graph * );

extern Triangulation *graph2manifold(FILE *fp, double order, int *n )
{
	Triangulation	*manifold;
	EdgeClass	*edge;
	int **code=NULL,nv;
	Graph *gamma = NULL;

	if ((nv = read_file(fp , &code)) == 0 )
		return NULL;

	*n = nv;

	gamma = plantri2locus( nv, code, order );

	if ((manifold = triangulate_graph_complement( gamma ))==NULL)
		return NULL;

	if (order > -1 )
	for(	edge = manifold->edge_list_begin.next;
		edge!=&manifold->edge_list_end;
		edge = edge->next )
		if (edge->is_singular)
		{
			edge->singular_order = order;
			edge->old_singular_order = order;
		}

	return manifold;
}

extern int read_file( FILE *fp, int ***code )
{
	int ne,v,e,nv;
	char	c, first_line[200];

	fgets( first_line, 199, fp );

	if (fscanf( fp, "%d", &nv )!=1 || nv == 0 )
	{
		fclose(fp);
		return 0;
	}

	(*code) = NEW_ARRAY( nv, int * );

	for( v=0; v<nv; v++)	
	{
		fscanf( fp, "%d", &ne );
		(*code)[v] = NEW_ARRAY( ne + 2, int );
		(*code)[v][0] = ne;

		do
			c=getc(fp);
		while (isspace(c));

		(*code)[v][1] = (c=='+') ? 1: ((c=='-') ? -1 : 0);

		for( e=2; e<ne+2; e++)	
			fscanf(fp,"%d",&(*code)[v][e] );
	}

	fclose(fp);

	return nv;
}

static Graph *plantri2locus( int nv, int **code, double order )
{
	Graph *gamma = NULL;
	GraphMeeting *m;
	int i,j;

	gamma = NEW_STRUCT( Graph );
	gamma->num_meetings = nv;
	gamma->num_components = 0;
	gamma->meeting = NEW_ARRAY( nv, GraphMeeting );

	for( i = 0; i < nv; i++ )
	{
		m = &gamma->meeting[i];
		m->type = (code[i][1]==0) ? Inter : Cross;
		m->num_strands = code[i][0];
		m->strand = NEW_ARRAY( code[i][0], int );
		m->label = NEW_ARRAY( code[i][0], int );
		m->component = NEW_ARRAY(code[i][0], int );
		m->neighbor = NEW_ARRAY(code[i][0], int );
		m->tet = NULL;

		for(j=0;j<code[i][0];j++)
			get_meeting_info( m, i, j, code );
	}

	if (order > -1 )
		create_labels_and_components( gamma );
	else	create_components( gamma );

	return gamma;
}

static void get_meeting_info( GraphMeeting *m, int index, int strand, int **code)
{
	int nbr,s,a,b,i,j;

	s = (code[index][1]==-1) ? (strand + 3) % 4 : strand;

	m->label[strand] = -1;
	m->component[strand] = -1;

	nbr = code[index][s + 2];
	m->neighbor[strand] = nbr; 

	for( i = 0; i < code[index][0]; i++)
	if (	code[index][i+2] == nbr &&
		code[index][((i+code[index][0]-1) % code[index][0]) + 2] != nbr)
		break;

	for( j = 0; j < code[nbr][0]; j++)
	if (	code[nbr][j+2] == index &&
		code[nbr][((j+1) % code[nbr][0]) + 2] != index)
		break;

	if (	code[index][0]==i ||
		code[nbr][0]==j  )
		uFatalError("get_meeting_info", "plantriBatch");

	a = (s-i+code[index][0]) % code[index][0];
	b = (j-a+code[nbr][0]) % code[nbr][0];

	m->strand[strand] = (code[nbr][1]==-1) ? (b+1) % 4 : b;

	return;
}

static void create_labels_and_components( Graph *gamma )
{
	int		i,j,l,s0,s1;
	GraphMeeting	*m,*m0,*m1;

	for(i=0;i<gamma->num_meetings;i++)
		for( j = 0; j < gamma->meeting[i].num_strands; j++)
		{
			gamma->meeting[i].label[j] = -2;
			gamma->meeting[i].component[j] = -1;
		}

	for(i=0,l=0;i<gamma->num_meetings;i++)
	{
		m = &gamma->meeting[i];

		if (m->type == Inter)
		{
			for( j = 0; j < m->num_strands; j++)
			if (m->label[j] == -2 )
			{
				m0 = m;
				s0 = j;

				do{
					m0->component[s0] = gamma->num_components;	
					m0->label[s0] = -1;
					m1 = &gamma->meeting[m0->neighbor[s0]];
					s1 = m0->strand[s0];
	
					if (m1->type == Inter)
					{
						m1->label[s1] = l++;
						break;
					}

					m1->component[s1] = gamma->num_components;
					m1->label[s1] = -1;
					m0 = m1;
					s0 = (s1 + 2) % 4;

				}while( TRUE );
			}
			else	m->component[j] = gamma->num_components;

			gamma->num_components++;
		}
	}


	for(i=0;i<gamma->num_meetings;i++)
	if (gamma->meeting[i].type == Inter )
		for(j=0;j<gamma->meeting[i].num_strands;j++)
		if (gamma->meeting[i].component[j] == -1 )
		uFatalError("create_labels_and_components","plantriBatch");

	label_torus_components( gamma, l );

	return;
}


static void label_torus_components( Graph *gamma, int l )
{
	int i,j,s0,s1;
	GraphMeeting *new,*m0,*m1;
	
	for(i=0;i<gamma->num_meetings;i++)
	if (gamma->meeting[i].type == Cross )
		for(j=0;j<gamma->meeting[i].num_strands;j++)
		if (gamma->meeting[i].component[j] == -1 )
		/* this must belong to a torus component */		
	{
		resize_meeting_array( gamma, gamma->num_meetings + 1 );

		new = &gamma->meeting[gamma->num_meetings];
		m0 = &gamma->meeting[i];
		s0 = j;
		m1 = &gamma->meeting[m0->neighbor[s0]]; 
		s1 = m0->strand[s0];
	
		if (m1->type == Inter)
			uFatalError("label_torus_components","plantriBatch");	
	
		new->strand[0] = s0;
		new->neighbor[0] = i;
		new->label[0] = l++;

		new->strand[1] = s1;
		new->neighbor[1] = m0->neighbor[s0];
		new->label[1] = -1;	

		m0->strand[s0] = 0;
		m0->neighbor[s0] = gamma->num_meetings;

		m1->strand[s1] = 1;
		m1->neighbor[s1] = gamma->num_meetings;

		m0 = new;
		s0 = 0;

		do{
			m0->component[s0] = gamma->num_components;

			m1 = &gamma->meeting[m0->neighbor[s0]];
			s1 = m0->strand[s0];

			m1->component[s1] = gamma->num_components;

			if (m1->type == Inter)
				break;

			s0 = (s1 + 2) % 4;
			m0 = m1;
	 	}while(TRUE);	

		gamma->num_components++;
		gamma->num_meetings++;		
	}


}

static void     resize_meeting_array(
        Graph *gamma,
        int new_array_size)
{
        /*
         * Resize the meeting array.
         * Do NOT change the num_meetings.
         */

 GraphMeeting   *old_array,
                *new_array;
 int            i,
                j;

 if ( new_array_size < gamma->num_meetings )
        uFatalError("resize_meeting_array", "graph_complement");

 old_array = gamma->meeting;
 new_array = NEW_ARRAY( new_array_size, GraphMeeting );

 for ( i = 0; i< gamma->num_meetings; i++ )
        new_array[i] = old_array[i];

 for ( i = gamma->num_meetings; i< new_array_size; i++)
 {
        new_array[i].type = Inter;
        new_array[i].num_strands = 2;
        new_array[i].component = NEW_ARRAY( 2 , int );

        new_array[i].label = NEW_ARRAY( 2 , int );
        for( j = 0; j < 4; j++ )
                new_array[i].label[j] = -1;

        new_array[i].strand = NEW_ARRAY( 2 , int );
        new_array[i].neighbor = NEW_ARRAY( 2, int );
 }

 my_free( old_array );
 gamma->meeting = new_array;
}



static void create_components( Graph *gamma )
{
	int i,j,a,c,p,q,q1,s0,s1;
	GraphMeeting *m,*m0,*m1;

	gamma->num_components = 0;

	for(i=0;i<gamma->num_meetings;i++)
	{
		for( j = 0; j < gamma->meeting[i].num_strands; j++)
		{
			gamma->meeting[i].label[j] = 0;
			gamma->meeting[i].component[j] = -1;
		}
	}

        for(i=0,a=-2;i<gamma->num_meetings;i++)
        {
                m = &gamma->meeting[i];

                if (m->type == Inter)
                {
                        for( j = 0; j < m->num_strands; j++)
                        if (m->label[j] == 0 )
                        {
                                m0 = m;
                                s0 = j;

                                do{
                                        m0->component[s0] = a;
					m0->label[s0] = -1;

                                        m1 = &gamma->meeting[m0->neighbor[s0]];
                                        s1 = m0->strand[s0];

					m1->component[s1] = a;
					m1->label[s1] = -1;

                                        if (m1->type == Inter)
                                                break;

                                        m0 = m1;
                                        s0 = (s1 + 2) % 4;

                                }while( TRUE );

				a--;
			}
		}
	}

	for(i=0,c=0;i<gamma->num_meetings;i++)
	{
		if (gamma->meeting[i].type == Inter && gamma->meeting[i].component[0] < -1 )
		{
			a = gamma->meeting[i].component[0];

			do{

				for(p = 0; p < gamma->num_meetings; p++)
					for(q=0;q<gamma->meeting[p].num_strands;q++)
					if (gamma->meeting[p].component[q]==a)
						gamma->meeting[p].component[q] = c;

				for(p = 0, a=0; a==0 && p < gamma->num_meetings; p++)
				if (gamma->meeting[p].type == Inter )
					for(q=0;a==0 && q<gamma->meeting[p].num_strands;q++)
					if (gamma->meeting[p].component[q]==c)
						for(q1=0;a==0 && q1<gamma->meeting[p].num_strands;q1++)
						if (gamma->meeting[p].component[q1] < -1 )
							a = gamma->meeting[p].component[q1];
			}
			while(a<0);

			c++;
		}
	}


        for(i=0;i<gamma->num_meetings;i++)
        {
                m = &gamma->meeting[i];

                if (m->type == Cross)
                {
                        for( j = 0; j < m->num_strands; j++)
                        if (m->label[j] == 0 )
                        {
                                m0 = m;
                                s0 = j;

                                do{
                                        m0->component[s0] = c;
                                        m0->label[s0] = -1;

                                        m1 = &gamma->meeting[m0->neighbor[s0]];
                                        s1 = m0->strand[s0];

                                        m1->component[s1] = c;
                                        m1->label[s1] = -1;

                                        m0 = m1;
                                        s0 = (s1 + 2) % 4;

					if (m0->component[s0] == c)
						break;

                                }while( TRUE );

				c++;
                        }
                }
        }

	gamma->num_components = c;

        for(i=0;i<gamma->num_meetings;i++)
        for(j=0;j<gamma->meeting[i].num_strands;j++)
        if ( gamma->meeting[i].component[j] < 0 )
                uFatalError("create_labels_and_components","plantriBatch");

	return;
}


