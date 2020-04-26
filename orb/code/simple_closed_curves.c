#include "kernel.h"
#include <stdlib.h>
#include <stdio.h>

struct CurveMark{
	FundamentalEdge	*edge;
	CurveMark	*next,
			*prev,
			*forwards,
			*backwards,
			*mate;
	Boolean		is_fake,
			curve_incident;

	CurveSegment	*curve;
	int		index; /* for debugging */
};





static void		set_up_fake_marks( FundamentalEdge *);
static CurveSegment	*is_scc( FundamentalEdge *, int, int *);
static Boolean		find_next_mark( CurveMark *, CurveMark **, CurveSegment *, CurveSegment **, int, int, int *,Boolean *);	
static void	follow_curve( CurveDirection, CurveMark *, CurveMark **, CurveSegment *, CurveSegment **,CurveMark *,Boolean *);
static Boolean	follow_domain( CurveMark *,CurveMark **, CurveSegment *, CurveSegment **, int, CurveMark *,Boolean *);
static void		link_curve_fields( CurveSegment *);
static Permutation	create_gluing( VertexIndex, VertexIndex );
static Boolean		swap_sides( CurveMark *,int , int *, int);
static void		free_marks( CurveMark *);
static int		reduce_word( int *, int );

int mark_index;

extern CurveSegment	*curve_is_simple(
	FundamentalEdge	*domain,
	int curve_length,
	int *curve )
{
 CurveSegment *curve_start;

 mark_index = 1;

 if (domain == NULL)
	uFatalError("curve_is_simple", "simple_closed_curves");

 if ((curve_length = reduce_word( curve, curve_length )) == 0)
	return NULL;

 set_up_fake_marks( domain );

 if ((curve_start = is_scc( domain, curve_length, curve))==NULL)
	return NULL;
 else
 {
	link_curve_fields(curve_start);
	return curve_start;
 }

}



static void	set_up_fake_marks( FundamentalEdge *domain )
{
 CurveMark	*fake1,
		*fake2;
 FundamentalEdge	*domain1;

 domain1 = domain;

 domain1->inner_curve = NULL;

 fake1 = NEW_STRUCT( CurveMark );
 fake2 = NEW_STRUCT( CurveMark );

 fake1->index = 0;
 fake2->index = 0;


 fake1->is_fake = TRUE;
 fake2->is_fake = TRUE;

 fake1->curve = NULL;
 fake2->curve = NULL;

 fake1->curve_incident=FALSE;
 fake2->curve_incident=FALSE;

 fake1->edge = domain1;
 fake2->edge = domain1;

 fake1->next = fake2;
 fake1->prev = NULL;
 fake1->forwards = NULL;
 fake1->backwards = NULL;
 fake1->mate = NULL;

 fake2->next = NULL;
 fake2->prev = fake1;
 fake2->forwards = NULL;
 fake2->backwards = NULL;
 fake2->mate = NULL; 

 domain1->head_mark = fake1;
 domain1->tail_mark = fake2;

 for( domain1 = domain1->next; domain1 != domain; domain1 = domain1->next)
 { 
	domain1->inner_curve = NULL;

	fake1 = NEW_STRUCT( CurveMark );
 	fake2 = NEW_STRUCT( CurveMark );

	fake1->index = 0;
	fake2->index = 0;

	fake1->is_fake = TRUE;
	fake2->is_fake = TRUE;

	fake1->curve = NULL;
	fake2->curve = NULL;

 	fake1->curve_incident=FALSE;
 	fake2->curve_incident=FALSE;

	fake1->edge = domain1;
 	fake2->edge = domain1;

	fake1->next = fake2;
	fake1->prev = NULL;
	fake1->mate = NULL;
	fake1->forwards = NULL;
	fake1->backwards = NULL;

	fake2->prev = fake1;
	fake2->next = NULL;
	fake2->mate = NULL;
	fake2->forwards = NULL;
	fake2->backwards = NULL;

	domain1->head_mark = fake1;
	domain1->tail_mark = fake2;
 }


 do{
	domain1->head_mark->prev = domain1->prev->tail_mark;
	domain1->tail_mark->next = domain1->next->head_mark;

	if (domain1->gluing_parity == orientation_reversing)
	{
		domain1->head_mark->mate = domain1->mate->head_mark;
		domain1->tail_mark->mate = domain1->mate->tail_mark;
	}
	else
	{
		domain1->head_mark->mate = domain1->mate->tail_mark;
		domain1->tail_mark->mate = domain1->mate->head_mark;
	}	

	domain1 = domain1->next;
 }
 while ( domain1 != domain );

}





static CurveSegment	*is_scc(
 FundamentalEdge	*domain,
 int			curve_length,
 int			*curve)
{

 CurveMark		*mark1,
			*mark2,
			*new_mark1,
			*new_mark2;
 CurveSegment		*prev_curve,
			*cur_curve;
 FundamentalEdge	*domain1;
 int			i;
 Boolean		curve_used;

 for( domain1 = domain;
	domain1->index != curve[0];
	domain1 = domain1->next);

 mark1 = NEW_STRUCT(CurveMark);
 mark2 = NEW_STRUCT(CurveMark);

 mark1->index = mark_index;
 mark2->index = mark_index++;

 mark1->is_fake = FALSE;
 mark2->is_fake = FALSE;

 mark1->curve = NULL;
 mark2->curve = NULL;

 mark1->curve_incident=FALSE;
 mark2->curve_incident=FALSE;

 mark1->edge = domain1;
 mark2->edge = domain1->mate;

 mark1->forwards = NULL;
 mark2->forwards = NULL;

 mark1->backwards = NULL;
 mark2->backwards = NULL;

 mark1->mate = mark2;
 mark2->mate = mark1;

 mark1->next = mark1->edge->tail_mark;
 mark1->prev = mark1->edge->head_mark;
 mark1->edge->head_mark->next = mark1;
 mark1->edge->tail_mark->prev = mark1;

 mark2->next = mark2->edge->tail_mark;
 mark2->prev = mark2->edge->head_mark;
 mark2->edge->head_mark->next = mark2;
 mark2->edge->tail_mark->prev = mark2;

 mark1 = mark2;
 prev_curve = NULL;


 for( i = 1; i<curve_length; i++)
 {
	curve_used = FALSE;

	if (!find_next_mark(mark1,&mark2,prev_curve,&cur_curve,i,curve_length,curve,&curve_used))
	{
		free_curve(cur_curve);
		free_marks(mark1);

		return NULL;
	}
	new_mark1 = NEW_STRUCT(CurveMark);
	new_mark2 = NEW_STRUCT(CurveMark);


	new_mark1->index = mark_index;
	new_mark2->index = mark_index++;

 	new_mark1->curve_incident=FALSE;
 	new_mark2->curve_incident=FALSE;

	new_mark1->curve = cur_curve;

	if (curve_used)
	{
		new_mark1->curve_incident = TRUE;
		cur_curve->mark[curve_forwards] = new_mark1;
	}
	new_mark1->is_fake = FALSE;
 	new_mark2->is_fake = FALSE;

	new_mark1->forwards = NULL;
	new_mark1->backwards = mark1;		/* haven't set marks at the end of curves */
	new_mark1->mate = new_mark2;
	new_mark1->edge = mark2->edge;

	new_mark2->forwards = NULL;
	new_mark2->backwards = NULL;
	new_mark2->curve = NULL;
	new_mark2->mate = new_mark1;
	new_mark2->edge = mark2->edge->mate;

	mark1->forwards = new_mark1;

	new_mark1->next = mark2->next;
	new_mark1->prev = mark2;
	mark2->next = new_mark1;
	new_mark1->next->prev = new_mark1;

	if (mark2->edge->gluing_parity == orientation_reversing )
	{
		new_mark2->next = mark2->mate->next;
		new_mark2->prev = mark2->mate;
		mark2->mate->next = new_mark2;
		new_mark2->next->prev = new_mark2;
	}
	else
	{
		new_mark2->prev = mark2->mate->prev;
		new_mark2->next = mark2->mate;
		mark2->mate->prev = new_mark2;
		new_mark2->prev->next = new_mark2;
	}

	mark1 = new_mark2;
	prev_curve = cur_curve;
 }

 curve_used = FALSE;

 if (!find_next_mark(mark1,&mark2,prev_curve,&cur_curve,0,curve_length,curve,&curve_used)
	|| mark2->next->is_fake == TRUE 
	|| mark2->next->forwards != NULL
	|| mark2->next->backwards != NULL )
 {
	free_curve( cur_curve );
	free_marks( mark1 );

	return NULL;
 }

 mark1->forwards = mark2->next;
 mark2->next->backwards = mark1;

 mark2->next->curve = cur_curve;

 if (curve_used)
 {
	mark2->next->curve_incident = TRUE;	
 	cur_curve->mark[curve_forwards] = mark2->next;
 }

 cur_curve->next[curve_forwards] = mark2->next->mate->curve;
 mark2->next->mate->curve->next[curve_backwards] = cur_curve;

 free_marks(mark1);

 return cur_curve; /* free the marks */
}



static Boolean	find_next_mark(
	CurveMark	*mark1,
	CurveMark	**mark2,
	CurveSegment	*prev_curve,
	CurveSegment	**cur_curve,
	int		i,
	int		curve_length,
	int		*curve,
	Boolean		*curve_used)
{
 CurveMark		*cur_mark,
			*next_mark;
 CurveSegment		*next_curve;

 cur_mark = mark1;
 *curve_used = FALSE;

 while(follow_domain(cur_mark,&next_mark,prev_curve,&next_curve,curve[i],mark1,curve_used))
 {

	cur_mark = next_mark;
	prev_curve = next_curve;

 	if (cur_mark->forwards != NULL)
		follow_curve(curve_forwards,cur_mark,&next_mark,prev_curve,&next_curve,mark1,curve_used);
	else if (cur_mark->backwards != NULL)
		follow_curve(curve_backwards,cur_mark,&next_mark,prev_curve,&next_curve,mark1,curve_used);
	else
		uFatalError("find_next_mark","simple_closed_curves");

	cur_mark = next_mark;
	prev_curve = next_curve;
/*fprintf(stderr, "fc2: m-%d e-%d\n",cur_mark->index,cur_mark->edge->index);*/
 }	

 if (next_mark == mark1)
	return FALSE;



 if ( i != 0 && !next_mark->next->is_fake
	&& next_mark->next->forwards == NULL
	&& next_mark->next->backwards == NULL
	&& swap_sides( next_mark, i, curve, curve_length ))
	next_mark = next_mark->next;


 *mark2 = next_mark; 
 *cur_curve = next_curve;

 return TRUE;
}



static void	follow_curve(
	CurveDirection	direction,
	CurveMark	*cur_mark,
	CurveMark	**next_mark,
	CurveSegment	*cur_curve,
	CurveSegment	**next_curve,
	CurveMark	*mark1,
	Boolean		*curve_used)
{
	CurveSegment	*curve,
			*curve1,
			*new_curve;
	Boolean		First_Time;

/*fprintf(stderr, "fc1: m-%d ,e-%d\n",cur_mark->index,cur_mark->edge->index);*/

	curve = cur_curve;
	curve1 = cur_mark->curve;

	First_Time = TRUE;

	if (!cur_mark->curve_incident)
	{
		*next_mark = (direction == curve_forwards) ?
				cur_mark->forwards : cur_mark->backwards;

		*next_curve = cur_curve;

		return;
	}



	do{

		if (First_Time)
			First_Time = FALSE;
		else
		{
			curve = new_curve;
			curve1 = curve1->next[direction];
		};

		new_curve = NEW_STRUCT(CurveSegment);

		new_curve->forwards = ((direction == curve_forwards && curve1->forwards )
					 || (direction == curve_backwards && !curve1->forwards )) ?
						TRUE : FALSE;

		new_curve->mark[curve_forwards] = NULL;
		new_curve->mark[curve_backwards] = NULL;

		new_curve->edge = curve1->edge;

		if (*curve_used == FALSE)
		{
			mark1->curve = new_curve;
			mark1->curve_incident = TRUE;

			new_curve->mark[curve_backwards] = mark1;
			*curve_used = TRUE;
		}

		if (new_curve->forwards == curve1->forwards )
		{
			new_curve->next[curve_left] = curve1->next[curve_left];
			new_curve->next[curve_right] = curve1;

			curve1->next[curve_left] = new_curve;
		}
		else
		{
			new_curve->next[curve_left] = curve1->next[curve_right];
			new_curve->next[curve_right] = curve1;

			curve1->next[curve_right] = new_curve;
		}


		new_curve->next[curve_forwards] = NULL;
		new_curve->next[curve_backwards] = curve;

		if (new_curve->next[curve_left] != NULL)
		{
			if (new_curve->next[curve_left]->forwards == new_curve->forwards)
				new_curve->next[curve_left]->next[curve_right] = new_curve;
			else
				new_curve->next[curve_left]->next[curve_left] = new_curve;
		}
		else if (!new_curve->forwards)
			new_curve->edge->inner_curve = new_curve;


		if (curve != NULL)
			curve->next[curve_forwards] = new_curve;


	}while( curve1->mark[direction] == NULL);

	*next_mark = curve1->mark[direction];
	*next_curve = new_curve;

	return;
}


static void	link_curve_fields( CurveSegment *curve )
{

 CurveSegment *cur_curve;

 cur_curve = curve;


 do{
	if (cur_curve->forwards)
	{
		cur_curve->left_vertex = cur_curve->edge->vertex;
		cur_curve->left_face = cur_curve->edge->face;
		cur_curve->left_tet = cur_curve->edge->tet;
		cur_curve->left_gluing = create_gluing(cur_curve->left_vertex,cur_curve->left_face);

		cur_curve->right_vertex = cur_curve->edge->mate->vertex;
		cur_curve->right_face = cur_curve->edge->mate->face;
		cur_curve->right_tet = cur_curve->edge->mate->tet;
		cur_curve->right_gluing = create_gluing(cur_curve->right_face,cur_curve->right_vertex);

		if ( cur_curve->next[curve_right] == NULL)
			cur_curve->next[curve_right] = cur_curve->edge->mate->inner_curve;
	}
	else
	{
		cur_curve->right_vertex = cur_curve->edge->vertex;
		cur_curve->right_face = cur_curve->edge->face;
		cur_curve->right_tet = cur_curve->edge->tet;
		cur_curve->right_gluing = create_gluing(cur_curve->right_face,cur_curve->right_vertex); /* ??? */

		cur_curve->left_vertex = cur_curve->edge->mate->vertex;
		cur_curve->left_face = cur_curve->edge->mate->face;
		cur_curve->left_tet = cur_curve->edge->mate->tet;
		cur_curve->left_gluing = create_gluing(cur_curve->left_vertex,cur_curve->left_face); /* ??? */

		if ( cur_curve->next[curve_left] == NULL)
			cur_curve->next[curve_left] = cur_curve->edge->mate->inner_curve;
	}

	cur_curve = cur_curve->next[curve_forwards];

 }while( cur_curve != curve);

 /* what about non-orientable cusps?? */

}

static Boolean follow_domain(
	CurveMark	*cur_mark,
	CurveMark 	**next_mark,
	CurveSegment 	*cur_curve,
	CurveSegment 	**next_curve,
	int 		edge_index,
	CurveMark 	*mark1,
	Boolean 	*curve_used)
{
 CurveMark	*mark;
 CurveSegment	*curve,
		*prev_curve;

/*fprintf(stderr,"fd\n");*/

 mark = cur_mark;
 prev_curve = cur_curve;

 while( TRUE )
 {
/*fprintf(stderr, "mark->index = %d, mark->edge->index = %d\n",mark->index,mark->edge->index);*/

	if (mark->edge->index == edge_index)
	{
		*next_curve = prev_curve;
		*next_mark = mark;
		return FALSE;
	}

	if (!mark->next->is_fake)
	{
		if (mark->next == mark1)
		{
			*next_curve = prev_curve;
			*next_mark = mark1;
			return FALSE;
		}

		if (mark->next->forwards != NULL || mark->next->backwards != NULL)
		{
			*next_curve = prev_curve;
			*next_mark = mark->next;
			return TRUE;
		}
		else
		{
			if (!mark->next->next->is_fake)
			{
				if (mark->next->next == mark1)
				{
					*next_curve = prev_curve;
					*next_mark = mark1;
					return FALSE;
				}
				else
				{
					*next_curve = prev_curve;
					*next_mark = mark->next->next;
					return TRUE;
				}
			}
			else
				mark = mark->next;
		}
	}

	if (mark->edge->index < 0)
	{
		curve = NEW_STRUCT( CurveSegment);

		curve->forwards = TRUE;

		curve->mark[curve_forwards] = NULL;
		curve->mark[curve_backwards] = NULL;

		if (*curve_used == FALSE)
		{
			mark1->curve = curve;
			mark1->curve_incident = TRUE;
			curve->mark[curve_backwards] = mark1;

			*curve_used = TRUE;
		}

		curve->edge = mark->edge;

		curve->next[curve_forwards] = NULL;
		curve->next[curve_left] = curve->edge->inner_curve;
		curve->next[curve_right] = NULL;
		curve->next[curve_backwards] = prev_curve;

		if (curve->edge->inner_curve != NULL)
		{
			if (curve->edge->inner_curve->forwards)
				curve->edge->inner_curve->next[curve_right] = curve;
			else
				curve->edge->inner_curve->next[curve_left] = curve;
		}

		curve->edge->inner_curve = curve;

		if ( prev_curve != NULL)
			prev_curve->next[curve_forwards] = curve;

		prev_curve = curve;

	}

	mark = mark->next->next;

	if (mark->edge->index > 0)
	{
		curve = NEW_STRUCT(CurveSegment);

		curve->forwards = TRUE;

		curve->mark[curve_forwards] = NULL;
		curve->mark[curve_backwards] = NULL;

		if (*curve_used == FALSE)
		{
			mark1->curve = curve;
			mark1->curve_incident = TRUE;
			curve->mark[curve_backwards] = mark1;

			*curve_used = TRUE;
		}

		curve->edge = mark->edge;

		curve->next[curve_forwards] = NULL;
		curve->next[curve_left] = curve->edge->inner_curve;
		curve->next[curve_right] = NULL;
		curve->next[curve_backwards] = prev_curve;

		if (curve->edge->inner_curve != NULL)
		{
			if (curve->edge->inner_curve->forwards)
				curve->edge->inner_curve->next[curve_right] = curve;
			else
				curve->edge->inner_curve->next[curve_left] = curve;
		}

		curve->edge->inner_curve = curve;

		if ( prev_curve != NULL)
			prev_curve->next[curve_forwards] = curve;

		prev_curve = curve;
	}

 }

 uFatalError("follow_domain","simple_closed_curves");

 return FALSE;

}



static Permutation create_gluing( VertexIndex v1, FaceIndex v2)
{
 Permutation	gluing1, gluing2;

 gluing1 = CREATE_PERMUTATION( 0,v2,
				1,v1,
				2,one_face_at_edge[edge_between_vertices[v1][v2]],
				3,other_face_at_edge[edge_between_vertices[v1][v2]]);
 gluing2 = CREATE_PERMUTATION( 0,v2,
				1,v1,
				2,other_face_at_edge[edge_between_vertices[v1][v2]],
				3,one_face_at_edge[edge_between_vertices[v1][v2]]);

 if ( parity[gluing1] == orientation_reversing && parity[gluing2] == orientation_reversing )
	uFatalError("create_gluing","simple_closed_curves");

 return (parity[gluing1] == orientation_preserving) ? gluing1 : gluing2;

}


static Boolean swap_sides( CurveMark *mark, int index, int *curve, int curve_length)
{
 int		i,
		edge_index;
 CurveMark	*cur_mark,
		*start,
		*end,
		*mark1;

 start = mark->next;
 cur_mark = start;
 i = index;

 while( i<curve_length &&
	cur_mark->edge->index == curve[i])
 {
	i++;

	cur_mark = cur_mark->mate;

	if (cur_mark->forwards == NULL)
	{
		end = cur_mark;
		cur_mark = start;
	}
	else
		cur_mark = cur_mark->forwards;
 }


 if (i == curve_length)
	edge_index = curve[0];
 else
	edge_index = curve[i];

 mark1 = cur_mark;

 
 if (cur_mark == start)
	cur_mark = end;
 else
	cur_mark = cur_mark->backwards;

 if (curve_length == i)
 {
	for(;	cur_mark != mark1 &&
		cur_mark != start;
		cur_mark = cur_mark->prev);

	if (cur_mark != mark1)
	{
		/*fprintf(stderr,"swap %d\n",index);*/
		return TRUE;

	}	
 }
 else
 {
 	for(;	cur_mark->edge->index != edge_index
		&& cur_mark != mark1;
		cur_mark = cur_mark->prev );

 	if (cur_mark != mark1)
 	{
		/*fprintf(stderr,"swap %d\n",index);*/
		return TRUE;
 	}
 }

 /*fprintf(stderr,"don't swap %d\n",index);*/
 return FALSE;
}





static void	free_marks( CurveMark *mark )
{
 CurveMark	*next,
		*end;

 end = mark->prev;

 while (mark != end)
 {
	next = mark->next;	
	my_free( mark );
	mark = next;
 }	

 my_free(end);

}


extern void	free_curve( CurveSegment *curve)
{
 CurveSegment *cur;

 while (curve != NULL)
 {
	cur = curve->next[curve_backwards];
	my_free(curve);
	curve = cur;
 }


}


static int reduce_word( int *curve, int curve_length)
{
 Boolean not_finished;
 int	problem_letter,
	i;


 if (curve_length == 0 || curve_length == 1)
	return curve_length;



 not_finished = FALSE;


 for (i=0;i<curve_length-1;i++)
	if (curve[i] == -curve[i+1])
	{
		problem_letter = i;
		not_finished = TRUE;
	}

 if (curve[0] == -curve[curve_length-1])
 {
	problem_letter = curve_length-1;
	not_finished = TRUE;
 }


 while (not_finished){

	if (curve_length != 2)
	{
		if ( problem_letter<curve_length-1)
			for( i = problem_letter; i < curve_length - 2; i++)
				curve[i] = curve[i+2];
		else
			for( i = 0; i < curve_length - 1; i++)
				curve[i] = curve[i+1]; 
	}
	
	curve_length -= 2;

	not_finished = FALSE;

	for (i=0;i<curve_length-1;i++)
		if (curve[i] == -curve[i+1])
		{
			curve_length = i;
			not_finished = TRUE;
		}

	if (curve_length > 0 && curve[0] == -curve[curve_length-1])
	{
		problem_letter = curve_length-1;
		not_finished = TRUE;
        }

 };


 return curve_length;
}

