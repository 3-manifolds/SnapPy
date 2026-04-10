#include "kernel.h"

static void compute_one_boundary_relation( FundamentalEdge *, int, int **);
static int      basic_reduction( FundamentalEdge *, int ***, int );
static Boolean	is_length_two_relation( int *);
static int	remove_dead_relations( int ***, int num_rels);
static Boolean	is_length_two_relation( int * );

extern int	find_boundary_relations(
	FundamentalEdge		*domain,
	int			***relations)
{
	FundamentalEdge	*cur;
	int		**temp1,
			**temp2,
			i,
			j;

	temp1 = NULL;
	temp2 = NULL;

	i=0;
	cur = domain;

	do{

		if (!cur->checked_end[0])
		{
			temp2 = NEW_ARRAY( i+1, int *);

			for( j=0; j < i; j++)
				temp2[j] = temp1[j];

			my_free(temp1);

			temp1 = temp2;

			compute_one_boundary_relation( cur, 0, &temp1[i++]);
		}

		if (!cur->checked_end[1])
		{
			temp2 = NEW_ARRAY( i+1, int *);

			for( j=0; j < i; j++)
				temp2[j] = temp1[j];

			my_free(temp1);

			temp1 = temp2;

			compute_one_boundary_relation( cur, 1, &temp1[i++]);
		}
		cur = cur->next;	

	}while( cur != domain );

	(*relations) = temp1;

	cur = domain;

	do{
		cur->simplified_generator = NEW_ARRAY( 2, int);

		cur->simplified_generator[0] = cur->index;
		cur->simplified_generator[1] = 0;

		cur = cur->next;

	}while(cur != domain);


	return basic_reduction(domain, relations, i);	
}


static void	compute_one_boundary_relation(
	FundamentalEdge		*start,
	int			start_end,
	int			**relation)
{

	FundamentalEdge		*cur;
	int			cur_end,
				i,
				j,
				temp[50];

	i = 0;
	cur = start;
	cur_end = start_end;

	while( !cur->checked_end[cur_end] )
	{
		temp[i++] = cur->index;
		cur->checked_end[cur_end] = TRUE;

		if ((cur->gluing_parity == orientation_preserving && cur_end == 1) ||
			( cur->gluing_parity == orientation_reversing && cur_end == 0	))
		{
			cur->mate->checked_end[0] = TRUE;

			cur = cur->mate->prev;
			cur_end = 1;
		}
		else
		{
			cur->mate->checked_end[1] = TRUE;

			cur = cur->mate->next;
			cur_end = 0;
		}
	}

	(*relation) = NEW_ARRAY( i+1, int );

	for ( j = 0; j < i; j++)
		(*relation)[j] = temp[j];
	(*relation)[i] = 0;

}



static int	basic_reduction(
	FundamentalEdge	*domain,
	int		***relations,
	int		num_rels)
{
	int	dead_generator,
		replacement_generator,
		i,
		j,
		k;
	FundamentalEdge	*cur;


	for( i=0; i < num_rels; i++)
	if ( is_length_two_relation( (*relations)[i] ))
	{
		dead_generator = (*relations)[i][0];
		replacement_generator = -(*relations)[i][1];

		(*relations)[i][0] = 0;

		for( j =0; j< num_rels; j++ )
		{
			k = 0;

			while( (*relations)[j][k] != 0)
			{

				if  ( (*relations)[j][k] == dead_generator )
					(*relations)[j][k] = replacement_generator;

				if ( (*relations)[j][k] ==  -dead_generator ) 
					(*relations)[j][k] = -replacement_generator;

				k++;
			}
		}

		cur = domain;

		do{
			if ( cur->index == dead_generator )
				cur->simplified_generator[0] = replacement_generator;

			if (cur->index == -dead_generator )
				cur->simplified_generator[0] = -replacement_generator;

			cur = cur->next;

		}while( cur != domain);

	}

	return remove_dead_relations( relations, num_rels );

}

	
static int	remove_dead_relations(
	int	***relations,
	int	num_rels)
{
	int	**temp,
		i,
		j,
		dead_rels;

	dead_rels = 0;

	for( i=0; i<num_rels; i++)
		if ( (*relations)[i][0] == 0)
			dead_rels++;

	temp = NEW_ARRAY( num_rels - dead_rels, int *);

	j = 0;

	for( i=0; i < num_rels; i++)
		if ( (*relations)[i][0] != 0 )
			temp[j++] = (*relations)[i];
		else
			my_free( (*relations)[i]);

	my_free( *relations );

	(*relations) = temp;

	return num_rels - dead_rels;

}


static Boolean is_length_two_relation( int *relation )
{
	int i;

	i = 0;

	while( relation[i] != 0 )
		i++;

	return (i == 2);

}


