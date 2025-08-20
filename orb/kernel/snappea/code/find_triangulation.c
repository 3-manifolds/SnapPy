#include "kernel.h"
#include <stdio.h>

typedef struct container Container;

struct container
{
	Triangulation	*manifold;
	Container	*next,
			*prev;
};

static Boolean new_triangulation( Container *head, Triangulation *manifold );

extern int find_triangulation(
		Triangulation *manifold,
		Triangulation ***list,
		int num_triangulations,
		int max_randomizations )
{
	Triangulation	*copy0,
			*copy1;
	Container	*head = NULL,
			*cur = NULL,
			*tail = NULL;
	int		i,
			t;
	Boolean		have_canon;

	if ((*list)!=NULL || manifold == NULL )
		uFatalError("find_triangulation","find_triangulation");

	copy_triangulation( manifold, &copy0 );

	have_canon = FALSE;
	i = 0;
	t = 0;
	do{
		t++;
		cur = NEW_STRUCT( Container );
		cur->manifold = copy0;
		cur->prev = tail;
		cur->next = NULL;

		if (head==NULL)
		{
			head = cur;
			tail = cur;
		}
		else
		{
			tail->next = cur;
			tail = cur;
		}

		if (t==num_triangulations)
			break;

		copy_triangulation( copy0, &copy1 );

		copy0 = copy1;

		for( ; i < max_randomizations; i++ )
		{
			if (have_canon==FALSE)
			{
				if (canonize(copy0)==func_failed)
					randomize_triangulation( copy0 );
				have_canon = TRUE;
			}
			else	randomize_triangulation( copy0 );

			if (new_triangulation(head,copy0))
				break;
		}
		
	}while(i<max_randomizations);


	(*list) = NEW_ARRAY( t, Triangulation * );	

	i = 0;
	cur = head;
	for( i = 0; i < t; i++ )
	{
		(*list)[i] = head->manifold;
		cur = head->next;
		my_free( head );
		head = cur;
	}

	return t;
}

static Boolean new_triangulation( Container *head, Triangulation *manifold )
{
	Container *cur = head;

	while(cur!=NULL)
	{
		if (same_triangulation(manifold,cur->manifold,NULL))
			return FALSE;
		cur = cur->next;
	}

	return TRUE;
}

extern void free_triangulation_list( Triangulation **list, int num )
{
	int i;

	for(i=0;i<num;i++)
		free_triangulation( list[i] );
	my_free( list );
}

