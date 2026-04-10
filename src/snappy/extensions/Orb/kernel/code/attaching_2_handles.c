#include "kernel.h"
#include "../headers/subdivide_extra.h"
#include <stdio.h>

static CurveSegment	*adjust_curve( CurveSegment * );
static void		glue_handle_on( Triangulation *, CurveSegment * );
static void		free_cusps( Triangulation *);
static void		free_extra( Triangulation *);



extern Triangulation *attach_2_handle(
	 Triangulation *manifold,
	 CurveSegment *curve,
	 char		*name )
{
 Triangulation	*new_manifold;
 CurveSegment	*new_curve;

 new_manifold = subdivide( manifold, name, FALSE ); 

 new_curve = adjust_curve( curve );

 free_extra( manifold );

 glue_handle_on( new_manifold, new_curve );

 return new_manifold;
}




static CurveSegment	*adjust_curve(
	CurveSegment *curve_start)
{

 CurveSegment 	*cur;

 cur = curve_start;

 do{
	cur->other = NEW_STRUCT( CurveSegment);
	cur->other->other = cur;	
	cur = cur->next[curve_forwards];
 }
 while( cur != curve_start ) ;



 do{
	cur->other->next[curve_forwards] = cur->next[curve_forwards]->other;
 	cur->other->next[curve_backwards] = cur->next[curve_backwards]->other;
 	cur->other->next[curve_left] = (cur->next[curve_left] == NULL) ? NULL : cur->next[curve_left]->other;
 	cur->other->next[curve_right] = (cur->next[curve_right] == NULL) ? NULL : cur->next[curve_right]->other;

 	cur->other->left_vertex = cur->left_vertex;
 	cur->other->right_vertex = cur->right_vertex;

 	cur->other->left_face = cur->left_face;
	cur->other->right_face = cur->right_face;

	cur->other->left_gluing = cur->left_gluing;
	cur->other->right_gluing = cur->right_gluing;

	cur->other->left_tet = cur->left_tet->extra->outer_vertex_tet[cur->left_vertex];
	cur->other->right_tet = cur->right_tet->extra->outer_vertex_tet[cur->right_vertex];

 	cur = cur->next[curve_forwards];
 }
 while( cur != curve_start );



 return cur->other;
}




static void free_extra(
        Triangulation   *manifold)
{
 Tetrahedron     *tet;

 for (tet = manifold->tet_list_begin.next;
	tet != &manifold->tet_list_end;
	tet = tet->next)
	 {
	my_free(tet->extra);

	tet->extra = NULL;
	}
}



static void	glue_handle_on( Triangulation *manifold, CurveSegment *curve_start )
{
 CurveSegment	*cur;
 Permutation	gluing1,
		gluing2;


 gluing1 = CREATE_PERMUTATION( 0, 0, 1, 1, 2, 3, 3, 2);
 gluing2 = CREATE_PERMUTATION( 0, 1, 1, 0, 2, 2, 3, 3);

 cur = curve_start;

 do{
 	cur->tet = NEW_STRUCT( Tetrahedron );

	initialize_tetrahedron( cur->tet );
        INSERT_BEFORE( cur->tet,
		&manifold->tet_list_end);

	manifold->num_tetrahedra++;

 	cur = cur->next[curve_forwards];
 }
 while(cur != curve_start);


 do{
	cur->tet->neighbor[2] = cur->next[curve_forwards]->tet;
 	cur->tet->gluing[2] = gluing1;

	cur->tet->neighbor[3] = cur->next[curve_backwards]->tet;
	cur->tet->gluing[3] = gluing1;

 	if (cur->next[curve_left] == NULL)
 	{
		cur->tet->neighbor[0] = cur->left_tet;
 		cur->tet->gluing[0] = cur->left_gluing;
 		cur->left_tet->neighbor[cur->left_face] = cur->tet;
		cur->left_tet->gluing[cur->left_face] = inverse_permutation[cur->left_gluing];;
 	}
 	else
 	{
		cur->tet->neighbor[0] = cur->next[curve_left]->tet;	
		cur->tet->gluing[0] = ( cur->next[curve_left]->next[curve_right] == cur ) ?
					gluing2 : gluing1 ;
	}

	if (cur->next[curve_right] == NULL)
	{
		cur->tet->neighbor[1] = cur->right_tet;
		cur->tet->gluing[1] = cur->right_gluing;
		cur->right_tet->neighbor[cur->right_face] = cur->tet;
		cur->right_tet->gluing[cur->right_face] = inverse_permutation[ cur->right_gluing ];
	}
	else 
 	{
		cur->tet->neighbor[1] = cur->next[curve_right]->tet;
		cur->tet->gluing[1] = (cur->next[curve_right]->next[curve_left] == cur) ?
					gluing2 : gluing1 ;
	}

	cur = cur->next[curve_forwards];
 }
 while( cur != curve_start); 

 replace_edge_classes( manifold );
 orient_edge_classes( manifold );

 free_cusps( manifold );

 create_cusps( manifold );

 mark_fake_cusps( manifold );

 remove_finite_vertices( manifold );

 peripheral_curves( manifold );

 count_cusps( manifold );


 return;
}


static void free_cusps(
	Triangulation   *manifold)
{
 Cusp		*dead_cusp;
 Tetrahedron	*tet;
 int		i;

 while (manifold->cusp_list_begin.next != &manifold->cusp_list_end)
 {
	dead_cusp = manifold->cusp_list_begin.next;
	REMOVE_NODE(dead_cusp);
	my_free(dead_cusp);
 }

 for ( tet = manifold->tet_list_begin.next;
	tet != &manifold->tet_list_end;
	tet = tet->next )
	for ( i = 0; i < 4; i++)
 		tet->cusp[i] = NULL;

 manifold->num_cusps = 0;
 manifold->num_nonor_cusps = 0;
 manifold->num_or_cusps = 0;

}

