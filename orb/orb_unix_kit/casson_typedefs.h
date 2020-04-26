#ifndef _casson_typedefs_
#define _casson_typedefs_

#include "kernel_typedefs.h"

// was gui/casson.h

// Probably not used anymore
#define LN(ch)   (ch=='u') ? 0 : ((ch=='v') ? 1 : ((ch=='w') ? 2 : 3))

extern const int vertex_at_faces[4][4];

typedef struct CassonFormat CassonFormat;
typedef struct EdgeInfo EdgeInfo;
typedef struct TetEdgeInfo TetEdgeInfo;

struct CassonFormat
{
	SolutionType	type;
  //	bool		vertices_known;
        Boolean         vertices_known;
	int		num_tet;
	EdgeInfo	*head;
};

struct EdgeInfo
{
	int		index,
			one_vertex,
			other_vertex,
			singular_index;
	double 		singular_order,
			e_inner_product,
			v_inner_product1,
			v_inner_product2;

	TetEdgeInfo	*head;
	EdgeInfo	*prev,
			*next;
};

struct TetEdgeInfo
{
	int		tet_index,f1,f2, curves[8];
	double		dihedral_angle;
	TetEdgeInfo	*prev,
			*next;
};


#endif
