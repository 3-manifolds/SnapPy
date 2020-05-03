
typedef int CurveDirection;
enum
{
        curve_forwards = 0,
        curve_backwards = 1,
        curve_left = 2,
	curve_right = 3 
};

struct CurveSegment
{
 Boolean	forwards;

 VertexIndex 	left_vertex,
 		right_vertex;

 FaceIndex 	left_face,
 		right_face;

 CurveSegment	*next[4],
 		*other;

 Tetrahedron 	*tet,
 		*left_tet,
 		*right_tet;

 Permutation	left_gluing,
 		right_gluing;

 CurveMark	*mark[2];

 FundamentalEdge	*edge;
};

extern CurveSegment     *curve_is_simple( FundamentalEdge *, int, int *);
extern Triangulation	*attach_2_handle( Triangulation *, CurveSegment *, char *);
extern void		free_curve( CurveSegment * );

