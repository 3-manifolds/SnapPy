#include <stdio.h>
#include "kernel.h"

static Triangulation	*read_italian_triangulation( FILE * , double ***);
static void		find_italian_triangulation( FILE *, int );

extern Triangulation	*read_italian_census( int n, double ***italian_angles, char *s )
{
	FILE		*fp;
	Triangulation	*manifold;

	fp = NULL;

	if ( (fp=fopen( s,"r")) == NULL )
	{
		fprintf(stderr,"Unable to open Italian Census.\n");
		uFatalError("read_italian_census","italian_census");
	}

	find_italian_triangulation( fp, n );

	manifold = read_italian_triangulation( fp, italian_angles );

	return manifold;
}

static void		find_italian_triangulation(
	FILE	*fp,
	int	n )
{
	int	cur;
	char	s[200];

	cur = 0;

	while( cur != n )
		if ( NULL == fgets( s, 200, fp ) )
			uFatalError("find_italian_triangulation","italian_census");
		else if ( s[0] == '-')
			cur++;	
}


static Triangulation	*read_italian_triangulation(
	FILE	*fp,
	double	***italian_angles )
{
	Triangulation	*manifold;
	Tetrahedron	**tet;
	int		num_tet,
			nbr,
			i,
			j,
			k,
			gluing[4],
			italian_2_snappea_edge[6] = {1,0,2,4,5,3};
	char		s[200];
	manifold = NULL;


	manifold = NEW_STRUCT( Triangulation );
	initialize_triangulation( manifold );

	fscanf(fp, "%d", &num_tet);
	tet = NEW_ARRAY( num_tet, Tetrahedron *);

	for( i=0; i<num_tet; i++)
	{
		tet[i] = NEW_STRUCT( Tetrahedron );
		initialize_tetrahedron( tet[i] );
		INSERT_BEFORE( tet[i], &manifold->tet_list_end );
	}
        manifold->num_tetrahedra = num_tet;
	(*italian_angles) = NEW_ARRAY( num_tet, double *);

	for( i=0; i<num_tet; i++ )
	{
		for(j=0;j<4;j++)
		{
			fscanf( fp, "%d", &nbr );
			if ( nbr < 0 || nbr >= num_tet )
				uFatalError("read_italian_triagulation","italian_census");
			tet[i]->neighbor[j] = tet[nbr];
		}

		for(j=0;j<4;j++)
		{
			for(k=0;k<4;k++)
			{
				fscanf( fp, "%1d", &gluing[k] );
				if (gluing[k] < 0 || gluing[k] > 3 )
					uFatalError("read_italian_triagulation","italian_census");
			}
			tet[i]->gluing[j] = CREATE_PERMUTATION(
						0,gluing[0],
						1,gluing[1],
						2,gluing[2],
						3,gluing[3]);
		}
	}

	while( fgetc(fp) != ':');

	for(i=0;i<num_tet;i++)
	{
		(*italian_angles)[i] = NEW_ARRAY( 6, double );
		while( fgetc(fp) != ':');
		for(j=0;j<6;j++)
			fscanf(fp, "%lf", &(*italian_angles)[i][italian_2_snappea_edge[j]] );

	}

	while( TRUE )
		if ( NULL == fgets( s, 200, fp ) )
			uFatalError("read_italian_triangulation","italian_census");
                else if ( s[0] == '-')
                        break;
/*		else puts(s); ??? out puts the rest of the data */
	

        create_edge_classes( manifold );
        orient_edge_classes( manifold );

        create_cusps( manifold );

        mark_fake_cusps( manifold );

        count_cusps( manifold );

	return manifold;

}

