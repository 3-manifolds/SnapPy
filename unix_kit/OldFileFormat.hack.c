/*
 *	OldFileFormat.hack.c
 *
 *	This is some old hacked together code for reading manifolds
 *	in the old file format.  I'm intentionally not providing code to write
 *	triangulations in the old format.  I'm hoping people will change over
 *	to the new one.  When Macintosh SnapPea reads a file in the old format
 *	it puts the contents into a new untitled window.
 *
 *	Technical note:  What makes reading old files hard is the
 *	combination of the way orientations and peripheral curves are handled.
 *	Snappea 1.3.3 didn't orient triangulations the way SnapPea 2.0 does.
 *	Instead, each tetrahedron kept a flag indicating whether its local
 *	orientation agreed with the global orientation or not.  As for the
 *	peripheral curves, snappea 1.x originally didn't handle nonorientable
 *	manifolds at all -- that capability was retrofitted.  It keeps all
 *	peripheral curves on the boundary component itself (not on the double
 *	cover, as SnapPea 2.0 does).  Meridians on Klein bottle cusps can't
 *	be stored correctly in this format (snappea 1.3.3 used an ad hoc
 *	kludge).  Even torus cusps in orientable manifolds are a nuisance
 *	because we have to figure out which curves should go onto which sheet
 *	of the double cover.  The code below simply recomputes the peripheral
 *	curves for nonorientable manifolds from scratch (the original curves
 *	are ignored).  It also fusses around with the peripheral curves on
 *	orientable manifolds to make them conform to the usual orientation
 *	conventions.
 */

/* MC 2013-5-9 fixed scanf formats to stop gcc warnings. */

#include "kernel.h"
#include <stdio.h>

#define DEFAULT_NAME	"untitled"

		FuncResult	read_old_manifold(FILE *fp, Triangulation **manifold);
static	FuncResult	read_the_file(FILE *fp, Triangulation **manifold);
static	void		fix_orientation(Triangulation *manifold);


FuncResult read_old_manifold(
	FILE			*fp,
	Triangulation	**manifold_ptr)
{
	char	keyword[100];
	double  ultimate_dbl, penultimate_dbl;

	/*
	 *	First get the basic information out of the file.
	 */
	if (read_the_file(fp, manifold_ptr) == func_failed)
		return func_failed;

	/*
	 *	Then add the bells and whistles.  Notice that the original
	 *	peripheral curves are preserved for an orientable manifold,
	 *	but replaced for a nonorientable one.
	 */

	/*
	 *	read_the_file() puts all peripheral curves onto the right handed
	 *	sheets, which means they won't match up where right handed
	 *	tetrahedra meets left handed ones.  Fortunately, when the manifold
	 *	is orientable, orient() transfers all curves to the right handed
	 *	sheets.  This has the pleasant side effect of making the peripheral
	 *	curves match up.
	 */
	orient(*manifold_ptr);
	orient_edge_classes(*manifold_ptr);

	if ((*manifold_ptr)->orientability == nonorientable_manifold)
		peripheral_curves(*manifold_ptr);
	else
		fix_orientation(*manifold_ptr);

	find_complete_hyperbolic_structure(*manifold_ptr);

	/*
	 *	peripheral_curves() will have trashed the original curves,
	 *	even on the orientable cusps.  Repair the damage as best we can.
	 */
	if ((*manifold_ptr)->orientability == nonorientable_manifold)
		install_shortest_bases(*manifold_ptr);

	(*manifold_ptr)->CS_value_is_known =
				(3 == fscanf(fp, "%99s %lf %lf",
					keyword,
					&ultimate_dbl,
					&penultimate_dbl)
				 && strcmp(keyword, "CS") == 0);
	if ( (*manifold_ptr)->CS_value_is_known ) {
	  (*manifold_ptr)->CS_value[ultimate] = (Real)ultimate_dbl;
	  (*manifold_ptr)->CS_value[penultimate] = (Real)penultimate_dbl;
	}
	compute_CS_fudge_from_value(*manifold_ptr);

	return func_OK;
}


FuncResult read_the_file(
	FILE			*fp,
	Triangulation	**manifold_ptr)
{
	int				i,
					c,
					cusp_index,
					vertex,
					face,
					edge,
					count,
					nbr_index,
					edge_index,
					digit,
					d;
	Triangulation	*manifold;
	Tetrahedron		**tal,	/* tetrahedron address list */
	                         *tet;
	EdgeClass		**eal;	/* edge class address list	*/
	Cusp			**cal;	/* cusp address list		*/

	*manifold_ptr = NULL;	/* just in case anything goes wrong . . . */

	manifold = NULL;
	tal = NULL;
	eal = NULL;
	cal = NULL;

	manifold = NEW_STRUCT(Triangulation);
	initialize_triangulation(manifold);

	manifold->name = NEW_ARRAY(strlen(DEFAULT_NAME) + 1, char);
	strcpy(manifold->name, DEFAULT_NAME);

	if (fscanf(fp, "%d%d%d%d%*f%*d",
			&manifold->num_tetrahedra,
			&manifold->num_cusps,
			&manifold->num_nonor_cusps,
			&manifold->orientability) != 4)
		goto bail;

	for (i = 0; i < manifold->num_tetrahedra; i++)
	  if ( fscanf(fp, "%*d") != 0)	/* ignore edge class sizes */
	    goto bail;

	manifold->num_or_cusps = manifold->num_cusps - manifold->num_nonor_cusps;

	if (manifold->num_tetrahedra < 1
	 || manifold->orientability < 0 || manifold->orientability > 2
	 || manifold->num_or_cusps + manifold->num_nonor_cusps != manifold->num_cusps)
		goto bail;

	/*
	 *	If the user sends us an old format link projection by
	 *	mistake we'll eventually catch it in the code below.
	 *	But let's try to catch it here, to save possible problems with
	 *	allocating large numbers of Tetrahedra, Cusps and EdgeClasses.
	 */
	if (manifold->num_cusps > 2 * manifold->num_tetrahedra)
		goto bail;

	cal = NEW_ARRAY(manifold->num_cusps, Cusp *);
	for (count = 0; count < manifold->num_cusps; count++)
	{
		cal[count] = NEW_STRUCT(Cusp);
		initialize_cusp(cal[count]);
		cal[count]->index = count;
		INSERT_BEFORE(cal[count], &manifold->cusp_list_end);
		cal[count]->topology = (count < manifold->num_nonor_cusps) ? Klein_cusp : torus_cusp;
	}

	tal = NEW_ARRAY(manifold->num_tetrahedra, Tetrahedron *);
	for (count = 0; count < manifold->num_tetrahedra; count++)
	{
		tal[count] = NEW_STRUCT(Tetrahedron);
		initialize_tetrahedron(tal[count]);
		INSERT_BEFORE(tal[count], &manifold->tet_list_end);
	}

	eal = NEW_ARRAY(manifold->num_tetrahedra, EdgeClass *);
	for (count = 0; count < manifold->num_tetrahedra; count++)
	{
		eal[count] = NEW_STRUCT(EdgeClass);
		initialize_edge_class(eal[count]);
		INSERT_BEFORE(eal[count], &manifold->edge_list_end);
	}

	for (count = 0; count < manifold->num_tetrahedra; count++)
	{
		tet = tal[count];
		tet->index = count;
		for (face = 0; face < 4; face++)
		{
			if (fscanf(fp, "%d", &nbr_index) != 1
			 || nbr_index < 0 || nbr_index >= manifold->num_tetrahedra)
				goto bail;
			tet->neighbor[face] = tal[nbr_index];
		}
		for (face = 0; face < 4; face++)
			for (digit = 4; --digit >= 0; )
			{
				if (fscanf(fp, "%1d", &d) != 1
				 || d < 0 || d > 3)
					goto bail;
				tet->gluing[face] = (tet->gluing[face] << 2) + d;
			}
		for (vertex = 0; vertex < 4; vertex++)
		{
			if (fscanf(fp, "%d", &cusp_index) != 1
			 || cusp_index < 0 || cusp_index >= manifold->num_cusps)
				goto bail;
			tet->cusp[vertex] = cal[cusp_index];
		}
		for (c = 0; c < 2; c++)
			for (vertex = 0; vertex < 4; vertex++)
				for (face = 0; face < 4; face++)
				{
					if (fscanf(fp, "%d", &tet->curve[c][right_handed][vertex][face]) != 1)
						goto bail;
					tet->curve[c][left_handed][vertex][face] = 0;
				}
		for (edge = 0; edge < 6; edge++)
		{
			if (fscanf(fp, "%d", &edge_index) != 1
			 || edge_index < 0 || edge_index >= manifold->num_tetrahedra)
				goto bail;
			tet->edge_class[edge] = eal[edge_index];
			tet->edge_class[edge]->order++;
			tet->edge_class[edge]->incident_tet			= tet;
			tet->edge_class[edge]->incident_edge_index	= edge;
		}
		/*  parse but ignore tet orientation and tet shape */
		if ( fscanf(fp, "%*d") != 0 || fscanf(fp, "%*f%*f") != 0 )
		  goto bail;
	}

	my_free(tal);
	my_free(eal);
	my_free(cal);

	*manifold_ptr = manifold;

	return func_OK;

bail:
	/*
	 *	The file could not be read correctly.
	 *	Clean up and go home.
	 */
	free_triangulation(manifold);
	if (tal != NULL)
		my_free(tal);
	if (eal != NULL)
		my_free(eal);
	if (cal != NULL)
		my_free(cal);
	return func_failed;
}


static void fix_orientation(
	Triangulation *manifold)
{
	int			standard;
	Cusp		*cusp;
	Tetrahedron	*tet;
	int			i, j, k;

	copy_curves_to_scratch(manifold, 0, FALSE);
	copy_curves_to_scratch(manifold, 1, FALSE);
	compute_intersection_numbers(manifold);

	/*
	 *	All the intersection numbers ought to be the same.
	 */
	standard = manifold->cusp_list_begin.next->intersection_number[L][M];
	for (cusp = manifold->cusp_list_begin.next;
		 cusp != &manifold->cusp_list_end;
		 cusp = cusp->next)
		if (cusp->intersection_number[L][M] != standard)
			uFatalError("fix_orientation", "OldFileFormat.hack");

	/*
	 *	If the curves are backwards, reorient the manifold.
	 */
	if (standard == -1)
	{
		reorient(manifold);

		/*
		 *	Now the manifold's fine, but reorient() has reversed the
		 *	meridians, on the assumption that they were correctly
		 *	oriented to begin with, which they weren't.  Fix 'em.
		 */
		for (tet = manifold->tet_list_begin.next;
			 tet != &manifold->tet_list_end;
			 tet = tet->next)
			for (i = 0; i < 2; i++)
				for (j = 0; j < 4; j++)
					for (k = 0; k < 4; k++)
						tet->curve[M][i][j][k] = - tet->curve[M][i][j][k];
	}
}
