/*
 *	unix_file_io.c
 *
 *	This hacked together file allows unix-style programs
 *	to read and save Triangulations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unix_file_io.h"

#define READ_OLD_FILE_FORMAT 0

#define FALSE	0
#define TRUE	1

/*
 *	gcc complains about the
 *
 *		use of `l' length character with `f' type character
 *
 *	in fprintf() calls.  Presumably it considers the 'l' unnecessary
 *	because even floats would undergo default promotion to doubles
 *	in the function call (see section A7.3.2 in Appendix A of K&R 2nd ed.).
 *	Therefore I've changed "%lf" to "%f" in all fprintf() calls.
 *	If this makes trouble on your system, change it back, and please
 *	let me know (weeks@northnet.org).
 */


static TriangulationData	*ReadNewFileFormat(FILE *fp);
static void my_fgets( char *string, int max, FILE *fp );
#if READ_OLD_FILE_FORMAT
extern FuncResult			read_old_manifold(FILE *fp, Triangulation **manifold);
#endif


Triangulation *get_triangulation(
	char	*file_name)
{
	FILE			*fp;
	Boolean			theNewFormat;
	Triangulation	*manifold;

	/*
	 *	If the file_name is nonempty, read the file.
	 *	If the file_name is empty, read from stdin.
	 */
	if (strlen(file_name) > 0)
	{
		fp = fopen(file_name, "r");
		if (fp == NULL)
			return NULL;

		/*
		 *	Take a peek at the first line to see whether this is
		 *	the new file format or the old one.
		 */
		theNewFormat = (getc(fp) == '%');
		rewind(fp);
	}
	else
	{
		fp = stdin;
		theNewFormat = TRUE;  /* read only the new format from stdin */
	}

	if (theNewFormat == TRUE)
	{
		TriangulationData	*theTriangulationData;

		theTriangulationData = ReadNewFileFormat(fp);

		if (theTriangulationData==NULL)
		{
			if (fp != stdin)
				fclose(fp);

			return NULL;
		}
	
		data_to_triangulation(theTriangulationData, &manifold);
	
		free(theTriangulationData->name);
		free(theTriangulationData->cusp_data);
		free(theTriangulationData->tetrahedron_data);
		free(theTriangulationData);
	}
	else
	{
#if READ_OLD_FILE_FORMAT
		read_old_manifold(fp, &manifold);
#else
		fprintf(stderr, "The manifold is in the old file format.\n");
		fprintf(stderr, "I recommend converting it to the new format.\n");
		fprintf(stderr, "If absolutely necessary, I can provide code for reading the old format.\n");
		fprintf(stderr, "Questions?  Contact me at weeks@northnet.org.\n");
		uFatalError("get_triangulation", "unix file io");
#endif
	}

	if (fp != stdin)
		fclose(fp);

	return manifold;
}


static TriangulationData *ReadNewFileFormat(
	FILE	*fp)
{
	char				theScratchString[100];
	TriangulationData	*theTriangulationData;
	int					theTotalNumCusps,
						i,
						j,
						k,
						v,
						f;

	/*
	 *	Read and ignore the header (% Triangulation).
	 */

	my_fgets(theScratchString, 100, fp);

	if (strcmp(theScratchString, "% Triangulation") != 0)
		return NULL;

	/*
	 *	Allocate the TriangulationData.
	 */
	theTriangulationData = (TriangulationData *) malloc(sizeof(TriangulationData));
	if (theTriangulationData == NULL)
		uFatalError("ReadNewFileFormat1", "unix file io");
	theTriangulationData->name				= NULL;
	theTriangulationData->cusp_data			= NULL;
	theTriangulationData->tetrahedron_data	= NULL;

	/*
	 *	Allocate and read the name of the manifold.
	 */
	theTriangulationData->name = (char *) malloc(1000 * sizeof(char));
	if (theTriangulationData->name == NULL)
		uFatalError("ReadNewFileFormat", "unix file io");
	/*
	 *	The name will be on the first nonempty line.
	 */

	my_fgets(theTriangulationData->name, 1000, fp);

	/*
	 *	Read the filled solution type.
	 */
	fscanf(fp, "%s", theScratchString);

	if (strcmp(theScratchString, "not_attempted") == 0)
		theTriangulationData->solution_type = not_attempted;
	else if (strcmp(theScratchString, "geometric_solution") == 0)
		theTriangulationData->solution_type = geometric_solution;
	else if (strcmp(theScratchString, "nongeometric_solution") == 0)
		theTriangulationData->solution_type = nongeometric_solution;
	else if (strcmp(theScratchString, "flat_solution") == 0)
		theTriangulationData->solution_type = flat_solution;
	else if (strcmp(theScratchString, "degenerate_solution") == 0)
		theTriangulationData->solution_type = degenerate_solution;
	else if (strcmp(theScratchString, "other_solution") == 0)
		theTriangulationData->solution_type = other_solution;
	else if (strcmp(theScratchString, "no_solution") == 0)
		theTriangulationData->solution_type = no_solution;
	else
		uFatalError("ReadNewFileFormat", "unix file io");

	/*
	 *	Read the volume.
	 */
	fscanf(fp, "%lf", &theTriangulationData->volume);

	/*
	 *	Read the orientability.
	 */
	fscanf(fp, "%s", theScratchString);
		 if (strcmp(theScratchString, "oriented_manifold") == 0)
		theTriangulationData->orientability = oriented_manifold;
	else if (strcmp(theScratchString, "nonorientable_manifold") == 0)
		theTriangulationData->orientability = nonorientable_manifold;
	else if (strcmp(theScratchString, "unknown_orientability") == 0)
		theTriangulationData->orientability = unknown_orientability;
	else
		uFatalError("ReadNewFileFormat", "unix file io");

	/*
	 *	Read the Chern-Simons invariant, if present.
	 */
	fscanf(fp, "%s", theScratchString);
		 if (strcmp(theScratchString, "CS_known") == 0)
		theTriangulationData->CS_value_is_known = TRUE;
	else if (strcmp(theScratchString, "CS_unknown") == 0)
		theTriangulationData->CS_value_is_known = FALSE;
	else
		uFatalError("ReadNewFileFormat", "unix file io");
	if (theTriangulationData->CS_value_is_known == TRUE)
		fscanf(fp, "%lf", &theTriangulationData->CS_value);
	else
		theTriangulationData->CS_value = 0.0;

	/*
	 *	Read the number of cusps, allocate an array for the cusp data,
	 *	and read the cusp data.
	 */
	fscanf(fp, "%d%d",
			&theTriangulationData->num_or_cusps,
			&theTriangulationData->num_nonor_cusps);
	theTotalNumCusps = theTriangulationData->num_or_cusps
					 + theTriangulationData->num_nonor_cusps;
	theTriangulationData->cusp_data = (CuspData *) malloc(theTotalNumCusps * sizeof(CuspData));
	if (theTriangulationData->cusp_data == NULL)
		uFatalError("ReadNewFileFormat", "unix file io");
	for (i = 0; i < theTotalNumCusps; i++)
	{
		if (fscanf(fp, "%s%lf%lf",
				theScratchString,
				&theTriangulationData->cusp_data[i].m,
				&theTriangulationData->cusp_data[i].l) != 3)
			uFatalError("ReadNewFileFormat", "unix file io");
		switch (theScratchString[0])
		{
			case 't':
			case 'T':
				theTriangulationData->cusp_data[i].topology = torus_cusp;
				break;

			case 'k':
			case 'K':
				theTriangulationData->cusp_data[i].topology = Klein_cusp;
				break;

			default:
				uFatalError("ReadNewFileFormat", "unix file io");
		}
	}

	/*
	 *	Read the number of tetrahedra, allocate an array for the
	 *	tetrahedron data, and read the tetrahedron data.
	 */
	fscanf(fp, "%d", &theTriangulationData->num_tetrahedra);
	theTriangulationData->tetrahedron_data = (TetrahedronData *) malloc(theTriangulationData->num_tetrahedra * sizeof(TetrahedronData));
	if (theTriangulationData->tetrahedron_data == NULL)
		uFatalError("ReadNewFileFormat", "unix file io");
	for (i = 0; i < theTriangulationData->num_tetrahedra; i++)
	{
		/*
		 *	Read the neighbor indices.
		 */
		for (j = 0; j < 4; j++)
		{
			fscanf(fp, "%d", &theTriangulationData->tetrahedron_data[i].neighbor_index[j]);
			if (theTriangulationData->tetrahedron_data[i].neighbor_index[j] < 0
			 || theTriangulationData->tetrahedron_data[i].neighbor_index[j] >= theTriangulationData->num_tetrahedra)
				uFatalError("ReadNewFileFormat", "unix file io");
		}

		/*
		 *	Read the gluings.
		 */
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
			{
				fscanf(fp, "%1d", &theTriangulationData->tetrahedron_data[i].gluing[j][k]);
				if (theTriangulationData->tetrahedron_data[i].gluing[j][k] < 0
				 || theTriangulationData->tetrahedron_data[i].gluing[j][k] > 3)
					uFatalError("ReadNewFileFormat", "unix file io");
			}

		/*
		 *	Read the cusp indices.
		 *
		 *	99/06/04  Allow an index of -1 on "cusps" that are
		 *	really finite vertices.
		 */
		for (j = 0; j < 4; j++)
		{
			fscanf(fp, "%d", &theTriangulationData->tetrahedron_data[i].cusp_index[j]);
			if (theTriangulationData->tetrahedron_data[i].cusp_index[j] < -1
			 || theTriangulationData->tetrahedron_data[i].cusp_index[j] >= theTotalNumCusps)
				uFatalError("ReadNewFileFormat", "unix file io");
		}

		/*
		 *	Read the peripheral curves.
		 */
		for (j = 0; j < 2; j++)			 /* meridian, longitude     */
			for (k = 0; k < 2; k++)		 /* righthanded, lefthanded */
				for (v = 0; v < 4; v++)
					for (f = 0; f < 4; f++)
						fscanf(fp, "%d", &theTriangulationData->tetrahedron_data[i].curve[j][k][v][f]);

		/*
		 *	Read the filled shape (which the kernel ignores).
		 */
		fscanf(fp, "%lf%lf",
			&theTriangulationData->tetrahedron_data[i].filled_shape.real,
			&theTriangulationData->tetrahedron_data[i].filled_shape.imag);
	}

	return theTriangulationData;
}


void save_triangulation(
	Triangulation	*manifold,
	char			*file_name)
{
	TriangulationData	*theTriangulationData;
	FILE				*fp;

	/*
	 *	If the file_name is nonempty, write the file.
	 *	If the file_name is empty, write to stdout.
	 */
	if (strlen(file_name) > 0)
	{
		fp = fopen(file_name, "w");
		if (fp == NULL)
		{
			printf("couldn't open %s\n", file_name);
			return;
		}
	}
	else
		fp = stdout;

	triangulation_to_data(manifold, &theTriangulationData);
	WriteNewFileFormat(fp, theTriangulationData);
	free_triangulation_data(theTriangulationData);

	if (fp != stdout)
		fclose(fp);
}


extern void WriteNewFileFormat(
	FILE				*fp,
	TriangulationData	*data)
{
	int	i,
		j,
		k,
		v,
		f;

	fprintf(fp, "%% Triangulation\n");

	if (data->name != NULL)
		fprintf(fp, "%s\n", data->name);
	else
		fprintf(fp, "untitled\n");

	switch (data->solution_type)
	{
		case not_attempted:
			fprintf(fp, "not_attempted");
			break;

		case geometric_solution:
			fprintf(fp, "geometric_solution");
			break;

		case nongeometric_solution:
			fprintf(fp, "nongeometric_solution");
			break;

		case flat_solution:
			fprintf(fp, "flat_solution");
			break;

		case degenerate_solution:
			fprintf(fp, "degenerate_solution");
			break;

		case other_solution:
			fprintf(fp, "other_solution");
			break;

		case no_solution:
			fprintf(fp, "no_solution");
			break;

		default:
			fprintf(fp, "not_attempted");
			break;
	}

	if (data->solution_type != not_attempted)
		fprintf(fp, "  %.8f\n", data->volume);
	else
		fprintf(fp, "  %.1f\n", 0.0);

	switch (data->orientability)
	{
		case oriented_manifold:
			fprintf(fp, "oriented_manifold\n");
			break;

		case nonorientable_manifold:
			fprintf(fp, "nonorientable_manifold\n");
			break;

		case oriented_orbifold: /* DJH */
			fprintf(fp, "oriented_orbifold\n");
			break;

		case nonorientable_orbifold: /* DJH */
			fprintf( fp, "nonorientable_orbifold\n");
			break;
	}

	if (data->CS_value_is_known == TRUE)
		fprintf(fp, "CS_known %.16f\n", data->CS_value);
	else
		fprintf(fp, "CS_unknown\n");

	fprintf(fp, "\n%d %d\n", data->num_or_cusps, data->num_nonor_cusps);
	for (i = 0; i < data->num_or_cusps + data->num_nonor_cusps; i++)
		fprintf(fp, "    %s %16.12f %16.12f\n",
			(data->cusp_data[i].topology == torus_cusp) ? "torus" : "Klein",
			data->cusp_data[i].m,
			data->cusp_data[i].l);
	fprintf(fp, "\n");

	fprintf(fp, "%d\n", data->num_tetrahedra);
	for (i = 0; i < data->num_tetrahedra; i++)
	{
		for (j = 0; j < 4; j++)
			fprintf(fp, "%4d ", data->tetrahedron_data[i].neighbor_index[j]);
		fprintf(fp, "\n");

		for (j = 0; j < 4; j++)
		{
			fprintf(fp, " ");
			for (k = 0; k < 4; k++)
				fprintf(fp, "%d", data->tetrahedron_data[i].gluing[j][k]);
		}
		fprintf(fp, "\n");

		for (j = 0; j < 4; j++)
			fprintf(fp, "%4d ", data->tetrahedron_data[i].cusp_index[j]);
		fprintf(fp, "\n");

		for (j = 0; j < 2; j++)			/* meridian, longitude     */
			for (k = 0; k < 2; k++)		/* righthanded, lefthanded */
			{
				for (v = 0; v < 4; v++)
					for (f = 0; f < 4; f++)
						fprintf(fp, " %2d", data->tetrahedron_data[i].curve[j][k][v][f]);
				fprintf(fp, "\n");
			}

		if (data->solution_type != not_attempted)
			fprintf(fp, "%16.12f %16.12f\n\n",
				data->tetrahedron_data[i].filled_shape.real,
				data->tetrahedron_data[i].filled_shape.imag);
		else
			fprintf(fp, "%3.1f %3.1f\n\n", 0.0, 0.0);
	}
}

void my_fgets( char *string, int max, FILE *fp )
{
	int c, ch;

        c = 0; 

	/* remove newlines at the beginning */

        while( (ch = getc( fp )) == '\n' || ch == '\r' )
                rewind( fp );
        ungetc(ch, fp);
        while( (ch = getc( fp )) != '\n' && ch != '\r' && c < max-1 )
                string[c++] = ch;
        string[c] = '\0';
}

