/*
 *  unix_file_io.c
 *
 *  This hacked together file allows unix-style programs
 *  to read and save Triangulations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unix_file_io.h"

#include "kernel_namespace.h"

#define READ_OLD_FILE_FORMAT 0

#define FALSE   0
#define TRUE    1

/*
 *  gcc complains about the
 *
 *      use of `l' length character with `f' type character
 *
 *  in fprintf() calls.  Presumably it considers the 'l' unnecessary
 *  because even floats would undergo default promotion to doubles
 *  in the function call (see section A7.3.2 in Appendix A of K&R 2nd ed.).
 *  Therefore I've changed "%lf" to "%f" in all fprintf() calls.
 *  If this makes trouble on your system, change it back, and please
 *  let me know (weeks@northnet.org).
 */

/* Modified 2010/09/20 by Marc Culler to avoid quadratic time behavior
 * caused by the fact that sscanf calls strlen.
 */

/* Size for our temporary string buffers.
 */
#define SBSIZE 256

/* Copy up to the specified number of words of the string str into the
 * specified buffer, which must hold at least SBSIZE chars.  Return
 * the number of chars copied.  This assumes that any line in a
 * triangulation file will have length < SBSIZE.
 */ 
#define whitespace(x) (x=='\040'||x=='\f'||x=='\n'||x=='\r'||x=='\t'||x=='\v')
 
static int read_head(char *headbuf, const char *str, int words) {
    char *ptr = headbuf;
    int n;
    for (n=0; n < words; n++) {
        /* copy whitespace */
        while ( whitespace(*str) && ptr-headbuf < SBSIZE-1 )
            *ptr++ = *str++;
	/* copy non-whitespace */
	while ( *str != '\0' && !whitespace(*str) && ptr-headbuf < SBSIZE-1)
             *ptr++ = *str++;
    if ( *str == '\0' ) break;
  }
  /* add trailing null char */
  *ptr = '\0';
  return (ptr - headbuf); 
}

/* Modified 04/23/09 by Marc Culler to allow reading a triangulation
 * from a string.
 */

static TriangulationData    *ReadNewFileFormat(const char *buffer);
static void                 WriteNewFileFormat(FILE *fp, TriangulationData *data);
static char                 *StringNewFileFormat(TriangulationData *data);

#if READ_OLD_FILE_FORMAT
extern FuncResult           read_old_manifold(FILE *fp, Triangulation **manifold);
#endif

/* Modified 2010/7/27 by NMD to allow for Classic Mac OS and Windows line endings, as well Unix ones */

static Boolean  is_eol_char(const char *buffer);

Boolean is_eol_char(const char *buffer){
  return (*buffer == '\n') ||  (*buffer == '\r');
}

/* Modified 04/24/09 by Marc Culler to allow extracting a triangulation from a C string. */

Triangulation *read_triangulation_from_string(
    const char    *file_data)
{
    TriangulationData   *theTriangulationData;
    Triangulation       *manifold;

    if ( strncmp("% Triangulation", file_data, 15) != 0)
      uFatalError("read_triangulation_from_string", "unix file io");
	  
    theTriangulationData = ReadNewFileFormat(file_data);
    
    data_to_triangulation(theTriangulationData, &manifold);
    
    free(theTriangulationData->name);
    free(theTriangulationData->cusp_data);
    free(theTriangulationData->tetrahedron_data);
    free(theTriangulationData);

    return manifold;
}

Triangulation *read_triangulation(
    const char    *file_name)
{
    FILE            *fp;
    Boolean         theNewFormat = FALSE;
    Triangulation   *manifold = NULL;

    /*
     *  If the file_name is nonempty, read the file.
     *  If the file_name is empty, read from stdin.
     */
    if (strlen(file_name) > 0)
    {
        fp = fopen(file_name, "rb");
        if (fp == NULL)
            return NULL;

        /*
         *  Take a peek at the first line to see whether this is
         *  the new file format or the old one.
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
        TriangulationData   *theTriangulationData;
	long filesize = 0;
	char *buffer;
	  
	if ( fseek(fp, 0, SEEK_END) != 0    ||
	     ( filesize = ftell(fp) ) == -1 ||
	     fseek(fp, 0, SEEK_SET) != 0     )
	    uFatalError("read_triangulation", "unix file io");

	buffer = (char *)malloc(filesize + 1);
	if ( buffer == NULL)
	    uFatalError("read_triangulation", "unix file io");
	if ( fread(buffer, filesize, 1, fp) != 1 )
	    uFatalError("read_triangulation", "unix file io");

	buffer[filesize] = '\0';
	theTriangulationData = ReadNewFileFormat(buffer);
	free(buffer);
    
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
        uFatalError("read_triangulation", "unix file io");
#endif
    }

    if (fp != stdin)
        fclose(fp);

    return manifold;
}


static TriangulationData *ReadNewFileFormat(
    const char *buffer)
{
    char                *ptr;
    char                theScratchString[SBSIZE];
    char                HeadBuffer[SBSIZE];
    int                 count;
    TriangulationData   *theTriangulationData;
    int                 theTotalNumCusps,
                        i,
                        j,
                        k,
                        v,
                        f;
    size_t              size = 0;
    double              temp, temp_m, temp_l, temp_r, temp_i;


    /*
     *  Read and ignore the header (% Triangulation).
     */
    while (!is_eol_char(buffer)) buffer++;

    /*
     *  Allocate the TriangulationData.
     */
    theTriangulationData = (TriangulationData *) malloc(sizeof(TriangulationData));
    if (theTriangulationData == NULL)
        uFatalError("ReadNewFileFormat 1", "unix file io");
    theTriangulationData->name              = NULL;
    theTriangulationData->cusp_data         = NULL;
    theTriangulationData->tetrahedron_data  = NULL;

    /*
     *  Allocate and read the name of the manifold.
     *  The name will be on the first nonempty line.
     */

   while (is_eol_char(buffer)) buffer++;
   i = 0;
   ptr = NULL;
   while (!is_eol_char(buffer)){
     if (i + 2 > (int) size){
	   size += 100;
	   ptr = (char *) realloc(ptr, size);
	   if (ptr == NULL){
	       /* We could have saved the old pointer, to free it here.
		* But why bother?  We are toast anyway.
		*/
	       uFatalError("ReadNewFileFormat 2", "unix file io");
	       return NULL;
	   }
       }
       ptr[i++] = *buffer++;
   }
   ptr[i] = 0;
   theTriangulationData->name = ptr;
    /*
     *  Read the filled solution type.
     */
    buffer += read_head(HeadBuffer, buffer, 1);
    sscanf(HeadBuffer, "%s", theScratchString);
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
    else if (strcmp(theScratchString, "externally_computed") == 0)
        theTriangulationData->solution_type = externally_computed;
    else
        uFatalError("ReadNewFileFormat 3", "unix file io");

    /*
     *  Read the volume.
     */
    buffer += read_head(HeadBuffer, buffer, 1);
    sscanf(HeadBuffer, "%lf", &temp);
    theTriangulationData->volume = temp;


    /*
     *  Read the orientability.
     */
    buffer += read_head(HeadBuffer, buffer, 1); 
    sscanf(HeadBuffer, "%s", theScratchString);
    if (strcmp(theScratchString, "oriented_manifold") == 0)
        theTriangulationData->orientability = oriented_manifold;
    else if (strcmp(theScratchString, "nonorientable_manifold") == 0)
        theTriangulationData->orientability = nonorientable_manifold;
    else if (strcmp(theScratchString, "unknown_orientability") == 0)
        theTriangulationData->orientability = unknown_orientability;
    else
        uFatalError("ReadNewFileFormat 4", "unix file io");

    /*
     *  Read the Chern-Simons invariant, if present.
     */
    buffer += read_head(HeadBuffer, buffer, 1); 
    sscanf(HeadBuffer, "%s", theScratchString);
    if (strcmp(theScratchString, "CS_known") == 0)
        theTriangulationData->CS_value_is_known = TRUE;
    else if (strcmp(theScratchString, "CS_unknown") == 0)
        theTriangulationData->CS_value_is_known = FALSE;
    else
        uFatalError("ReadNewFileFormat 5", "unix file io");
    if (theTriangulationData->CS_value_is_known == TRUE) {
      buffer += read_head(HeadBuffer, buffer, 1);
      sscanf(HeadBuffer, "%lf", &temp);
      theTriangulationData->CS_value = temp;

      }
    else
        theTriangulationData->CS_value = 0.0;

    /*
     *  Read the number of cusps, allocate an array for the cusp data,
     *  and read the cusp data.
     */
    buffer += read_head(HeadBuffer, buffer, 2); 
    sscanf(HeadBuffer, "%d%d",
	   &theTriangulationData->num_or_cusps,
           &theTriangulationData->num_nonor_cusps);
    theTotalNumCusps = theTriangulationData->num_or_cusps
                     + theTriangulationData->num_nonor_cusps;
    theTriangulationData->cusp_data = (CuspData *) malloc(theTotalNumCusps * sizeof(CuspData));
    if (theTriangulationData->cusp_data == NULL)
        uFatalError("ReadNewFileFormat 6", "unix file io");
    for (i = 0; i < theTotalNumCusps; i++)
    {
        buffer += read_head(HeadBuffer, buffer, 3);
        if (sscanf(HeadBuffer, "%s%lf%lf",
		   theScratchString,
		   &temp_m,
		   &temp_l) != 3)
            uFatalError("ReadNewFileFormat 7", "unix file io");
	theTriangulationData->cusp_data[i].m = temp_m;
	theTriangulationData->cusp_data[i].l = temp_l;
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
                uFatalError("ReadNewFileFormat 8", "unix file io");
        }
    }

    /*
     *  Read the number of tetrahedra, allocate an array for the
     *  tetrahedron data, and read the tetrahedron data.
     */
    buffer += read_head(HeadBuffer, buffer, 1); 
    sscanf(HeadBuffer, "%d", &theTriangulationData->num_tetrahedra);
    theTriangulationData->tetrahedron_data = (TetrahedronData *) malloc(theTriangulationData->num_tetrahedra * sizeof(TetrahedronData));
    if (theTriangulationData->tetrahedron_data == NULL)
        uFatalError("ReadNewFileFormat 9", "unix file io");
    for (i = 0; i < theTriangulationData->num_tetrahedra; i++)
    {
        /*
         *  Read the neighbor indices.
         */
        for (j = 0; j < 4; j++)
        {
          buffer += read_head(HeadBuffer, buffer, 1); 
	  sscanf(HeadBuffer, "%d", &theTriangulationData->tetrahedron_data[i].neighbor_index[j]);
	  if (theTriangulationData->tetrahedron_data[i].neighbor_index[j] < 0
             || theTriangulationData->tetrahedron_data[i].neighbor_index[j] >= theTriangulationData->num_tetrahedra)
                uFatalError("ReadNewFileFormat 10", "unix file io");
        }

        /*
         *  Read the gluings.
         */
	/* This assumes that the gluings are groups of 4 digits with no
           whitespace between, but with whitespace between groups. */
	buffer += read_head(HeadBuffer, buffer, 4);
	ptr = HeadBuffer;
        for (j = 0; j < 4; j++)
            for (k = 0; k < 4; k++)
            {
	      sscanf(ptr, "%1d%n", &theTriangulationData->tetrahedron_data[i].gluing[j][k], &count);
	      ptr += count;
	      if (theTriangulationData->tetrahedron_data[i].gluing[j][k] < 0
                 || theTriangulationData->tetrahedron_data[i].gluing[j][k] > 3)
                    uFatalError("ReadNewFileFormat 11", "unix file io");
            }

        /*
         *  Read the cusp indices.
         *
         *  99/06/04  Allow an index of -1 on "cusps" that are
         *  really finite vertices.
         */
        for (j = 0; j < 4; j++)
        {
	  buffer += read_head(HeadBuffer, buffer, 1); 
	  sscanf(HeadBuffer, "%d", &theTriangulationData->tetrahedron_data[i].cusp_index[j]);
	  if (theTriangulationData->tetrahedron_data[i].cusp_index[j] < -1
             || theTriangulationData->tetrahedron_data[i].cusp_index[j] >= theTotalNumCusps)
                uFatalError("ReadNewFileFormat 12", "unix file io");
        }

        /*
         *  Read the peripheral curves.
         */
        for (j = 0; j < 2; j++)          /* meridian, longitude     */
            for (k = 0; k < 2; k++)      /* righthanded, lefthanded */
                for (v = 0; v < 4; v++)
		  for (f = 0; f < 4; f++){
		      buffer += read_head(HeadBuffer, buffer, 1);
		      sscanf(HeadBuffer, "%d", &theTriangulationData->tetrahedron_data[i].curve[j][k][v][f]);
		  }

        /*
         *  Read the filled shape (which the kernel ignores).
         */
	buffer += read_head(HeadBuffer, buffer, 2);
        sscanf(HeadBuffer, "%lf%lf",
	       &temp_r,
	       &temp_i);
	theTriangulationData->tetrahedron_data[i].filled_shape.real = temp_r;
	theTriangulationData->tetrahedron_data[i].filled_shape.imag = temp_i;
    }

    return theTriangulationData;
}


Boolean write_triangulation(
    Triangulation   *manifold,
    const char            *file_name)
{
    TriangulationData   *theTriangulationData;
    FILE                *fp;

    /*
     *  If the file_name is nonempty, write the file.
     *  If the file_name is empty, write to stdout.
     */
    if (strlen(file_name) > 0)
    {
        fp = fopen(file_name, "w");
        if (fp == NULL)
        {
            printf("couldn't open %s\n", file_name);
            return FALSE;
        }
    }
    else
        fp = stdout;

    triangulation_to_data(manifold, &theTriangulationData);
    WriteNewFileFormat(fp, theTriangulationData);
    free_triangulation_data(theTriangulationData);

    if (fp != stdout)
        fclose(fp);
    return TRUE;
}


static void WriteNewFileFormat(
    FILE                *fp,
    TriangulationData   *data)
{
    int i,
        j,
        k,
        v,
        f;

    fprintf(fp, "%% Triangulation\n");

    if (data->name != NULL)
        fprintf(fp, "%s\n", data->name);
    else
        fprintf(fp, "untitled");

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

        case externally_computed:
            fprintf(fp, "externally_computed");
            break;
    }

    if (data->solution_type != not_attempted && data->solution_type != externally_computed)
        fprintf(fp, "  %.8f\n", (double)data->volume);
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
            
        case unknown_orientability:
            /* This manifold is garbage */
            fprintf(fp, "ERROR: orientability not computed!\n");
            break;
    }

    if (data->CS_value_is_known == TRUE)
        fprintf(fp, "CS_known %.16f\n", (double)data->CS_value);
    else
        fprintf(fp, "CS_unknown\n");

    fprintf(fp, "\n%d %d\n", data->num_or_cusps, data->num_nonor_cusps);
    for (i = 0; i < data->num_or_cusps + data->num_nonor_cusps; i++)
        fprintf(fp, "    %s %16.12f %16.12f\n",
            (data->cusp_data[i].topology == torus_cusp) ? "torus" : "Klein",
		(double)data->cusp_data[i].m,
		(double)data->cusp_data[i].l);
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

        for (j = 0; j < 2; j++)         /* meridian, longitude     */
            for (k = 0; k < 2; k++)     /* righthanded, lefthanded */
            {
                for (v = 0; v < 4; v++)
                    for (f = 0; f < 4; f++)
                        fprintf(fp, " %2d", data->tetrahedron_data[i].curve[j][k][v][f]);
                fprintf(fp, "\n");
            }

        if (data->solution_type != not_attempted)
            fprintf(fp, "%16.12f %16.12f\n\n",
		    (double)data->tetrahedron_data[i].filled_shape.real,
		    (double)data->tetrahedron_data[i].filled_shape.imag);
        else
            fprintf(fp, "%3.1f %3.1f\n\n", 0.0, 0.0);
    }
}

/* Added by Marc Culler 2010-12-17 to allow writing a triangulation
 * to a string.  Memory is malloc'ed for the string.  Caller must free.
 */

char *string_triangulation(
    Triangulation   *manifold)
{
    TriangulationData   *theTriangulationData;
    char                *result;

    triangulation_to_data(manifold, &theTriangulationData);
    result = StringNewFileFormat(theTriangulationData);
    free_triangulation_data(theTriangulationData);
    return result;
}

static char *StringNewFileFormat(
    TriangulationData   *data)
{
    int i,
        j,
        k,
        v,
        f,
        size;
    char *buffer;
    char *p;

    size = 100*(10 + data->num_or_cusps + data->num_nonor_cusps + 8*data->num_tetrahedra);
    buffer = (char *)malloc(size);
    if ( buffer == NULL)
      uFatalError("StringNewFileFormat", "unix file io");
    p = buffer;

    p += sprintf(p, "%% Triangulation\n");

    if (data->name != NULL)
        p += sprintf(p, "%s\n", data->name);
    else
        p += sprintf(p, "untitled\n");

    switch (data->solution_type)
    {
        case not_attempted:
            p += sprintf(p, "not_attempted");
            break;

        case geometric_solution:
            p += sprintf(p, "geometric_solution");
            break;

        case nongeometric_solution:
            p += sprintf(p, "nongeometric_solution");
            break;

        case flat_solution:
            p += sprintf(p, "flat_solution");
            break;

        case degenerate_solution:
            p += sprintf(p, "degenerate_solution");
            break;

        case other_solution:
            p += sprintf(p, "other_solution");
            break;

        case no_solution:
            p += sprintf(p, "no_solution");
            break;

        case externally_computed:
            p += sprintf(p, "externally_computed");
            break;

    }

    if (data->solution_type != not_attempted && data->solution_type != externally_computed)
        p += sprintf(p, "  %.8f\n", (double)data->volume);
    else
        p += sprintf(p, "  %.1f\n", 0.0);

    switch (data->orientability)
    {
        case oriented_manifold:
            p += sprintf(p, "oriented_manifold\n");
            break;

        case nonorientable_manifold:
            p += sprintf(p, "nonorientable_manifold\n");
            break;
        case unknown_orientability:
            /* This manifold is garbage */
            p += sprintf(p, "ERROR: orientability not computed!\n");
            break;
    }

    if (data->CS_value_is_known == TRUE)
        p += sprintf(p, "CS_known %.16f\n", (double)data->CS_value);
    else
        p += sprintf(p, "CS_unknown\n");

    p += sprintf(p, "\n%d %d\n", data->num_or_cusps, data->num_nonor_cusps);
    for (i = 0; i < data->num_or_cusps + data->num_nonor_cusps; i++)
        p += sprintf(p, "    %s %16.12f %16.12f\n",
            (data->cusp_data[i].topology == torus_cusp) ? "torus" : "Klein",
		     (double)data->cusp_data[i].m,
		     (double)data->cusp_data[i].l);
    p += sprintf(p, "\n");

    p += sprintf(p, "%d\n", data->num_tetrahedra);
    for (i = 0; i < data->num_tetrahedra; i++)
    {
        for (j = 0; j < 4; j++)
            p += sprintf(p, "%4d ", data->tetrahedron_data[i].neighbor_index[j]);
        p += sprintf(p, "\n");

        for (j = 0; j < 4; j++)
        {
            p += sprintf(p, " ");
            for (k = 0; k < 4; k++)
                p += sprintf(p, "%d", data->tetrahedron_data[i].gluing[j][k]);
        }
        p += sprintf(p, "\n");

        for (j = 0; j < 4; j++)
            p += sprintf(p, "%4d ", data->tetrahedron_data[i].cusp_index[j]);
        p += sprintf(p, "\n");

        for (j = 0; j < 2; j++)         /* meridian, longitude     */
            for (k = 0; k < 2; k++)     /* righthanded, lefthanded */
            {
                for (v = 0; v < 4; v++)
                    for (f = 0; f < 4; f++)
                        p += sprintf(p, " %2d", data->tetrahedron_data[i].curve[j][k][v][f]);
                p += sprintf(p, "\n");
            }

        if (data->solution_type != not_attempted && data->solution_type != externally_computed)
            p += sprintf(p, "%16.12f %16.12f\n\n",
		 (double)data->tetrahedron_data[i].filled_shape.real,
		 (double)data->tetrahedron_data[i].filled_shape.imag);
        else
            p += sprintf(p, "%3.1f %3.1f\n\n", 0.0, 0.0);
   }
   return buffer;
}

#include "end_namespace.h"
