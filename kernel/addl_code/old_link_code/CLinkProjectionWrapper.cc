/*
 *	CLinkProjectionWrapper.cp
 *
 *	This isn't really a "wrapper" in the sense that CTriangulationWrapper
 *	is a wrapper, but it seemed like a good idea for all my
 *	itsContents files to have similar names.
 */

extern "C"{
  #include "kernel.h"
}

#include "CLinkProjectionWrapper.h"

CLinkProjectionWrapper::CLinkProjectionWrapper(void)
{
	itsNumComponents	= 0;
	itsFirstVertices	= NULL;
	itsLastVertices		= NULL;

	itsNumVertices		= 0;
	itsHCoordinates		= NULL;
	itsVCoordinates		= NULL;

	itsNumEdges			= 0;
	itsBackwardVertices	= NULL;
	itsForwardVertices	= NULL;

	itsNumCrossings		= 0;
	itsUnderstrands		= NULL;
	itsOverstrands		= NULL;

	itsHotVertex		 = -1;

}


CLinkProjectionWrapper::~CLinkProjectionWrapper(void)
{

	ClearContents();
}


void CLinkProjectionWrapper::ClearContents(void)
{
	itsNumComponents = 0;

	itsNumVertices = 0;

	itsNumEdges = 0;

	itsNumCrossings = 0;

	itsHotVertex = -1;
}


void CLinkProjectionWrapper::ReadFile(
	FILE	*fp)
{
	char	theFirstChar;

	//	Typically we expect the current file format,
	//	but we want to read (but not write) the old format as well.

	fscanf(fp, "%c", &theFirstChar);
	rewind(fp);
	if (theFirstChar == '%')
		ReadNewFileFormat(fp);
	else
		ReadOldFileFormat(fp);
}


void CLinkProjectionWrapper::ReadNewFileFormat(
	FILE	*fp)
{
	int		i;
	char	theIgnoredString[100];

	//	ReadFile() should be called only for a "fresh" itsContents, so
	//	in principle the following ClearContents() call is unnecessary.
	ClearContents();

	//	Skip the header "% Link Projection".
	fgets(theIgnoredString, 100, fp);

	fscanf(fp, "%d", &itsNumComponents);
	itsFirstVertices = (int *) my_malloc(itsNumComponents * sizeof(int));
	itsLastVertices  = (int *) my_malloc(itsNumComponents * sizeof(int));
	for (i = 0; i < itsNumComponents; i++)
		fscanf(fp, "%d%d", &itsFirstVertices[i],  &itsLastVertices[i]);

	fscanf(fp, "%d", &itsNumVertices);
	itsHCoordinates = (int *) my_malloc(itsNumVertices * sizeof(int));
	itsVCoordinates = (int *) my_malloc(itsNumVertices * sizeof(int));
	for (i = 0; i < itsNumVertices; i++)
		fscanf(fp, "%d%d", &itsHCoordinates[i],  &itsVCoordinates[i]);

	fscanf(fp, "%d", &itsNumEdges);
	itsBackwardVertices = (int *) my_malloc(itsNumEdges * sizeof(int));
	itsForwardVertices  = (int *) my_malloc(itsNumEdges * sizeof(int));
	for (i = 0; i < itsNumEdges; i++)
		fscanf(fp, "%d%d", &itsBackwardVertices[i],  &itsForwardVertices[i]);

	fscanf(fp, "%d", &itsNumCrossings);
	itsUnderstrands = (int *) my_malloc(itsNumCrossings * sizeof(int));
	itsOverstrands  = (int *) my_malloc(itsNumCrossings * sizeof(int));
	for (i = 0; i < itsNumCrossings; i++)
		fscanf(fp, "%d%d", &itsUnderstrands[i],  &itsOverstrands[i]);

	fscanf(fp, "%d", &itsHotVertex);
}


void CLinkProjectionWrapper::ReadOldFileFormat(
	FILE	*fp)
{
	char	theBuffer[100];
	int		i,
			theIndex,
			theNextIndex,
			theHPosition,
			theVPosition;

	//	ReadFile() should be called only for a "fresh" itsContents, so
	//	in principle the following ClearContents() call is unnecessary.
	ClearContents();

	//	snappea 1.3.x used the num_free_loops variable, while at one point
	//	"snap pea" did not.  Most likely all files we encounter will include
	//	it (it will be zero), but let's take the robust approach of allowing
	//	for either possibility.  (theBuffer may contain four ints or it may
	//	contain only three -- we don't care which.)
	fgets(theBuffer, 100, fp);
	sscanf(theBuffer, "%d%d%d", &itsNumComponents, &itsNumEdges, &itsNumCrossings);

	//	The old format supported only closed (circular) link components,
	//	so the number of vertices equals the number of edges.
	itsNumVertices = itsNumEdges;

	//	Allocate all our arrays.

	itsFirstVertices = (int *) my_malloc(itsNumComponents * sizeof(int));
	itsLastVertices  = (int *) my_malloc(itsNumComponents * sizeof(int));

	itsHCoordinates = (int *) my_malloc(itsNumVertices * sizeof(int));
	itsVCoordinates = (int *) my_malloc(itsNumVertices * sizeof(int));

	itsBackwardVertices = (int *) my_malloc(itsNumEdges * sizeof(int));
	itsForwardVertices  = (int *) my_malloc(itsNumEdges * sizeof(int));

	itsUnderstrands = (int *) my_malloc(itsNumCrossings * sizeof(int));
	itsOverstrands  = (int *) my_malloc(itsNumCrossings * sizeof(int));

	//	The old format stored only edges, not vertices.
	//	We make the convention that each vertex inherit the
	//	index of its forwardEdge.

	//	For each component, read itsFirstVertex == itsLastVertex.
	for (i = 0; i < itsNumComponents; i++)
	{
		fscanf(fp, "%d", &itsFirstVertices[i]);
		itsLastVertices[i] = itsFirstVertices[i];
	}

	//	Read the edge and vertex information.
	//	(Use "assignment supression" to ignore the fourth number in each row.)
	for (i = 0; i < itsNumEdges; i++)
	{
		fscanf(fp, "%d%d%d%*d%d", &theIndex, &theHPosition, &theVPosition, &theNextIndex);
		itsBackwardVertices[theIndex] = theIndex;
		itsForwardVertices [theIndex] = theNextIndex;
		itsHCoordinates[theIndex] = theHPosition;
		itsVCoordinates[theIndex] = theVPosition;
	}

	//	Read the understrand and overstrand at each crossing,
	//	and ignore all the other information in each row.
	for (i = 0; i < itsNumCrossings; i++)
		fscanf(fp, "%*d%*d%*d%*d%d%d%*d%*d%*d%*d%*d",  &itsOverstrands[i], &itsUnderstrands[i]);

	//	The old format supported only closed link components,
	//	so there can be no hot vertex.
	itsHotVertex = -1;
}


void CLinkProjectionWrapper::WriteFile(
	FILE	*fp)
{
	int	i;

	fprintf(fp, "%% Link Projection\n");

	fprintf(fp, "\n%d\n", itsNumComponents);
	for (i = 0; i < itsNumComponents; i++)
		fprintf(fp, " %3d %3d\n", itsFirstVertices[i],  itsLastVertices[i]);

	fprintf(fp, "\n%d\n", itsNumVertices);
	for (i = 0; i < itsNumVertices; i++)
		fprintf(fp, " %3d %3d\n", itsHCoordinates[i],  itsVCoordinates[i]);

	fprintf(fp, "\n%d\n", itsNumEdges);
	for (i = 0; i < itsNumEdges; i++)
		fprintf(fp, " %3d %3d\n", itsBackwardVertices[i],  itsForwardVertices[i]);

	fprintf(fp, "\n%d\n", itsNumCrossings);
	for (i = 0; i < itsNumCrossings; i++)
		fprintf(fp, " %3d %3d\n", itsUnderstrands[i],  itsOverstrands[i]);

	fprintf(fp, "\n%d\n", itsHotVertex);
}
