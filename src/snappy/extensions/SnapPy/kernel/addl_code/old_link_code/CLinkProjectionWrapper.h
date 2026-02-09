/*
 *	CLinkProjectionWrapper.h
 *
 *	This isn't really a "wrapper" in the sense that CTriangulationWrapper
 *	is a wrapper, but it seemed like a good idea for all my
 *	itsContents files to have similar names.
 */

#include <stdio.h>

class CLinkProjectionWrapper
{
public:

	//	Let its data members be public, so CLinkProjection::ContentsToWindow()
	//	and CLinkProjection::WindowToContents() can access them
	//	without a lot of fussing around.

	int	itsNumComponents,
		*itsFirstVertices,
		*itsLastVertices;

	int	itsNumVertices,
		*itsHCoordinates,	//	left-to-right
		*itsVCoordinates;	//	top-to-bottom

	int	itsNumEdges,
		*itsBackwardVertices,
		*itsForwardVertices;

	int	itsNumCrossings,
		*itsUnderstrands,
		*itsOverstrands;

	int	itsHotVertex;	// = -1 if none

					CLinkProjectionWrapper(void);
	virtual			~CLinkProjectionWrapper(void);

	virtual void	ReadFile(FILE *fp);
	virtual void	WriteFile(FILE *fp);

	virtual void	ClearContents(void);

protected:

	virtual void	ReadNewFileFormat(FILE *fp);
	virtual void	ReadOldFileFormat(FILE *fp);
};
