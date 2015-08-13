/******************************************************************************
 CLinkProjectionPanorama.cp

				CLinkProjectionPanorama Panorama/EditText Class
	
	Copyright © 1995 Jeff Weeks, Geometry Center. All rights reserved.

	Modified 7/12/98 by Nathan Dunfield

 ******************************************************************************/

#include "CLinkProjectionPanorama.h"

extern "C"{
   #include "kernel.h"
}

#define LINE_WIDTH			2
#define HIT_EPSILON			4
#define CROSSING_EPSILON	4

//	We allow piecewise entry of links.  That is, the user can create
//	a few strands here and a few strands there, and eventually connect
//	them up.  Each strand is directed.  When the user connects two
//	strands with opposing directions, one direction or the other
//	must change.

//	At the highest level we keep track of LPComponents only ("LP" stands
//	for "Link Projection").  Each LPComponent keeps track of its
//	constituent LPVertices and LPEdges.  In the final link projection
//	all components will of course be circular.  But at the intermediate
//	stages some may be linear.  An LPEdge always has (non-NULL) pointers
//	to the LPVertices which come before and after it.  An LPVertex has
//	pointers to the LPEdges which come before and after it, except for
//	the first and last LPVertices in a linear component, whose backwardEdge
//	and forwardEdge fields respectively are NULL.

struct LPComponent
{
	//	The index determines the color.
	int			index;

	//	For a linear component, firstVertex and lastVertex are what
	//	you would expect.  For a circular component, firstVertex ==
	//	lastVertex may be an vertex in the component.
	LPVertex	*firstVertex,
				*lastVertex;

	//	The various components are kept on a NULL-terminated
	//	singly linked list.
	LPComponent	*next;
};

struct LPVertex
{
	//	Where is this vertex?
	Point		position;

	//	Which component does it belong to?
	LPComponent	*component;

	//	What are the incident edges?  "Backward" and "forward" are
	//	relative to the direction of the component.  A value of NULL
	//	means we're at the beginning or end of a linear component, and
	//	the corresponding edge is not yet present.
	LPEdge		*backwardEdge,
				*forwardEdge;

	//	The index is *not* maintained globally, but is used only
	//	temporarily in WindowToContents().
	int			index;
};

struct LPEdge
{
	//	What are the incident vertices?  "Backward" and "forward" are
	//	relative to the direction of the component.  These will never
	//	be NULL.
	LPVertex	*backwardVertex,
				*forwardVertex;

	//	Let this edge's endpoints be pt1 and pt2.  The equation of the line
	//	containing them is
	//
	//					y - pt1.v     pt2.v - pt1.v
	//					---------  =  -------------
	//					x - pt1.h     pt2.h - pt1.h
	//
	//	or
	//
	//	x*(pt2.v - pt1.v) + y*(pt1.h - pt2.h) = pt1.h*pt2.v - pt1.v*pt2.h
	//
	//	Write this last equation as
	//
	//						a*x + b*y = c
	//	where
	//					a = pt2.v - pt1.v
	//					b = pt1.h - pt2.h
	//					c = pt1.h * pt2.v  -  pt1.v * pt2.h
	//
	//	[There's also a nice dot product interpretation of this equation.]
	//
	//	The coordinates of pt1 and pt2 will typically be in the range
	//	0 to 1000, so we can check this equation using exact (32-bit)
	//	integer arithmetic.
	int			a,
				b,
				c;

	//	The index is *not* maintained globally, but is used only
	//	temporarily in WindowToContents().
	int			index;
};

struct LPCrossing
{
	//	Where is this crossing, rounded to integer coordinates?
                   Point		position;

	//	Keep track of a neighborhood of the crossing for
	//	(1) drawing it nicely, and
	//	(2) enforcing general position.
                   Rect		neighborhood;

	//	Which edges cross here?
	LPEdge		*understrand,
				*overstrand;

	//	The KLPCrossing is *not* maintained globally, but is set and used
	//	only locally within the code which converts to the KLPProjection
	//	format.
	KLPCrossing	*itsKLPCrossing;

	//	The index is *not* maintained globally, but is used only
	//	temporarily in WindowToContents().
	int			index;

	//	The visited_ fields are *not* maintained globally, but are used
	//	only temporarily in DoAlternate().
	Boolean		visited_once,
				visited_twice;

	//	The LPCrossings are kept on a NULL-terminated, singly linked list.
	LPCrossing	*next;
};

void SetRect(Rect* r, int a, int b, int c, int d){
  r->left = a;
  r->top = b;
  r->right = a + c;
  r->bottom = a + d;
}


CLinkProjectionPanorama::CLinkProjectionPanorama(void){
  itsComponentList = NULL;
  itsNumComponents = 0;
  itsHotVertex = NULL;
  itsCrossingList = NULL;
}

CLinkProjectionPanorama::~CLinkProjectionPanorama(void){
  FreeLinkProjection();
}


KLPProjection *CLinkProjectionPanorama::CreateKLPProjection(void)
{
	KLPProjection	*theKLPProjection;

	//	The link projection should be nonempty.
	if (itsComponentList == NULL)
	{
		uAcknowledge("Can't triangulate the complement of the empty link.  Click the Help button for instructions on drawing a link.");
		return NULL;
	}

	//	If some components aren't yet complete, return NULL.
	if (AllComponentsComplete() == FALSE)
	{
		uAcknowledge("All link components must be complete loops before you can triangulate the complement.");
		return NULL;
	}

	//	Allocate the basic KLPProjection.
	theKLPProjection = (KLPProjection *) my_malloc(sizeof(KLPProjection));

	//	Count the crossings.
	theKLPProjection->num_crossings = CountCrossings();

	//	Count the free loops.
	//	For a hyperbolic link, there won't be any.
	theKLPProjection->num_free_loops = CountFreeLoops();

	//	Copy the number of components.
	theKLPProjection->num_components = itsNumComponents;
	
	//	Allocate the array of crossings.
	theKLPProjection->crossings = (KLPCrossing *) my_malloc(theKLPProjection->num_crossings * sizeof(KLPCrossing));

	//	Assign a KLPCrossing to each LPCrossing.
	AssignKLPCrossings(theKLPProjection->crossings);

	//	Just to be safe, initialize all fields in all KLPCrossings.
	InitKPCCrossings();

	//	Transfer the data describing the link projection from our private
	//	LP format to the kernel's KLP format.
	TransferCrossingData();

	//	Just to be safe, restore the itsKLPCrossing pointers to NULL.
	EraseKLPCrossingPointers();

	return theKLPProjection;
}


Boolean CLinkProjectionPanorama::AllComponentsComplete(void)
{
	LPComponent	*theComponent;

	for (	theComponent = itsComponentList;
			theComponent != NULL;
			theComponent = theComponent->next)
	{
		if (theComponent->firstVertex->backwardEdge == NULL)
			return FALSE;
	}

	return TRUE;
}


int CLinkProjectionPanorama::CountCrossings(void)
{
	int			theCount;
	LPCrossing	*theCrossing;

	theCount = 0;

	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)

		theCount++;

	return theCount;
}


int CLinkProjectionPanorama::CountFreeLoops(void)
{
	//	This implementation of CountFreeLoops() is somewhat inefficient,
	//	in that it retraverses the crossing list once for each component.
	//	However, it is simple and robust, and will most likely require
	//	a negligible amount of time in any case.

	int			theNumFreeLoops;
	LPComponent	*theComponent;
	LPCrossing	*theCrossing;
	Boolean		theComponentHasACrossing;

	theNumFreeLoops = 0;

	for (	theComponent = itsComponentList;
			theComponent != NULL;
			theComponent = theComponent->next)
	{
		theComponentHasACrossing = FALSE;

		for (	theCrossing = itsCrossingList;
				theCrossing != NULL;
				theCrossing = theCrossing->next)
		{
			if (theCrossing-> overstrand->backwardVertex->component == theComponent
			 || theCrossing->understrand->backwardVertex->component == theComponent)
			{
				theComponentHasACrossing = TRUE;
				break;
			}
		}

		if (theComponentHasACrossing == FALSE)
			theNumFreeLoops++;
	}

	return theNumFreeLoops;
}


void CLinkProjectionPanorama::AssignKLPCrossings(
	KLPCrossing	*theKLPCrossings)
{
	LPCrossing	*theCrossing;
	int			theIndex;

	for (	theCrossing = itsCrossingList, theIndex = 0;
			theCrossing != NULL;
			theCrossing = theCrossing->next, theIndex++)

		theCrossing->itsKLPCrossing = &theKLPCrossings[theIndex];
}


void CLinkProjectionPanorama::InitKPCCrossings(void)
{
	LPCrossing	*theCrossing;
	int			i,
				j;

	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)
	{
		for (i = 0; i < 2; i++)
		{
			for (j = 0; j < 2; j++)
			{
				theCrossing->itsKLPCrossing->neighbor[i][j]	= NULL;
				theCrossing->itsKLPCrossing->strand[i][j]	= KLPStrandUnknown;
			}

			theCrossing->itsKLPCrossing->component[i] = -1;
		}

		theCrossing->itsKLPCrossing->handedness = KLPCrossingTypeUnknown;
	}
}


void CLinkProjectionPanorama::TransferCrossingData(void)
{
	LPCrossing	*theCrossing;

	//	Compute the handedness first, because we'll need it
	//	to compute the strand and component data.
	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)
	{
		TransferHandednessData(theCrossing);
	}

	//	Figure out where the x- and y-strands go.
	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)
	{
		TransferStrandData(theCrossing);
	}

	//	Note the components' indices.
	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)
	{
		TransferComponentData(theCrossing);
	}
}


void CLinkProjectionPanorama::TransferHandednessData(
	LPCrossing	*aCrossing)
{
	//	Think of the understrand and the overstrand as vectors.
	//	If the z-component of the cross product
	//
	//				understrand X overstrand
	//
	//	is positive, the KLPCrossingType is KLPHalfTwistCL.
	//	Otherwise it's KLPHalfTwistCCL.
	//
	//	Technical note:  The Mac's (h,v) screen coordinates use
	//	a vertical axis which points down, so the above formula
	//	might be the opposite of what you were expecting.

	Point	underVector,
			overVector;

	underVector.h = aCrossing->understrand-> forwardVertex->position.h
				  - aCrossing->understrand->backwardVertex->position.h;
	underVector.v = aCrossing->understrand-> forwardVertex->position.v
				  - aCrossing->understrand->backwardVertex->position.v;

	overVector.h  = aCrossing->overstrand -> forwardVertex->position.h
				  - aCrossing->overstrand ->backwardVertex->position.h;
	overVector.v  = aCrossing->overstrand -> forwardVertex->position.v
				  - aCrossing->overstrand ->backwardVertex->position.v;

	aCrossing->itsKLPCrossing->handedness =
		(underVector.h*overVector.v - underVector.v*overVector.h > 0) ?
		KLPHalfTwistCL :
		KLPHalfTwistCCL;
		
}


void CLinkProjectionPanorama::TransferStrandData(
	LPCrossing	*aCrossing)
{
	LPCrossing	*theNextCrossing;
	LPEdge		*theNextEdge;

	//	We keep track of a strand by a pair (LPCrossing, LPEdge), so
	//	we can deduce whether it's the understrand or the overstrand.

	//	Find the first crossing in the direction of the positive understrand.
	FindNextCrossing(aCrossing, aCrossing->understrand, &theNextCrossing, &theNextEdge);
	RecordNeighbors (aCrossing, aCrossing->understrand,  theNextCrossing,  theNextEdge);

	//	Find the first crossing in the direction of the positive overstrand.
	FindNextCrossing(aCrossing, aCrossing->overstrand, &theNextCrossing, &theNextEdge);
	RecordNeighbors (aCrossing, aCrossing->overstrand,  theNextCrossing,  theNextEdge);
}


void CLinkProjectionPanorama::FindNextCrossing(
	LPCrossing	*aFromCrossing,
	LPEdge		*aFromEdge,
	LPCrossing	**aToCrossing,
	LPEdge		**aToEdge)
{
	int	theCityBlockDistance;

	//	To measure the relative positions of the crossings along an edge,
	//	we compare the city block distance from each crossing to the
	//	edge's backwardVertex.
	theCityBlockDistance =
		ABS(aFromCrossing->position.h - aFromEdge->backwardVertex->position.h)
	  + ABS(aFromCrossing->position.v - aFromEdge->backwardVertex->position.v);

	//	Check whether there is a crossing further along aFromEdge.
	//	If there isn't one, FindCrossingJustAfter() will return NULL.
	*aToCrossing	= FindCrossingJustAfter(aFromEdge, theCityBlockDistance);
	*aToEdge		= aFromEdge;

	//	If there are no further crossings on this edge, start looking
	//	at subsequent edges until we find one.  Just to be safe, pass
	//	-1 for the minimum city block distance instead of 0.  (Edges aren't
	//	allowed to pass through vertices, but the coordinates of a crossing
	//	are rounded to integers, and might conceivably coincide with a
	//	vertex.  Fortunately edges are excluded from whole neighborhoods
	//	of previous crossings, so we needn't worry about roundoff error
	//	changing the order of crossings.)
	//
	//	Note that we're sure to eventually find aToCrossing.
	//	In the worst case, we follow the component around until we get
	//	back to aFromCrossing.
	while (*aToCrossing == NULL)
	{
		*aToEdge		= (*aToEdge)->forwardVertex->forwardEdge;
		*aToCrossing	= FindCrossingJustAfter(*aToEdge, -1);
	}
}


LPCrossing *CLinkProjectionPanorama::FindCrossingJustAfter(
	LPEdge		*anEdge,
	int			aMinimumDistance)
{
	LPCrossing	*theFirstCrossing,
				*theCrossing;
	int			theFirstDistance,
				theDistance;

	//	Find the first crossing on anEdge which is a distance strictly
	//	greater than aMinimumDistance from the initial vertex.
	//	If none exists, return NULL. 

	theFirstCrossing = NULL;
	theFirstDistance = 0x7FFF;	//	actually we're using 32-bit ints,
								//	but this will do

	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)
	{
		if (theCrossing->understrand != anEdge
		 && theCrossing->overstrand  != anEdge)
			continue;

		theDistance =
			ABS(theCrossing->position.h - anEdge->backwardVertex->position.h)
		  + ABS(theCrossing->position.v - anEdge->backwardVertex->position.v);

		if (theDistance <= aMinimumDistance)
			continue;

		if (theDistance < theFirstDistance)
		{
			theFirstCrossing = theCrossing;
			theFirstDistance = theDistance;
		}
	}

	return theFirstCrossing;
}


void CLinkProjectionPanorama::RecordNeighbors(
	LPCrossing	*aFromCrossing,
	LPEdge		*aFromEdge,
	LPCrossing	*aToCrossing,
	LPEdge		*aToEdge)
{
	KLPStrandType	theFromStrandType,
					theToStrandType;

	theFromStrandType = GetStrandType(aFromCrossing, aFromEdge);
	theToStrandType   = GetStrandType(aToCrossing,   aToEdge);

	aFromCrossing->itsKLPCrossing->neighbor[theFromStrandType][KLPForward ] = aToCrossing->itsKLPCrossing;
	aFromCrossing->itsKLPCrossing->strand  [theFromStrandType][KLPForward ] = theToStrandType;

	aToCrossing  ->itsKLPCrossing->neighbor[ theToStrandType ][KLPBackward] = aFromCrossing->itsKLPCrossing;
	aToCrossing  ->itsKLPCrossing->strand  [ theToStrandType ][KLPBackward] = theFromStrandType;
}


KLPStrandType CLinkProjectionPanorama::GetStrandType(
	LPCrossing	*aCrossing,
	LPEdge		*anEdge)
{
	if (aCrossing->understrand == anEdge)
	{
		switch (aCrossing->itsKLPCrossing->handedness)
		{
			case KLPHalfTwistCL:
				return KLPStrandY;

			case KLPHalfTwistCCL:
				return KLPStrandX;

			default:
				uFatalError("GetStrandType", "CLinkProjectionPanorama");
		}
	}
	else if (aCrossing->overstrand == anEdge)
	{
		switch (aCrossing->itsKLPCrossing->handedness)
		{
			case KLPHalfTwistCL:
				return KLPStrandX;

			case KLPHalfTwistCCL:
				return KLPStrandY;

			default:
				uFatalError("GetStrandType", "CLinkProjectionPanorama");
		}
	}
	else
		uFatalError("GetStrandType", "CLinkProjectionPanorama");

	return KLPStrandUnknown;	//	keep the compiler happy
}


void CLinkProjectionPanorama::TransferComponentData(
	LPCrossing	*aCrossing)
{
	int	theUnderComponentIndex,
		theOverComponentIndex;

	theUnderComponentIndex = aCrossing->understrand->backwardVertex->component->index;
	theOverComponentIndex  = aCrossing-> overstrand->backwardVertex->component->index;

	switch (aCrossing->itsKLPCrossing->handedness)
	{
		case KLPHalfTwistCL:
			aCrossing->itsKLPCrossing->component[KLPStrandX] = theOverComponentIndex;
			aCrossing->itsKLPCrossing->component[KLPStrandY] = theUnderComponentIndex;
			break;

		case KLPHalfTwistCCL:
			aCrossing->itsKLPCrossing->component[KLPStrandX] = theUnderComponentIndex;
			aCrossing->itsKLPCrossing->component[KLPStrandY] = theOverComponentIndex;
			break;

		default:
			uFatalError("TransferComponentData", "CLinkProjectionPanorama");
	}
}


void CLinkProjectionPanorama::EraseKLPCrossingPointers(void)
{
	LPCrossing	*theCrossing;

	for (	theCrossing = itsCrossingList;
			theCrossing != NULL;
			theCrossing = theCrossing->next)

		theCrossing->itsKLPCrossing = NULL;
}


void CLinkProjectionPanorama::FreeKLPProjection(
	KLPProjection	*aKLPProjection)
{
	if (aKLPProjection != NULL)
	{
		if (aKLPProjection->crossings != NULL)
			my_free( aKLPProjection->crossings);

		my_free( aKLPProjection);
	}
}


void CLinkProjectionPanorama::FreeLinkProjection(void)
{
	LPComponent	*theDeadComponent;
	LPCrossing	*theDeadCrossing;

	itsNumComponents = 0;

	while (itsComponentList != NULL)
	{
		theDeadComponent = itsComponentList;
		itsComponentList = itsComponentList->next;
		FreeComponent(theDeadComponent);
	}

	while (itsCrossingList != NULL)
	{
		theDeadCrossing = itsCrossingList;
		itsCrossingList = itsCrossingList->next;
		my_free(theDeadCrossing);
	}

	itsHotVertex = NULL;
}


void CLinkProjectionPanorama::FreeComponent(
	LPComponent	*aComponent)
{
	LPVertex	*theDeadVertex;
	LPEdge		*theDeadEdge;

	//	aComponent may be linear or circular.
	//	If it's circular, remove an LPEdge to make it linear.
	//	(Note:  the second part of the "if" statement makes sure
	//	aComponent is really a loop, and not just a single vertex.)
	if (aComponent->firstVertex == aComponent->lastVertex
	 && aComponent->firstVertex->forwardEdge != NULL)
	{
		theDeadEdge = aComponent->firstVertex->forwardEdge;
		aComponent->firstVertex = theDeadEdge->forwardVertex;
		aComponent->firstVertex->backwardEdge	= NULL;
		aComponent->lastVertex ->forwardEdge	= NULL;
		my_free( theDeadEdge);
	}

	//	The component is now linear.  Remove segments one at a time
	//	until only a single LPVertex remains.
	while (aComponent->firstVertex != aComponent->lastVertex)
	{
		theDeadVertex			= aComponent->firstVertex;
		theDeadEdge				= theDeadVertex->forwardEdge;
		aComponent->firstVertex	= theDeadEdge->forwardVertex;
		aComponent->firstVertex->backwardEdge = NULL;
		my_free( theDeadVertex);
		my_free( theDeadEdge);
	}

	//	Free the last LPVertex.
	my_free( aComponent->firstVertex);

	//	Free the LPComponent itself.
	my_free( aComponent);
}

void CLinkProjectionPanorama::CalculateEquation(
	LPEdge	*aNewEdge)
{
	//	Recall from the comments in struct LPEdge that the equation
	//	for the line containing aNewEdge is
	//
	//						a*x + b*y = c
	//	where
	//					a = pt2.v - pt1.v
	//					b = pt1.h - pt2.h
	//					c = pt1.h * pt2.v  -  pt1.v * pt2.h

	Point	pt1,
			pt2;

	pt1 = aNewEdge->backwardVertex->position;
	pt2 = aNewEdge->forwardVertex ->position;

	aNewEdge->a = pt2.v - pt1.v;
	aNewEdge->b = pt1.h - pt2.h;
	aNewEdge->c = pt1.h * pt2.v  -  pt1.v * pt2.h;
}

void CLinkProjectionPanorama::MaybeCreateOneCrossing(
	LPEdge	*anOverstrand,
	LPEdge	*anUnderstrand)
{
	//	Proposition.  Two edges E1 and E2 intersect iff the vertices of E1
	//	lie on opposite sides of the line containing E2 and the vertices
	//	of E2 lie on opposite sides of the line containing E1.
	//
	//	Proof.  Exercise for the reader.

	LPEdge		*e1,
				*e2;
	Point		p1,
				p2,
				q1,
				q2;
	LPCrossing	*theNewCrossing;
	double		theDet;

	e1 = anOverstrand;
	e2 = anUnderstrand;

	//	We assume the general position code has done its thing before
	//	we check for crossings, so we may safely assume that one edge
	//	cannot contain a vertex of another except at a common endpoint.
	if (e1->backwardVertex	== e2->backwardVertex
	 || e1->backwardVertex	== e2->forwardVertex
	 || e1->forwardVertex	== e2->backwardVertex
	 || e1->forwardVertex	== e2->forwardVertex)
		return;

	p1 = e1->backwardVertex->position;
	q1 = e1->forwardVertex ->position;

	p2 = e2->backwardVertex->position;
	q2 = e2->forwardVertex ->position;

	//	Recall that the coordinates of the vertices will lie (approximately)
	//	in the range 0 to 1000.  The following comparisions all involve
	//	binominal expressions in the original coordinates, so there is no
	//	danger of overflowing the 32-bit integer arithmetic.

	if
	(
		//	p2 and q2 lie on the same side of e1
		(
			(e1->a * p2.h  +  e1->b * p2.v  >  e1->c)
		 ==
			(e1->a * q2.h  +  e1->b * q2.v  >  e1->c)
		)
	 ||
		//	p1 and q1 lie on the same side of e2
		(
			(e2->a * p1.h  +  e2->b * p1.v  >  e2->c)
		 ==
			(e2->a * q1.h  +  e2->b * q1.v  >  e2->c)
		)
	)
		return;

	//	The edges intersect.
	//	Find the point of intersection and create the LPCrossing.
	//	To find the coordinates of the intersection, solve
	//
	//						a1*x + b1*y = c1
	//						a2*x + b2*y = c2
	//
	//	by writing it as a matrix equation
	//
	//					( a1  b1 ) ( x )     ( c1 )
	//					(        ) (   )  =  (    )
	//					( a2  b2 ) ( y )     ( c2 )
	//
	//	and inverting the matrix.
	//
	//					( x )      1  ( b2 -b1 ) ( c1 )
	//					(   )  =  --- (        ) (    )
	//					( y )     det (-a2  a1 ) ( c2 )

	theNewCrossing = (LPCrossing *) my_malloc(sizeof(LPCrossing));

	//	Note:  a and b are "first order" expressions in the coordinates
	//	of the endpoints, so we may safely multiply them together without
	//	fear of overflowing the 32-bit integer arithmetic.
	theDet = (double) (e1->a * e2->b  -  e1->b * e2->a);

	//	Now we swtich over to floating point arithmetic, to avoid computing
	//	trinomials which might overflow the 32-bit integer arithmetic.
	theNewCrossing->position.h =
		(	((double) e2->b * (double) e1->c)
		  - ((double) e1->b * (double) e2->c) )  /  theDet;
	theNewCrossing->position.v =
		(	((double) e1->a * (double) e2->c)
		  - ((double) e2->a * (double) e1->c) )  /  theDet;

	SetRect(&theNewCrossing->neighborhood,
			(int)theNewCrossing->position.h - CROSSING_EPSILON,
			(int)theNewCrossing->position.v - CROSSING_EPSILON,
			(int)theNewCrossing->position.h + CROSSING_EPSILON,
			(int)theNewCrossing->position.v + CROSSING_EPSILON);

	theNewCrossing->overstrand	= anOverstrand;
	theNewCrossing->understrand	= anUnderstrand;

	theNewCrossing->itsKLPCrossing = NULL;

	theNewCrossing->next	= itsCrossingList;
	itsCrossingList			= theNewCrossing;
}

void CLinkProjectionPanorama::ContentsToWindow(
	CLinkProjectionWrapper	*aContents)
{
	LPComponent	**theComponents;
	LPVertex	**theVertices;
	LPEdge		**theEdges;
	LPVertex	*theVertex;
	int			i;

	FreeLinkProjection();

	if (aContents == NULL)
		return;

	//	Create temporary arrays to hold pointers to the components,
	//	vertices, edges.  (The crossings will be handled differently.)
	theComponents	= (LPComponent **) my_malloc(aContents->itsNumComponents * sizeof(LPComponent *));
	theVertices		= (LPVertex **)    my_malloc(aContents->itsNumVertices   * sizeof(LPVertex *));
	theEdges		= (LPEdge **)      my_malloc(aContents->itsNumEdges      * sizeof(LPEdge *));

	//	Allocate the data structures used to hold the LP data.
	for (i = 0; i < aContents->itsNumComponents; i++)
		theComponents[i] = (LPComponent *) my_malloc(sizeof(LPComponent));
	for (i = 0; i < aContents->itsNumVertices;   i++)
		theVertices[i]   = (LPVertex *)    my_malloc(sizeof(LPVertex));
	for (i = 0; i < aContents->itsNumEdges;      i++)
		theEdges[i]      = (LPEdge *)      my_malloc(sizeof(LPEdge));

	//	Fill in the data.

	itsComponentList	= (aContents->itsNumComponents > 0) ? theComponents[0] : NULL;
	itsNumComponents	= aContents->itsNumComponents;
	itsHotVertex		= (aContents->itsHotVertex >= 0) ? theVertices[aContents->itsHotVertex] : NULL;
	itsCrossingList		= NULL;

	for (i = 0; i < aContents->itsNumComponents; i++)
	{
		theComponents[i]->index			= i;
		theComponents[i]->firstVertex	= theVertices[aContents->itsFirstVertices[i]];
		theComponents[i]->lastVertex	= theVertices[aContents->itsLastVertices [i]];
		theComponents[i]->next			= (i+1 < aContents->itsNumComponents) ? theComponents[i+1] : NULL;
	}

	for (i = 0; i < aContents->itsNumVertices; i++)
	{
		theVertices[i]->position.h		= aContents->itsHCoordinates[i];
		theVertices[i]->position.v		= aContents->itsVCoordinates[i];
		theVertices[i]->component		= NULL;		//	will be set below
		theVertices[i]->backwardEdge	= NULL;		//	will be set by edge (if any)
		theVertices[i]->forwardEdge		= NULL;		//	will be set by edge (if any)
		theVertices[i]->index			= -1;		//	used only locally
	}

	for (i = 0; i < aContents->itsNumEdges; i++)
	{
		theEdges[i]->backwardVertex	= theVertices[aContents->itsBackwardVertices[i]];
		theEdges[i]->forwardVertex	= theVertices[aContents->itsForwardVertices [i]];
		CalculateEquation(theEdges[i]);
		theEdges[i]->index			= -1;		//	used only locally

		theEdges[i]->backwardVertex->forwardEdge = theEdges[i];
		theEdges[i]->forwardVertex->backwardEdge = theEdges[i];
	}

	//	Let MaybeCreateOneCrossing() do the work of reconstructing
	//	the crossings.  In this case "Maybe" is a misnomer -- the
	//	two edges are known to cross.
	for (i = 0; i < aContents->itsNumCrossings; i++)
		MaybeCreateOneCrossing(	theEdges[aContents->itsOverstrands [i]],
								theEdges[aContents->itsUnderstrands[i]]);

	//	Tell each vertex which component it belongs to.
	for (i = 0; i < aContents->itsNumComponents; i++)
	{
		theVertex = theComponents[i]->firstVertex;
		while (TRUE)
		{
			theVertex->component = theComponents[i];
			if (theVertex->forwardEdge == NULL)
				break;
			theVertex = theVertex->forwardEdge->forwardVertex;
			if (theVertex->component != NULL)
				break;
		}
	}

	//	Free the temporary arrays.
	my_free(theComponents);
	my_free(theVertices);
	my_free(theEdges);

}
