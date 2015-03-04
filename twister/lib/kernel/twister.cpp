#include <vector>
#include <queue>
#include <iostream>
#include <sstream>
#include "global.h"
#include "twister.h"

///////////// Permutation code //////////////////

perm::perm()
{
	// By default, a permutation is the identity permutation.
	for (int i = 0; i < 4; i++)
		image[i] = i;
}

perm::perm(const perm &to_copy)
{
	for (int i = 0; i < 4; i++)
		image[i] = to_copy.image[i];
}

perm::perm(const int* image_in)
{
	for (int i = 0; i < 4; i++)
		image[i] = image_in[i];
}

int perm::operator[](int input) const
{
	return image[input];
}

perm perm::inverse() const
{
	perm outp;
	
	for (int i = 0; i < 4; i++)
		outp.image[image[i]] = i;
	
	return outp;
}

// Composition of permutations
perm perm::of(const perm &other) const
{
	perm outp;
	
	for (int i=0; i<4; i++)
		outp.image[i] = image[other.image[i]];
	
	return outp;
}

///////////// Tetrahedron code //////////////////

tetra::tetra(manifold *M, Category mycategory, Position myposition, Side myside, int mylayer, tetra *my_parent)
{
	output_debugging("tet", 3);
	
	for (int i = 0; i < 4; i++)
		gluedto[i] = NULL;
	
	// There's nothing canonical to put in gluing[],
	// so we put in the identity permutation.
	for (int i = 0; i < 4; i++)
		gluing[i] = perm(1,0,2,3);
	
	next = NULL;
	prev = NULL;
	
	parent = my_parent;
	if (my_parent != NULL) my_parent->child = this;
	child = NULL;
	
	M->onemore(this);  // Add this new tetra to the doubly linked list of tetras.
	
	category = mycategory;
	position = myposition;
	side = myside;
	layer = mylayer;
	
	snap_index = -1;  // Index only determined when snap_print is called.
	cusp_number = -1;  // The cusp number of vertex 3 because vertices 0, 1 and 2 are always material.
	temp = -1;
	for (int i = 0; i < 4; i++)
	{
		peripheral_curves[longitude][i] = 0;
		peripheral_curves[meridian][i] = 0;
	}
}

tetra::~tetra()
{
	output_debugging("detet", 3);
}

// The function gluesym() glues together two tetrahedra
// along faces - which it checks are free.
void tetra::gluesym(tetra *whereglue, int whichface, const perm &how)
{
	if ( (gluedto[whichface] != NULL) || (whereglue->gluedto[how[whichface]] != NULL) )
		output_error("Invalid gluing. Possible causes; self-intersecting curves, intersecting 2-handles.");
	
	this->set_gluedto(whichface, whereglue);
	this->set_gluing(whichface,how);
	
	whereglue->set_gluedto(how[whichface], this);
	whereglue->set_gluing(how[whichface], how.inverse());
}

// The function ungluesym() breaks appart two tetrahedra that 
// are glued together via face whichface of this.
void tetra::ungluesym(int whichface)
{
	if (gluedto[whichface] == NULL) 
		return;  // This is correct as it simplifies later loops.
	
	if (gluedto[whichface]->get_gluedto(gluing[whichface][whichface]) == NULL)
		output_error("Non-symmetric gluing detected.");
	
	gluedto[whichface]->set_gluedto(gluing[whichface][whichface], NULL);  // Unglue him from me.
	this->set_gluedto(whichface, NULL);  // Then unglue me from him.
}

// The function subbedby() breaks open a gluing and replaces
// one of the tetrahedra (A below, glued to B above) by another 
// (D below, now glued to B above like A used to be).
// 
// For oriented manifolds, the permutation "how" should be even.
// 
// It's okay if some of the tetrahedra A, B, and D coincide, but
// the three faces involved must be distinct.
void tetra::subbedby(tetra *whereglue, int whichface, const perm &how)
{
	// before:    A-B   D free
	// after:     A free   D-B
	
	// A: this
	// D: whereglue
	// B: *gluedto[whichface]
	
	// "how" maps vertices of A to vertices of D
	// (Identity means glue D to B just like A was glued)
	// "how" should be an even permutation.
	
	if (gluedto[whichface] == NULL)
		return;
	if (whereglue->gluedto[how[whichface]] != NULL)
		output_error("Invalid subbedby. Attempted to substitute non-free face.");
	
	perm AtoB = gluing[whichface];
	
	tetra *B = gluedto[whichface];
	ungluesym(whichface);  // Separate A & B.
	B->gluesym(whereglue, AtoB[whichface], how.of(AtoB.inverse()));  // Glue B to D.
}

// Given a non-3 fromface, walk_about the edge meeting face fromface and face 3 and then glue.
static const perm perm_walk_about_list[3]= {perm(3,1,2,0), perm(0,3,2,1), perm(0,1,3,2)};

void tetra::walk_about(int fromface)
{
	output_debugging("walk", 3);
	
	if (gluedto[fromface] != NULL)
		return;
	
	// Here's how we start.
	tetra *currenttetra = this;
	perm init_how = perm_walk_about_list[fromface];
	perm how = init_how;
	
	// Now walk around the edge.
	while ((currenttetra->gluedto[how[fromface]]) != NULL)
	{
		int dir = how[fromface];
		how = (currenttetra->get_gluing(dir)).of(how.of(init_how));  // Composition of 3 odd permutations is also odd.
		currenttetra = currenttetra->get_gluedto(dir);
	}
	
	// We've found the last tetra next to the edge
	// so glue it to the first one.
	gluesym(currenttetra, fromface, how);
}

///////////// Cube code //////////////////

cube::cube(manifold &home_in, int layer_number, cube *parent_cube)
{
	output_debugging("cube", 3);
	
	home = &home_in;
	status = true;
	
	is_glued[0] = is_glued[1] = false;
	
	topleft = new tetra(home, cubie, top, left, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(top, left) : NULL);
	midleft = new tetra(home, cubie, mid, left, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(mid, left) : NULL);
	lowleft = new tetra(home, cubie, low, left, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(low, left) : NULL);
	
	topright = new tetra(home, cubie, top, right, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(top, right) : NULL);
	midright = new tetra(home, cubie, mid, right, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(mid, right) : NULL);
	lowright = new tetra(home, cubie, low, right, layer_number, (parent_cube != NULL) ? parent_cube->get_tetra(low, right) : NULL);
	
	// Assemble these tetrahedra into a cube.
	// We'll glue these into two triangular prisms, and then a cube.
	// See the diagram of the cube.
	
	// Gluings inside left prism.
	topleft->gluesym(midleft, 3, perm(0,1,3,2));
	midleft->gluesym(lowleft, 3, perm(3,1,2,0));
	
	// Gluings inside right prism.
	topright->gluesym(midright, 3, perm(0,3,2,1));
	midright->gluesym(lowright, 3, perm(3,1,2,0));
	
	// Gluings between the two prisms.
	midleft->gluesym(topright, 1, perm(1,0,2,3));
	lowleft->gluesym(midright, 1, perm(1,0,2,3));
}

// Glues two cubes together in the standard direction, orienations given by u1 and u2.
void glue_cubes(cube *c1, cube *c2, bool u1, bool u2)
{
	// Depending on the orientations of c1 and c2, glue some
	// vertical face of cube c1 to some vertical face of cube c2.
	int toface = u1 ? 1 : 0;
	perm howto = u1 ? perm(0,2,1,3) : perm(2,1,0,3);
	
	if (u1 && u2)
	{
		(c1->get_tetra(low, right))->gluesym(c2->get_tetra(low, left), toface, howto);
		(c1->get_tetra(top, right))->gluesym(c2->get_tetra(top, left), toface, howto);
	}
	else if (u1 && !u2)
	{
		(c1->get_tetra(low, right))->gluesym(c2->get_tetra(low, right), toface, howto);
		(c1->get_tetra(top, right))->gluesym(c2->get_tetra(mid, right), toface, howto);
	}
	else if (!u1 && u2)
	{
		(c1->get_tetra(mid, left))->gluesym(c2->get_tetra(low, left), toface, howto);
		(c1->get_tetra(top, left))->gluesym(c2->get_tetra(top, left), toface, howto);
	}
	else if (!u1 && !u2)
	{
		(c1->get_tetra(mid, left))->gluesym(c2->get_tetra(low, right), toface, howto);
		(c1->get_tetra(top, left))->gluesym(c2->get_tetra(mid, right), toface, howto);
	}
	
	return;
}

// Replaces the child cube of a given cube with a gadget. This is identical 
// to the original cube on 4 of its sides while the remaining two (opposite) sides 
// have a tunnel connecting them, in the direction specified by upright. 
// A gadget is constructed from: the original cube, an exact copy of 
// that cube and two triangular prisms called gadget-pups (pup as in pup-tent).

//  #-----------------#   #-----------------#
//  |\---\       /---/|   |\---\       /---/|
//  |     \--|--/   / |   | \   \--|--/     |
//  |       /|     /  |   |  \     |\       |
//  |      /{X}   /   |   |   \   {X}\      |
//  |     /{XXX} /    |   |    \ {XXX}\     |
//  |    / {XXX}/     |   |     \{XXX} \    |
//  |   /   {X}/      |   |      \{X}   \   |
//  |  /     |/       |   |       \|     \  |
//  | /   /--|--\     |   |     /--|--\   \ |
//  |/---/       \---\|   |/---/       \---\|
//  #-----------------#   #-----------------#
//         Front                  Back

// This is what the front and back of the gadget look like. The X's denote the
// tunnel running throught the center of the cube. The gadget-pups can be see at
// the top and bottom of each picture. The original cube is on the left (from the 
// front view) while the clone of it is on the right.

cube *convert_cube_to_gadget(cube *C, bool upright)
{
	output_debugging("convert", 3);
	
	manifold *M = C->home;
	int layer_number = M->get_num_layers();  // We will need this a lot, so better remember it.
	
	cube *C2 = new cube(*M, layer_number);
	
	tetra *topleft = C->get_tetra(top, left)->get_child();
	tetra *midleft = C->get_tetra(mid, left)->get_child();
	tetra *lowleft = C->get_tetra(low, left)->get_child();
	tetra *topright = C->get_tetra(top, right)->get_child();
	tetra *midright = C->get_tetra(mid, right)->get_child();
	tetra *lowright = C->get_tetra(low, right)->get_child();
	
	tetra *topleft2 = C2->get_tetra(top, left);
	tetra *midleft2 = C2->get_tetra(mid, left);
	tetra *lowleft2 = C2->get_tetra(low, left);
	tetra *topright2 = C2->get_tetra(top, right);
	tetra *midright2 = C2->get_tetra(mid, right);
	tetra *lowright2 = C2->get_tetra(low, right);
	
	topleft2->set_parent(topleft);
	midleft2->set_parent(midleft);
	lowleft2->set_parent(lowleft);
	topright2->set_parent(topright);
	midright2->set_parent(midright);
	lowright2->set_parent(lowright);
	
	for (int i = 0; i < 6; i++)
	{
		C->get_tetra(i)->get_child()->set_category(gadget_left);  // Remember, the old pieces are now part of the gadget too.
		C2->get_tetra(i)->set_category(gadget_right);  // Remember, the new pieces are now part of the gadget too.
	}
	
	// Subbedby the RHS of the duplicant in.
	int toglue = (upright) ? 2 : 1;  // Which face is on the right depends on the orientation.
	tetra *side_tetra = (upright) ? midright : topright;
	tetra *side_tetra2 = (upright) ? midright2 : topright2;
	side_tetra->subbedby(side_tetra2, toglue, perm(0,1,2,3));
	
	side_tetra = (upright) ? lowright : lowright;
	side_tetra2 = (upright) ? lowright2 : lowright2;
	side_tetra->subbedby(side_tetra2, toglue, perm(0,1,2,3));
	
	// These will form the gadget_pups above & below.
	tetra *newfront = new tetra(M, gadget_pup, top, neither, layer_number, topleft2);
	tetra *newmid = new tetra(M, gadget_pup, top, neither, layer_number, midleft2);
	tetra *newback = new tetra(M, gadget_pup, top, neither, layer_number, lowleft2);
	tetra *newfront2 = new tetra(M, gadget_pup, low, neither, layer_number, topright2);
	tetra *newmid2 = new tetra(M, gadget_pup, low, neither, layer_number, midright2);
	tetra *newback2 = new tetra(M, gadget_pup, low, neither, layer_number, lowright2);
	
	/*
	---================= Travellers Beware ==================---
	---=== You are now entering the Land of Hard - Coding ===---
	
	Over the next 100 lines or so we will build and insert the 
	upper and lower pup-gadgets. These individually hand 
	crafted pup-gadgets have 32+ hard-coded gluings and permutations.
	
	Due to its extreme nature, this ride is NOT suitable for 
	guests with heart conditions, vertigo or high blood pressure.
	
	Please keep your arms and legs inside the vehicle at ALL times.
	---================ You have been warned =================---
	*/
	
	// Pup-gadgets have been designed this way for a good reason. See (*1*).
	if (upright)
	{
		// First the upper gadget_pup.
		
		// Glue up the gadget pup.
		newfront->gluesym(newmid, 1, perm(0,1,3,2));
		newmid->gluesym(newback, 2, perm(1,2,3,0));
		
		// Subbedby the gadget pup into the top.
		topleft->subbedby(newfront, 1, perm(2,3,0,1));
		topright->subbedby(newback, 2, perm(3,0,2,1));
		
		// Glue down gadget pup left side (it's a good job that subbedby gave us some free faces).
		newfront->gluesym(topleft, 2, perm(2,3,1,0));
		newmid->gluesym(topright, 3, perm(3,0,1,2));
		
		// Glue down gadget pup right side.
		newmid->gluesym(topleft2, 0, perm(1,2,3,0));
		newback->gluesym(topright2, 1, perm(1,2,3,0));
		
		// Now for the bottom gadget_pup.
		
		// Glue up the gadget pup.
		newfront2->gluesym(newmid2, 2, perm(0,2,1,3));
		newmid2->gluesym(newback2, 3, perm(3,2,0,1));
		
		// Subbedby the gadget pup into the bottom.
		lowleft->subbedby(newfront2, 3, perm(1,2,0,3));
		lowright->subbedby(newback2, 3, perm(1,0,3,2));
		
		// Glue down gadget pup left side (it's a good job that subbedby gave us some free faces).
		newfront2->gluesym(lowleft, 1, perm(2,3,1,0));
		newmid2->gluesym(lowright, 2, perm(2,0,3,1));
		
		// Glue down gadget pup right side.
		newmid2->gluesym(lowleft2, 0, perm(3,2,0,1));
		newback2->gluesym(lowright2, 3, perm(1,0,2,3));
	}
	else
	{
		// First the upper gadget_pup.
		
		// Glue up the gadget pup.
		newfront->gluesym(newmid, 2, perm(0,2,1,3));
		newmid->gluesym(newback, 3, perm(3,2,0,1));
		
		// Subbedby the gadget pup into the top.
		topleft->subbedby(newback, 1, perm(0,2,3,1));
		topright->subbedby(newfront, 2, perm(2,1,3,0));
		
		// Glue down gadget pup left side (it's a good job that subbedby gave us some free faces).
		newback->gluesym(topleft, 3, perm(0,3,2,1));
		newmid->gluesym(topright, 0, perm(2,3,1,0));
		
		// Glue down gadget pup right side.
		newmid->gluesym(topleft2, 2, perm(2,3,1,0));
		newfront->gluesym(topright2, 1, perm(3,2,0,1));
		
		// Now for the bottom gadget_pup.
		
		// Glue up the gadget pup.
		newfront2->gluesym(newmid2, 1, perm(0,1,3,2));
		newmid2->gluesym(newback2, 2, perm(1,2,3,0));
		
		// Subbedby the gadget pup into the bottom.
		lowleft->subbedby(newback2, 3, perm(0,3,1,2));
		lowright->subbedby(newfront2, 3, perm(1,2,0,3));
		
		// Glue down gadget pup left side (it's a good job that subbedby gave us some free faces).
		newback2->gluesym(lowleft, 1, perm(0,3,2,1));
		newmid2->gluesym(lowright, 0, perm(3,2,0,1));
		
		// Glue down gadget pup right side.
		newmid2->gluesym(lowleft2, 3, perm(2,1,0,3));
		newfront2->gluesym(lowright2, 2, perm(2,0,3,1));
	}
	
	/*
	---=== You are now leaving the Land of Hard - Coding ===---
		 Please ensure you have all your belongings with you.
				We hope you enjoyed your stay.
	*/
	
	return C2;
}

///////////// Annulus code //////////////////

// Builds a thickend annulus.

annulus::annulus(const std::vector<cube*> &sq_in, const std::vector<bool> &upright_in)
{
	length = int(sq_in.size());
	
	sq = new cube*[length];
	upright = new bool[length];
	
	for (int i = 0; i < length; i++)
	{
		sq[i] = sq_in[i];
		upright[i] = upright_in[i];
	}
	
	// Check for self-intersection.
	for (int i = 0; i < length; i++)
		for (int j = i + 1; j < length; j++)
			if (sq[i] == sq[j])
				output_error("An annulus must not intersect itself.");
	
	// Now glue the cubes together.
	for (int i = 0; i < length; i++)
		glue_cubes(sq[i], sq[(i+1) % length], upright[i], upright[(i+1) % length]);
}

annulus::annulus(const annulus &a)
{
	length = a.length;
	
	sq = new cube*[length];
	upright = new bool[length];
	
	for (int i = 0; i < length; i++)
	{
		sq[i] = a.sq[i];
		upright[i] = a.upright[i];
	}
}

annulus::~annulus()
{
	delete[] sq;
	delete[] upright;
}

// Constructs a 2-handle, and glues it
// above or below this thickened annulus.
void annulus::twohandle(bool is_above)
{
	output_debugging("handle", 2);
	
	manifold *M = sq[0]->home;
	
	if (M->get_manifold_type() == bundle)
		output_error("Cannot attach two-handles to a bundle.");
	
	if (M->get_num_layers() != 0)
		output_error("Handles must be attached before any twisting or drilling occurs.");
	
	tetra *far=NULL, *mid=NULL, *near=NULL;
	
	for(int i = 0; i < length; i++)
	{
		// Glue together a bunch of prisms and
		// attach them to the tops of the cubes.
		
		// construct three tetras
		far = new tetra(M, handle);
		mid = new tetra(M, handle);
		near = new tetra(M, handle);
		
		// Build the prism
		
		// XOR since up/down inverts uprightness
		far->gluesym(mid, 3, (upright[i] ^ is_above) ? perm(0,3,2,1) : perm(0,1,3,2));
		mid->gluesym(near, 3, perm(3,1,2,0));
		
		// Glue prism to the floor / ceiling
		if (is_above)
		{
			far->gluesym(sq[i]->topleft, 0, perm(1,0,2,3));
			mid->gluesym(sq[i]->topright, 0, (upright[i]) ? perm(2,1,0,3) : perm(2,0,3,1));
		}
		else
		{
			far->gluesym(sq[i]->lowleft, 0, perm(3,2,0,1));
			mid->gluesym(sq[i]->lowright, 0, (upright[i]) ? perm(3,0,1,2) : perm(3,2,0,1));
		}
	}
	
	// The next bit of code depends on the values of
	// "far", "mid", and "near" left over from the
	// end of the previous 'for' loop.
	
	tetra *farprev, *midprev, *nearprev;
	int toface = (is_above) ? 2 : 1;
	for (int i = 0; i < length; i++)
	{
		int prev = (i + length - 1) % length;
		
		// Play pass the pointer.
		farprev = far;
		midprev = mid;
		nearprev = near;
		
		far = (is_above) ? sq[i]->topleft->get_gluedto(1) : sq[i]->lowleft->get_gluedto(3);
		mid = far->get_gluedto(3);
		near = mid->get_gluedto(3);
		
		// Glue the prisms to each other.
		near->gluesym(nearprev, toface, perm(0,2,1,3));
		
		(upright[i] ? far : mid)->gluesym((upright[prev] ? midprev : farprev), toface, perm(0,2,1,3));
	}
	
	return;
}

// Recalling how we parse words (see main in twister_main.cpp) 
// a positive twist below corresponds to a left Dehn twist.

// This inserts a "speed bump", i.e. solid 
// torus, above the surface cross interval.
void annulus::twist(Sign whichway)
{
	// There's a much worse way to do this, which for
	// obscenely large values of "length" (>6?)
	// actually uses fewer tetrahedra.
	// (The idea is: retriangulate the annulus so that 
	// it has length 1 and then do the twist(s).  Then
	// undo the retriangulation.  But this is pointless, 
	// SnapPea will totally change the triangulation, 
	// anyway.   A more general notion is to perform a long 
	// twist by conjugating it to a short one...)
	
	// Every speed bump consists of "length" identical layers.
	// Don't peel them apart.  
	// Never look inside a speed bump.  Trust us.
	output_debugging("twist", 2);
	
	manifold *M = sq[0]->home;
	
	for (int i = 0; i < length; i++)
		sq[i]->set_status(false);
	
	std::vector< std::vector<tetra *> > tet;
	int i, j;
	
	// Create a bunch of tetrahedra. We'll refer
	// to tet[i][j] as the new tetrahedron in the
	// i^th row (from the bottom) and j^th column
	// (from the left). 
	
	for( i = 0; i < length; i++ )
	{
		std::vector<tetra *> next_row;
		for( j = 0; j < length; j++ )
			next_row.push_back(new tetra(M, speedbump));
		
		tet.push_back(next_row);
	}
	
	// First thing to do is glue all tetras in a column
	// along their "blue faces".  (For terminology, refer
	// to a certain paper tablecloth in the new Thai
	// restaurant on Bancroft above Telegraph.)
	// Thus loop over j first (sorry 'bout that). 
	// While we're at it we'll attach the last 
	// and first blue face to the ceiling and floor.
	//
	// Keep gluing 'til you're blue in the face.
	
	for( j = 0; j < length; j++ ) 
	{
		for( i = 0; i < length - 1; i++ )
			tet[i][j]->gluesym(tet[i+1][j], 2, perm(0,1,3,2));
		
		// Attach the top of the column to the ceiling and the bottom of the column to the floor.
		sq[j]->topright->subbedby(tet[length-1][j], 2, (upright[j]) ? perm(0,1,2,3) : perm(1,3,2,0));
		sq[j]->topright->gluesym(tet[0][j], 2, (upright[j]) ? perm(0,1,3,2) : perm(1,2,3,0));
	}
	
	// That finishes off the blue faces.  
	// Now to glue neighboring columns together
	// along the "red" faces. While we're at it
	// we'll attach the left-over red faces to the 
	// ceiling and floor. This is a bit tricky, so 
	// watch carefully. Follow the red card...
	
	for( j = 0; j < length; j++ )
	{
		int next = (j+1) % length;
		int begin, end;
		
		// Assuming here that true == 1 and false == 0.
		int offset = whichway*(!upright[j] + upright[next]);
		
		if (whichway == plus)
		{
			begin = offset;
			end = length;
		}
		else
		{
			begin = 0;
			end = length + offset;
		}
		
		// Gluing the columns together.
		for( i = begin; i < end; i++ )
			tet[i][j]->gluesym(tet[i - offset][next], 1, perm(1,0,2,3));
		
		// Now glue the left over red faces on top to the 
		// ceiling using subbedby. Then glue remaining red
		// faces (on the bottom) to the newly available floor.
		
		if ( upright[j] && upright[next] )
		{
			// The offset is whichway*1.
			if (whichway == plus)
			{
				sq[next]->topleft->subbedby(tet[length - 1][next], 1, perm(1,0,3,2));
				sq[next]->topleft->gluesym(tet[0][j], 1, perm(0,1,3,2));
			}
			else
			{
				sq[next]->topleft->subbedby(tet[length - 1][j], 1, perm(0,1,2,3));
				sq[next]->topleft->gluesym(tet[0][next], 1, perm(1,0,2,3));
			}
		}      
		if ( upright[j] && !upright[next] )
		{
			// The offset is zero, so there are 
			// no left-over red faces.  :)
		}      
		
		if ( !upright[j] && upright[next] )
		{
			// The offset is whichway*2.
			if (whichway == plus)
			{
				sq[next]->topleft->subbedby(tet[length - 1][next], 1, perm(1, 0, 3, 2));
				sq[next]->topleft->gluesym(tet[1][j], 1, perm(0,1,3,2));
				
				sq[j]->topleft->subbedby(tet[length - 2][next], 1, perm(2, 0, 1, 3));
				sq[j]->topleft->gluesym(tet[0][j], 1, perm(2,1,0,3));
			}
			else
			{
				sq[j]->topleft->subbedby(tet[length - 1][j], 1, perm(3, 1, 0, 2));
				sq[j]->topleft->gluesym(tet[1][next], 1, perm(3, 0, 1, 2));
				
				sq[next]->topleft->subbedby(tet[length - 2][j], 1, perm(0,1,2,3));
				sq[next]->topleft->gluesym(tet[0][next], 1, perm(1,0,2,3));
			}
		}      
		if ( !upright[j] && !upright[next] )
		{
			// The offset is whichway*1.
			if (whichway == plus)
			{
				sq[j]->topleft->subbedby(tet[length - 1][next], 1, perm(2, 0, 1, 3));
				sq[j]->topleft->gluesym(tet[0][j], 1, perm(2,1,0,3));
			}
			else
			{
				sq[j]->topleft->subbedby(tet[length - 1][j], 1, perm(3, 1, 0, 2));
				sq[j]->topleft->gluesym(tet[0][next], 1, perm(3, 0, 1, 2));
			}
		}      
		// Isn't that all perfectly clear?
		// No?  Well, is my face is red!
	}
}    

// drill() inserts a new layer and 'drills' a tunnel 
// down the core of this thickened annulus' youngest child.
//
// Requires 12*length extra tetra and produces 4*length 
// additional free faces.
void annulus::drill()
{
	output_debugging("drill", 2);
	
	manifold *M = sq[0]->home;
	
	// First find out if we need a new layer.
	bool need_new_layer = false;
	for (int i = 0; i < length; i++)
		if (! sq[i]->get_status())
		{
			need_new_layer = true;
			break;
		}
	
	if (need_new_layer) M->insert_layer();
	
	for (int i = 0; i < length; i++)
		sq[i]->set_status(false);
	
	// For each cube in the thickened annulus replace its child (in the new layer) by a gadget.
	std::vector<cube*> new_cubes;
	for (int i = 0; i < length; i++)
		new_cubes.push_back(convert_cube_to_gadget(sq[i], upright[i]));
	
	// Gadget pups can't get horizontal gluing information from their parents.
	// So we must provide it now.
	for (int i = 0; i < length; i++)
	{
		int j = (i + 1) % length;
		// Glue pup_gadgets together.
		tetra *current_top_back = new_cubes[i]->get_tetra(low, left)->get_child();
		tetra *next_top_front = new_cubes[j]->get_tetra(top, left)->get_child();
		tetra *current_bottom_back = new_cubes[i]->get_tetra(low, right)->get_child();
		tetra *next_bottom_front = new_cubes[j]->get_tetra(top, right)->get_child();
		
		// (*1*) Gadgets were designed so this was true. :)
		current_top_back->gluesym(next_top_front, 0, perm(0,1,3,2));
		current_bottom_back->gluesym(next_bottom_front, 0, perm(0,1,3,2));
		
		glue_cubes(new_cubes[i], new_cubes[j], upright[i], upright[j]);
	}
	
	for (int i = 0; i < length; i++)
		delete new_cubes[i];
	
	return;
}

///////////// Rectangle code ///////////////////////

// Builds a thickened rectangle.

// Really annuli and rectangles should inherit from a strip class.
// But this works so there's no need.

rectangle::rectangle(const std::vector<cube*> &sq_in, const std::vector<bool> &upright_in)
{
	front = -1;
	back = -1;
	one_ended = false;
	
	length = int(sq_in.size());
	
	sq = new cube*[length];
	upright = new bool[length];
	
	for (int i = 0; i < length; i++)
	{
		sq[i] = sq_in[i];
		upright[i] = upright_in[i];
	}
	
	// Check for self-intersection.
	for (int i = 0; i < length; i++)
		for (int j = i + 1; j < length; j++)
			if (sq[i] == sq[j])
				output_error("A rectangle must not intersect itself.");
	
	// Now glue the cubes together.
	for (int i = 0; i < length - 1; i++)
		glue_cubes(sq[i], sq[i+1], upright[i], upright[i+1]);
}

rectangle::rectangle(const rectangle &r)
{
	length = r.length;
	
	sq = new cube*[length];
	upright = new bool[length];
	
	for (int i = 0; i < length; i++)
	{
		sq[i] = r.sq[i];
		upright[i] = r.upright[i];
	}
}

rectangle::~rectangle()
{
	delete[] sq;
	delete[] upright;
}

// drill() 'drills' a tunnel down the core of this 
// rectangle's youngest child, connecting the 
// two boundary components (that this rectangle linked 
// between) together.
//
// Requires 12*length extra tetra and produce 4*length+8 
// additional free faces.
void rectangle::drill()
{
	output_debugging("drill", 2);
	
	manifold *M = sq[0]->home;
	
	// First find out if we need a new layer.
	bool need_new_layer = false;
	for (int i = 0; i < length; i++)
		if (! sq[i]->get_status())
		{
			need_new_layer = true;
			break;
		}
	
	if ((M->marked_points[front] == half_twisted) || (M->marked_points[back] == half_twisted))
		need_new_layer = true;
	
	if (need_new_layer) M->insert_layer();
	
	for (int i = 0; i < length; i++)
		sq[i]->set_status(false);
	
	M->marked_points[front] = drilled;
	M->marked_points[back] = drilled;
	
	// Convert each child cube to a gadget.
	std::vector<cube*> new_cubes;
	for (int i = 0; i < length; i++)
		new_cubes.push_back(convert_cube_to_gadget(sq[i], upright[i]));
	
	// Gadget pups can't get horizontal gluing information from their parents.
	// So we must provide it now.
	for (int i = 0; i < length - 1; i++)
	{
		int j = i + 1;
		// Glue pup_gadgets together.
		tetra *current_top_back = new_cubes[i]->get_tetra(low, left)->get_child();
		tetra *next_top_front = new_cubes[j]->get_tetra(top, left)->get_child();
		tetra *current_bottom_back = new_cubes[i]->get_tetra(low, right)->get_child();
		tetra *next_bottom_front = new_cubes[j]->get_tetra(top, right)->get_child();
		
		// (*1*) Gadgets were designed so this was true. :)
		current_top_back->gluesym(next_top_front, 0, perm(0,1,3,2));
		current_bottom_back->gluesym(next_bottom_front, 0, perm(0,1,3,2));
		
		glue_cubes(new_cubes[i], new_cubes[j], upright[i], upright[j]);
	}
	
	for (int i = 0; i < length; i++)
		delete new_cubes[i];
	
	return;
}

// We now define 3 types of walking move: swirl, tunnel_walk & boundary_walk.
// These will make it easier to attach the two-handle that realises a half-twist.
// All these walking moves obey the golden rule of walking:
// ALWAYS stay in the same layer.

// This table contains the permutations needed to
// subby in a new filler tetrahedron.
static const perm perm_walk_list[3] = {perm(3,0,2,1), perm(1,3,2,0), perm(2,1,3,0)};

// Two steps: Walk to a gadget_pup then to a NOT gadget_pup.
tetra *swirl(manifold *M, tetra *current, int layer_number)
{
	for (int i = 0; i < 3; i++)
	{
		if (current->get_gluedto(i) != NULL)
			continue;  // Stop us going backwards.
		
		current->walk_about(i);
		
		if ((current->get_gluedto(i))->get_layer() != layer_number) 
			current->ungluesym(i);  // Nope, don't want to go out of the layer.
		else if (current->get_gluedto(i)->get_category() == gadget_pup)
		{
			tetra *next_tetra = new tetra(M, filler);
			current->subbedby(next_tetra, i, perm_walk_list[i]);
			current->walk_about(i);
			current = next_tetra;
			break;
		}
		else
			current->ungluesym(i);  // Nope, don't want to go that way.
	}
	
	for (int i = 0; i < 3; i++)
	{
		if (current->get_gluedto(i) != NULL)
			continue;  // Stop us going backwards.
		
		current->walk_about(i);
		
		if (current->get_gluedto(i)->get_layer() != layer_number)
			current->ungluesym(i);  // Nope, don't want to go out of the layer.
		else if (current->get_gluedto(i)->get_category() != gadget_pup)
		{
			tetra *next_tetra = new tetra(M, filler);
			current->subbedby(next_tetra, i, perm_walk_list[i]);
			current->walk_about(i);
			current = next_tetra;
			break;
		}
		else
			current->ungluesym(i);  // Nope, don't want to go that way.
	}
	return current;
}

// One step: Walk to one that's the same category.
tetra *tunnel_walk(manifold *M, tetra *current, int layer_number)
{
	for (int i = 0; i < 3; i++)
	{
		if (current->get_gluedto(i) != NULL)
			continue;  // Stop us going backwards.
		
		current->walk_about(i);
		
		if (current->get_gluedto(i)->get_layer() != layer_number) 
			current->ungluesym(i);  // Nope, don't want to go out of the layer.
		else if (current->get_gluedto(i)->get_category() == current->get_gluedto(3)->get_category())
		{
			tetra *next_tetra = new tetra(M, filler);
			current->subbedby(next_tetra, i, perm_walk_list[i]);
			current->walk_about(i);
			current = next_tetra;
			break;
		}
		else
			current->ungluesym(i);  // Nope, don't want to go that way.
	}
	return current;
}

// One step: Walk to anything (in the layer) except a gadget_pup. Look out for markers.
tetra *boundary_walk(manifold *M, tetra *current, int layer_number)
{
	for (int i = 0; i < 3; i++)
	{
		if (current->get_gluedto(i) != NULL)
			continue;  // Stop us going backwards.
		
		current->walk_about(i);
		
		if (current->get_gluedto(i)->get_category() == marker)  // Good, we hit a marker we put in at the start.
		{
			return current->get_gluedto(i);  // We've already glued current to the marker, so just return it.
		}
		else if (current->get_gluedto(i)->get_layer() != layer_number)
			current->ungluesym(i);  // Nope, don't want to go out of the layer.
		else if (current->get_gluedto(i)->get_category() != gadget_pup)
		{
			tetra *next_tetra = new tetra(M, filler);
			current->subbedby(next_tetra, i, perm_walk_list[i]);
			current->walk_about(i);
			current = next_tetra;
			break;
		}
		else
			current->ungluesym(i);  // Nope, don't want to go that way.
	}
	return current;
}

// Let R denote this rectangle, connecting boundary components A & B
// not necessarily distinct. Let N be a neighborhood of R \cup A \cup B.
// Note that N is homeomorphic to a pair of pants.

// We will now discribe a half-twist on R. Recalling how we parse words 
// (see main in twister_main.cpp) a positive twist below corresponds to
// a left half-twist.

// Suppose that A and B are distinct. In this case the rectangle links
// together the 'cuffs' of the pants.

// I.e. If === is the rectangle joining boundary components A
// & B the half-twist has support on:

/*
//   /---------------\
//   | /-\       /-\ |
//   | |A|=======|B| |
//   | \-/       \-/ |
//   \---------------/
*/

// The second possibility is that A = B. In this case we say that the 
// rectangle is 'one ended'.

// I.e. If === is the rectangle from A to itself then the half-twist 
// has support on:

/*
//   /-------------\
//   | //=======\\ |
//   | ||       || |
//   | ||  /-\  || |
//   | ||  | |  || |
//   | ||  \-/  || |
//   | ||       || |
//   | ||  /-\  || |
//   | \\==|A|==// |
//   |     \-/     |
//   \-------------/
*/

// In the one ended case, a 'half-twist' on R is a twist of the same sign 
// on both components of the frontier of N.

void rectangle::half_twist(Sign whichway)
{
	output_debugging("half-twist", 2);
	
	// Let's remember some frequently used information.
	manifold *M = sq[0]->home;
	
	// First find out if we need a new layer.
	bool need_new_layer = false;
	for (int i = 0; i < length; i++)
		if (! sq[i]->get_status())
		{
			need_new_layer = true;
			break;
		}
	
	if ((M->marked_points[front] != standard) || (M->marked_points[back] != standard))
		need_new_layer = true;
	
	if (need_new_layer) M->insert_layer();
	
	drill();  // First drill a tunnel down the core of the rectangle.
	// Then attach a 2-handle in such a way as to perform a half-twist.
	
	for (int i = 0; i < length; i++)
		sq[i]->set_status(false);
	
	M->marked_points[front] = half_twisted;
	M->marked_points[back] = half_twisted;
	
	// We do this in 2 stages:
	// 1) Glue in a singular 2-handle.
	// 2) Thicken it to a regular 2-handle.
	
	// More frequently used information.
	int layer_number = M->get_num_layers();
	
	// Get some marker tetrahedra.
	tetra *first_marker = new tetra(M, marker);  // This will be our starting marker.
	tetra *mid_marker = new tetra(M, marker);  // This will be our midpoint marker.
	
	// Glue the marker in.
	first_marker->gluesym(sq[0]->get_tetra((upright[0]) ? 0 : 4)->get_child(), 3, (upright[0]) ? perm(1,3,0,2) : perm(1,3,0,2));
	
	if (one_ended)
	{
		mid_marker->gluesym(sq[0]->get_tetra((upright[0]) ? 2 : 5)->get_child()->get_child(), 3, (upright[0]) ? perm(0,1,3,2) : perm(0,1,3,2));
	}
	else
	{
		mid_marker->gluesym(sq[length-1]->get_tetra((upright[length-1]) ? 3 : 0)->get_child(), 3, (upright[length-1]) ? perm(0,3,2,1) : perm(3,1,2,0));
	}
	// ---===================================================---
	
	// Now begin working our way around, extending the singular 2-handle as we go.
	tetra *current = first_marker;  // Start at the first_marker.
	
	// Go through the tunnel for the first time (pre / post swirl as necessary).
	if (whichway == plus)  // Pre-swirl
		current = swirl(M, current, layer_number);  // See figure above convert_cube_to_gadget to see how this moves across the boundary faces.
	
	for (int c = 0; c < 2 * length + 2; c++)  // Work our way down the tunnel.
		current = tunnel_walk(M, current, layer_number);
	
	if (whichway == minus)  // Post-swirl
		current = swirl(M, current, layer_number);
	
	// Work our way around the 'far' boundary component until we get back to a marker.
	while (current != mid_marker)
		current = boundary_walk(M, current, layer_number);
	
	current->set_category(filler);  // Make our mid marker part of the singular 2-handle.
	
	// Go through the tunnel a second time (pre / post swirl as necessary).
	if ((whichway == minus) ^ (one_ended))  // Pre-swirl
		current = swirl(M, current, layer_number);
	
	for (int c = 0; c < 2 * length + 2; c++)  // Work our way through the tunnel again.
		current = tunnel_walk(M, current, layer_number);
	
	if ((whichway == plus) ^ (one_ended))  // Post-swirl
		current = swirl(M, current, layer_number);
	
	// Work our way around the 'near' boundary component until we hit our starting marker.
	while (current != first_marker)
		current = boundary_walk(M, current, layer_number);
	
	current->set_category(filler);  // Make our starting marker part of the singular 2-handle.
	// We've now gone all the way around once and added in a singular 2-handle.
	
	// We'd better go around again thickening the singular 2-handle into a regular one.
	// Note: Folding-off without thickening may produce a non-manifold.
	
	tetra *thicken_marker = new tetra(M, marker);  // This will be our thickening start point marker.
	// We glue thicken_marker down via its 3 face with 0 glued to the singularity of the 2-handle.
	// The edge of the 2-handle meeting the 1 & 2 vertices is 'top'.
	
	perm p = (current->get_gluedto(0) == NULL) ? perm(3,1,2,0) : ((current->get_gluedto(1) == NULL) ? perm(3,2,0,1) : perm(3,0,1,2));
	thicken_marker->gluesym(current, 3, p);
	current = thicken_marker;
	
	current->walk_about(1);
	while (current->get_gluedto(1) != thicken_marker) 
	{
		tetra *next_tetra = new tetra(M, filler);
		current->subbedby(next_tetra, 1, perm(0,3,1,2));
		current->walk_about(1);  // Re-close the faces.
		current = next_tetra;  // Move forwards
		current->walk_about(1);
	}
	
	current = current->get_gluedto(1);
	current->set_category(filler);
	
	// All done. Our singular 2-handle is now a regular one.
	
	return;
}

///////////// Manifold code //////////////////

manifold::manifold(std::string name_in, Manifold_type mytype)
{
	output_debugging("manifold", 2);
	
	name = name_in;
	manifold_type = mytype;
	num_layers = 0;
	num_cusps = 0;
	
	first_tetra = NULL;
	last_tetra = NULL;
}

manifold::~manifold()
{
	output_debugging("demanifold", 2);
	
	// Repeatedly remove the first tetrahedra until there's nothing left.
	while (first_tetra)
		oneless(first_tetra);
	
	for (int i = 0; i < int(cubes.size()); i++)
		delete cubes[i];
	
	return;
}

// Add a tetrahedra to the doubly linked list.
void manifold::onemore(tetra *newguy)
{
	if (last_tetra == NULL)  // If our manifold contains no tetra to begin with, use newguy to start the list.
	{
		first_tetra = newguy;
		last_tetra = newguy;
	}
	else
	{
		last_tetra->set_next(newguy);
		newguy->set_prev(last_tetra);
		last_tetra = newguy;
	}
	
	return;
}

// Remove a tetrahedra from the doubly linked list.
void manifold::oneless(tetra *oldguy)
{
	// 1) Unglue each face.
	for (int i = 0; i < 4; i++)
		oldguy->ungluesym(i);
	
	// Maintain the doubly linked list.
	// 2a) Go to prev and fix next, go to next and fix prev.
	// 2b) Fix first and last (if needed).
	if (oldguy->get_prev() != NULL) 
		oldguy->get_prev()->set_next(oldguy->get_next()); 
	else 
		first_tetra = oldguy->get_next();
	
	if (oldguy->get_next() != NULL) 
		oldguy->get_next()->set_prev(oldguy->get_prev());
	else 
		last_tetra = oldguy->get_prev();
	
	// 3) Delete the oldguy.
	delete oldguy;
	
	return;
}

// Create a cone over each component of the boundary of the 
// manifold and glue each cone on.
// Note: SnapPea will detect the topology of each vertex link.
tetra *manifold::capoff()
{
	output_debugging("capoff", 2);
	
	if (!last_tetra)  // Empty manifold?
		return NULL;
	
	tetra *old_last = last_tetra;
	
	// This describes how to attach newguy to the boundary face.
	perm perm_capoff_list[4] = {perm(3,1,2,0), perm(0,3,2,1), perm(0,1,3,2), perm(0,2,1,3)};
	
	// 1) Glue an extra tetrahedron onto each free face in the manifold.
	for (tetra *current = last_tetra; current != NULL; current = current->get_prev())  // Must go backwards, otherwise loop :(
		for (int i = 0; i < 4; i++)
		{
			if (current->get_gluedto(i) != NULL)
				continue;
			
			// Capping tetra have the same properties as the tetra they are connected to.
			tetra *newguy = new tetra(this, cap, current->get_position(), current->get_side(), current->get_layer());  
			newguy->gluesym(current, 3, perm_capoff_list[i]);
		}
	
	// 2) Glue the three free faces of the tetrahedra created in 1) together using walk_about().
	for (tetra *current = old_last->get_next(); current != NULL; current = current->get_next())
		for (int i = 0; i < 3; i++)
		{
			if (current->get_gluedto(i) != NULL)
				continue;
			
			current->walk_about(i);
		}
	
	return old_last->get_next();  // Return a pointer to the first cappoff tetra, if NULL then there's no boundary.
}

// Given vertices a and b of tetrahedron "basepoint" returns the degree of the corresponding edge.
int edge_degree(tetra *basepoint, int a, int b)
{
	// c and d are the indices of the other two vertices.
	int c = (a+1) % 4;
	if (c == b) c = (c+1) % 4;
	int d = 6 - (a + b + c);
	
	tetra *current = basepoint;
	int current_c = c;
	int current_d = d;
	int count = 0;
	do
	{
		perm p = current->get_gluing(current_c);
		current = current->get_gluedto(current_c);
		int old_d = current_d;
		current_d = p[current_c];
		current_c = p[old_d];
		count++;
	} while (!((current == basepoint) && (current_c == c) && (current_d == d)));
	
	return count;
}

// Installs a cusp number on the tetrahedra in the capoff
// Non-toroidal cusps are assigned a cusp number -1 so that
// canonical_cusps knows to avoid them later.
void manifold::identify_cusps(tetra *capoff_tetra)
{
	output_debugging("cusps", 2);
	
	// Clear the temp value everywhere.
	for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
		current->set_temp(-1);
	
	// Write theoretical cusp numbers to the temp field on every cusp.
	int theoretical_num_cusps = 0;  // This is how many cusps we find, although some might not be toroidal.
	for (tetra *basepoint = capoff_tetra; basepoint != NULL; basepoint = basepoint->get_next())
	{
		if (basepoint->get_temp() != -1) continue;  // Has already been assigned a theoretical cusp number.
		basepoint->set_temp(theoretical_num_cusps);
		
		std::queue<tetra *> tetra_list;
		tetra_list.push(basepoint);
		while (!(tetra_list.empty()))
		{
			tetra *current = tetra_list.front();
			tetra_list.pop();
			
			for (int i = 0; i < 3; i++)
			{
				tetra *adjacent = current->get_gluedto(i);
				if (adjacent->get_temp() == -1)
				{
					adjacent->set_temp(theoretical_num_cusps);
					tetra_list.push(adjacent);
				}
			}
		}
		theoretical_num_cusps++;
	}
	
	num_cusps = 0;  // The actual number of cusps that are toroidal. Klein bottles aren't possible as our manifolds are orientable.
	
	// Work on each cusp in turn.
	for (int cusp = 0; cusp < theoretical_num_cusps; cusp++)
	{
		std::vector<tetra *> cusp_tetra;  // All the tetra in this cusp.
		for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
			if (current->get_temp() == cusp) cusp_tetra.push_back(current);
		
		int num_faces = (int) cusp_tetra.size();
		
		// We've capped off, so the link of any point must be a closed surface.
		// Since Euler characteristic = V - E + F = V - (F/2), F should be even.
		if (num_faces % 2 != 0)
			output_error("Bad vertex link detected. Possible cause; not capping off before computing cusp information.");
		
		int num_vertices = num_faces * 3;
		std::vector<int> degree_frequency(num_vertices + 1, 0);  // The frequency of each degree.
		
		for (int i = 0; i < num_faces; i++)
			for (int j = 0; j < 3; j++)
				degree_frequency[edge_degree(cusp_tetra[i], j, 3)]++;
		
		if (degree_frequency[0] != 0)
			output_error("There should be no vertices of degree 0.");
		
		int Euler_characteristic = -(num_faces / 2);  // Now we work out the Euler characteristic of this cusp.
		
		for (int i = 1; i < num_vertices + 1; i++)
		{
			if (degree_frequency[i] % i != 0)
				output_error("Bad vertex link detected. A vertex of degree k should appear a multiple of k times.");
			
			Euler_characteristic += (degree_frequency[i] / i);
		}
		
		// Since the manifold is orientable, Euler characteristic must be even and <= 2.
		if ((Euler_characteristic > 2) || (Euler_characteristic % 2 != 0))
			output_error("Invalid cusp Euler characteristic.");
		
		// We could install the Euler characteristic to the cusp, but we wont need it again.
		int cusp_number = -1;  // The cusp number we're going to assign to this cusp.
		if (Euler_characteristic == 2)  // Spherical cusp.
		{
			if (this->get_manifold_type() == bundle) output_error("Spherical cusps are not possible in bundles.");
		}
		else if (Euler_characteristic == 0)  // Torus cusp.
		{
			cusp_number = num_cusps++;
		}
		else  // Other cusp.
		{
			output_warning("Non-torus Cusp.");  // Non-torus cusps. SnapPy wont like this.
		}
		
		// Write the real cusp number into the cusp tetra.
		for (int i = 0; i < num_faces; i++)
			cusp_tetra[i]->set_cusp_number(cusp_number);
	}
	
	return;
}

// Each cube has 12 triangular faces partitioned into 6 square faces.
// Suppose that T and T' make up a square face and capoff tetra A is 
// attached to T. We find the direction to move from A to B, the capoff
// tetra attached to T'.

// To understand this, see the diagram of the cube.
static const int dual_moves[6][4] = { {3,-1,3,-1},
									  {2,-1,-1,-1},
									  {-1,-1,0,-1},
									  {-1,3,-1,-1},
									  {-1,-1,3,-1},
									  {-1,0,0,-1} };

int dual_direction(tetra *A)
{
	tetra *T = A->get_gluedto(3);
	int T_index = 3 * (T->get_side()) + T->get_position();
	int T_face = A->get_gluing(3)[3];
	
	int move = dual_moves[T_index][T_face];
	if (move == -1) output_error("Face gluing error.");
	
	return T->get_gluing(T_face)[move];
}

// Install canonical meridians and longitudes on all torus cusps.

// L & M are orientated according to the convention in peripheral_curves.c (L:1299).
	// "The definition of the standard orientation for peripheral curves on
	// a torus is that when the fingers of your right hand point in the
	// direction of the meridian and your thumb points in the direction
	// of the longitude, the palm of your hand should face the cusp and
	// the back of your hand should face the fat part of the manifold."

void manifold::canonical_peripheral_curves(tetra *capoff_tetra)
{
	output_debugging("L&M", 2);
	
	if (!last_tetra)  // Empty manifold?
		return;
	
	if (capoff_tetra == NULL)  // We don't want to calculate this or there's no boundary faces.
		return;
	
	// Clear everything.
	for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
		current->set_temp(-1);
	
	// Start by assuming that we're building a bundle.
	// First we deal with the cusps coming from marked points of the surface.
	// Here the longitude is canonical (coming from the surface) and the 
	// dual 'meridian' slope is chosen to be shortest in the triangulation.
	// The longitude runs clockwise around the marked point (as seen on the surface
	// from above) while the meridian runs upwards (into your eye).
	for (tetra *basepoint = capoff_tetra; basepoint != NULL; basepoint = basepoint->get_next())
	{
		// Non-toroidal cusps have cusp number -1; avoid them.
		if ((basepoint->get_temp() != -1) || (basepoint->get_cusp_number() == -1)) continue;
		
		if (basepoint->get_layer() != num_layers) continue;
		
		// Work out the longitude.
		int bad_dir = 0;  // This will be the direction to move in to get out of this layer.
		for (bad_dir = 0; bad_dir < 3; bad_dir++)  // Work out the starting direction (actually the direction NOT to go in).
			if ((basepoint->get_gluedto(bad_dir))->get_layer() != num_layers) break;
		
		int back_dir = (bad_dir + 2) % 3;  // This will always hold the direction pointing back the way we just moved. Never go back this way.
		
		// Starting at the basepoint, install longitude until we get back to the basepoint.
		// NEVER exit the current layer.
		for (tetra *current = basepoint, *adjacent = NULL; adjacent != basepoint; current = adjacent)
		{
			current->set_temp(-2);
			for (int dir = 0; dir < 3; dir++)
			{
				if (dir == back_dir) continue;  // Don't move backwards.
				adjacent = current->get_gluedto(dir);
				if (adjacent->get_layer() != num_layers) continue;  // Don't move out of this layer.
				
				back_dir = (current->get_gluing(dir))[dir];  // Find out which way will take us backwards.
				current->set_peripheral_curves(longitude, dir, -1);  // Install the longitude in this direction.
				adjacent->set_peripheral_curves(longitude, back_dir, 1);
				break;
			}
		}
		
		// Now work out the meridian.
		// We look for a start and end on opposite sides of the longitude.
		int start_basepoint_dir = bad_dir, end_basepoint_dir = dual_direction(basepoint);
		
		tetra *start_basepoint = basepoint->get_gluedto(start_basepoint_dir);
		tetra *end_basepoint = basepoint->get_gluedto(end_basepoint_dir);
		
		// Start the meridian information.
		basepoint->set_peripheral_curves(meridian, start_basepoint_dir, -1);
		basepoint->set_peripheral_curves(meridian, end_basepoint_dir, 1);
		back_dir = (basepoint->get_gluing(end_basepoint_dir))[end_basepoint_dir];
		end_basepoint->set_peripheral_curves(meridian, back_dir, -1);
		back_dir = (basepoint->get_gluing(start_basepoint_dir))[start_basepoint_dir];
		start_basepoint->set_peripheral_curves(meridian, back_dir, 1);
		
		// Now find a path from start to end without crossing the longitude.
		// Remember, ones in the longitude have a temp of -2.
		
		// Flow distance forwards.
		int old_start_basepoint_temp = start_basepoint->get_temp();
		start_basepoint->set_temp(0);  // Seed flow at start_basepoint.
		int old_end_basepoint_temp = end_basepoint->get_temp();
		end_basepoint->set_temp(-1);  // Must allow the flow into end_basepoint.
		
		std::queue<tetra *> tetra_list;
		tetra_list.push(start_basepoint);
		while ((end_basepoint->get_temp() == -1) && !(tetra_list.empty()))
		{
			tetra *current = tetra_list.front();
			tetra_list.pop();
			for (int i = 0; i < 3; i++)
			{
				tetra *adjacent = current->get_gluedto(i);
				if (adjacent->get_temp() == -1)
				{
					adjacent->set_temp(current->get_temp() + 1);
					tetra_list.push(adjacent);
				}
			}
		}
		
		if (end_basepoint->get_temp() == -1)  // If the longitude was separating then clean up.
		{
			output_debugging("clean", 3);
			start_basepoint->set_temp(old_start_basepoint_temp);
			end_basepoint->set_temp(old_end_basepoint_temp);
			
			for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
			{
				if (current->get_cusp_number() != basepoint->get_cusp_number()) continue;
				if (current->get_temp() != -2) current->set_temp(-1);  // Leave -2's in as this longitude was separating.
				for (int i = 0; i < 4; i++)
				{
					current->set_peripheral_curves(longitude, i, 0);
					current->set_peripheral_curves(meridian, i, 0);
				}
			}
			continue;
		}
		
		// Flow meridian backwards.
		for (tetra *current = end_basepoint; current != start_basepoint; )
			for (int i = 0; i < 3; i++)
			{
				tetra *adjacent = current->get_gluedto(i);
				if (adjacent->get_temp() != (current->get_temp() - 1)) continue;
				
				current->set_peripheral_curves(meridian, i, 1);
				back_dir = (current->get_gluing(i))[i];
				adjacent->set_peripheral_curves(meridian, back_dir, -1);
				current = adjacent;
				break;
			}
		
		// Clear the temp data now that we're done.
		for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
			if (current->get_cusp_number() == basepoint->get_cusp_number())
				current->set_temp(-2);
		
		// If the manifold is a splitting we interchange the longitude and meridian. In order
		// to maintain the orientation convention, we flip the orientation of the meridian.
		if (this->get_manifold_type() == splitting)
			for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
				if (current->get_cusp_number() == basepoint->get_cusp_number())
					for (int i = 0; i < 3; i++)
					{
						int temp = current->get_peripheral_curves(longitude, i);
						current->set_peripheral_curves(longitude, i, current->get_peripheral_curves(meridian, i));
						current->set_peripheral_curves(meridian, i, -temp);
					}
	}
	
	// Now we repeat the whole process filling in the cusps that come from drillings.
	// Here the longitude is horizontal while the meridian goes up, around, down and back.
	for (tetra *basepoint = capoff_tetra; basepoint != NULL; basepoint = basepoint->get_next())
	{
		// Non-toroidal cusps have cusp number -1; avoid them.
		if ((basepoint->get_temp() != -1) || (basepoint->get_cusp_number() == -1)) continue;
		
		// We look for an even more specific type of basepoint.
		// To understand this, see the diagram of the cube.
		// 
		if (((basepoint->get_gluedto(3))->get_category() != gadget_left) || 
			((basepoint->get_gluedto(3))->get_position() != low) || 
			((basepoint->get_gluedto(3))->get_side() != right) ) continue;
		
		// Don't use basepoints that are glued to the end of rectangles (check for adjacent gadget_pups).
		bool good = true;
		for (int i = 0; i < 3; i++)
			if (((basepoint->get_gluedto(i))->get_gluedto(3))->get_category() == gadget_pup)
				good = false;
		
		if (!good) continue;
		
		// Now we work out the meridian starting from this basepoint.
		int back_dir = 0;
		bool upright = (basepoint->get_gluing(3)[3] == 2);  // Determine if the cube we attached to was a plus cube.
		int movements[4] = {3, upright ? 2 : 0, 1, 3};  // Array of movement directions.
		
		perm mount_to_basepoint = basepoint->get_gluing(3).inverse();
		tetra *current = basepoint;
		for (int i = 0; i < 4; i++)
		{
			int real_movement = (current->get_gluing(3).inverse())[movements[i]];
			current->set_temp(-2);
			current->set_peripheral_curves(meridian, real_movement, -1);
			back_dir = (current->get_gluing(real_movement))[real_movement];
			current = current->get_gluedto(real_movement);
			current->set_peripheral_curves(meridian, back_dir, 1);
		}
		if (current != basepoint) output_error("Misaligned meridian on a drilling.");
		
		// Now work out the longitude.
		// We look for a start and end on opposite sides of the meridian.
		
		int start_basepoint_dir = upright ? mount_to_basepoint[0] : mount_to_basepoint[2];
		int end_basepoint_dir = upright ? mount_to_basepoint[1] : mount_to_basepoint[0];
		
		tetra *start_basepoint = basepoint->get_gluedto(start_basepoint_dir);
		tetra *end_basepoint = basepoint->get_gluedto(end_basepoint_dir);
		
		// Start the longitude information.
		basepoint->set_peripheral_curves(longitude, start_basepoint_dir, -1);
		basepoint->set_peripheral_curves(longitude, end_basepoint_dir, 1);
		back_dir = (basepoint->get_gluing(end_basepoint_dir))[end_basepoint_dir];
		end_basepoint->set_peripheral_curves(longitude, back_dir, -1);
		back_dir = (basepoint->get_gluing(start_basepoint_dir))[start_basepoint_dir];
		start_basepoint->set_peripheral_curves(longitude, back_dir, 1);
		
		// Now find a path from start to end without crossing the meridian.
		// Remember, ones in the meridian have temp == -2.
		
		// Flow distance forwards.
		int old_start_basepoint_temp = start_basepoint->get_temp();
		start_basepoint->set_temp(0);  // Seed flow at start_basepoint.
		int old_end_basepoint_temp = end_basepoint->get_temp();
		end_basepoint->set_temp(-1);  // Better allow flow into end_basepoint.
		
		std::queue<tetra *> tetra_list;
		tetra_list.push(start_basepoint);
		while ((end_basepoint->get_temp() == -1) && !(tetra_list.empty()))
		{
			tetra *current = tetra_list.front();
			tetra_list.pop();
			for (int i = 0; i < 3; i++)
			{
				tetra *adjacent = current->get_gluedto(i);
				if (adjacent->get_temp() == -1)
				{
					adjacent->set_temp(current->get_temp() + 1);
					tetra_list.push(adjacent);
				}
			}
		}
		
		if (end_basepoint->get_temp() == -1)  // If the meridian was separating then clean up.
		{
			output_debugging("clean", 3);
			start_basepoint->set_temp(old_start_basepoint_temp);
			end_basepoint->set_temp(old_end_basepoint_temp);
			for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
			{
				if (current->get_cusp_number() != basepoint->get_cusp_number()) continue;
				if (current->get_temp() != -2) current->set_temp(-1);  // Leave -2's in as this meridian was separating.
				for (int i = 0; i < 4; i++)
				{
					current->set_peripheral_curves(longitude, i, 0);
					current->set_peripheral_curves(meridian, i, 0);
				}
			}
			continue;
		}
		
		// Flow longitude backwards.
		for (tetra *current = end_basepoint, *adjacent = NULL; current != start_basepoint; current = adjacent)
		{
			for (int i = 0; i < 3; i++)
			{
				adjacent = current->get_gluedto(i);
				if (adjacent->get_temp() != (current->get_temp() - 1)) continue;
				
				current->set_peripheral_curves(longitude, i, 1);
				back_dir = (current->get_gluing(i))[i];
				adjacent->set_peripheral_curves(longitude, back_dir, -1);
				break;
			}
		}
		
		// Clear the temp data now that we're done.
		for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
			if (current->get_cusp_number() == basepoint->get_cusp_number())
				current->set_temp(-2);
	}
	
	return;
}

// Ignoring for a moment, all the capoff tetrahedra,
// where it won't change the manifold (topologically) 
// glue adjacent boundary faces together. Doing so
// reduces the number of free faces, and thus the
// number of tetrahedra needed in capping off, making
// the manifold quicker to load in SnapPy.

// Warning: It is grossly wrong to try and perform 
// other actions on the manifold after doing this.
tetra *manifold::foldoff(tetra *capoff_tetra)
{
	output_debugging("foldoff", 2);
	
	if (!last_tetra)  // Empty manifold?
		return NULL;
	
	if (capoff_tetra == NULL)  // There's no boundary faces.
		return capoff_tetra;
	
	// Definition: Suppose that two boundary triangles T_1 = ABC and T_2 = A'B'C'
	// are glued along the edge BC ~ B'C' then T is said to have a valid fold in 
	// direction A iff A and A' are not in the same vertex class.
	
	// Definition: With T_1, T_2 as above in a manifold M, the manifold obtained 
	// by performing the 'close the book' move on T_1 in direction A is M / ~. 
	// Here ~ is the gluing of T_1 to T_2 by A ~ A', B ~ B' and C ~ C' extended 
	// linearly over T_1 and T_2.
	
	// Lemma: If T is a boundary triangle of the 3-manifold M which has a valid 
	// fold in direction A then the manifold obtained by performing the 'close the book' 
	// move on T in the direction A is homeomorphic to M.
	
	// Proof: See L122 of close_cusps.c in the SnapPea kernel. //
	
	// Lemma: If 2 boundary triangles each have no valid folds in any direction and,
	// after doing a valid fold elsewhere, they are now are adjacent to each other then 
	// they both still do not have a valid fold.
	
	// Proof: Omitted. //
	
	// Corollary: Only 1 pass through the list of tetra is ever needed.
	
	// We only foldoff capping tetra when it will preserve the meridian or longitude information.
	for (tetra *current = capoff_tetra; current != NULL; current = current->get_next())
	{
		// Let's find a direction in which we can foldoff.
		for (int dir = 0; dir < 3; dir++)  // Direction 3 is down into the manifold so don't bother trying that way.
		{
			tetra *adjacent = current->get_gluedto(dir);
			perm current_to_adjacent = current->get_gluing(dir);
			
			// Check we're not going to mess up the longitudes or meridians.
			bool no_peripheral_curves = true;  // Peripheral curve free?
			bool side_to_side = ((current->get_peripheral_curves(longitude, dir) != 0) && (current->get_peripheral_curves(meridian, dir) != 0));
			for (int i = 0; i < 2; i++)
			{
				Curves c = (Curves) i;  // Which type of curve we're looking at.
				for (int j = 0; j < 3; j++)
				{
					int k = current_to_adjacent[j];
					if ((current->get_peripheral_curves(c, j) != 0) || (adjacent->get_peripheral_curves(c, k) != 0)) no_peripheral_curves = false;
					if (current->get_peripheral_curves(c, j) != -adjacent->get_peripheral_curves(c, k)) side_to_side = false;
				}
			}
			
			if (!no_peripheral_curves && !side_to_side)
				continue;
			
			if (adjacent == current)
				continue;  // This would be caught by the A == A' test later, but it's an easy catch here.
			
			int A = dir;  // These are the 'opposite vertices'.
			int A_prime = current_to_adjacent[A];
			
			// Now we check A != A'. If this is the case then we have sufficent grounds to foldoff this way.
			tetra *walking = current;
			
			// Let's label the vertices of the current triangle.
			int curr_a = A, curr_b = (curr_a == 0) ? 1 : 0, curr_c = 3 - curr_a - curr_b;
			
			// Walk around the edge between A and vertex 3.
			do
			{
				perm current_gluing = walking->get_gluing(curr_b);
				
				walking = walking->get_gluedto(curr_b);
				
				curr_a = current_gluing[curr_a];
				curr_b = current_gluing[curr_b];
				curr_c = current_gluing[curr_c];
				
				curr_b ^= curr_c; curr_c ^= curr_b; curr_b ^= curr_c;  // Interchange curr_b & curr_c.
			} while (((walking != current) || (curr_a != A)) && ((walking != adjacent) || (curr_a != A_prime)));
			// Stop only when we get back to current or adjacent with the CORRECT orientation.
			
			if (walking == adjacent)  // A == A'. So give up and try a different direction.
				continue;
			
			// A != A'. All's good, we definitely can foldoff this way.
			
			// Better record the previous one so we can take a step forward once we delete current.
			tetra *prev = (current->get_prev() != adjacent) ? current->get_prev() : adjacent->get_prev();
			
			// Unglue pairs of corresponding faces (one pair from current & one pair from adjacent)
			// and glue them together to allow us to remove current and adjacent from the manifold,
			// reducing the number of free faces by 2.
			// Note: This process maintains the cappoff, so no recapping is required.
			for (int i = 0; i < 4; i++)
			{
				if (i == dir)
					continue;  // Now it would be a waste of time to even think about gluing that way.
				
				tetra *from = current->get_gluedto(i);
				tetra *to = adjacent->get_gluedto(current_to_adjacent[i]);
				int glueout = current->get_gluing(i)[i];
				
				perm from_to_current = (current->get_gluing(i)).inverse();
				perm adjacent_to_to = adjacent->get_gluing(current_to_adjacent[i]);
				perm composed = adjacent_to_to.of(current_to_adjacent.of(from_to_current));
				
				adjacent->ungluesym(current_to_adjacent[i]);
				current->ungluesym(i);
				
				from->gluesym(to, glueout, composed);
			}
			
			// Delete current & adjacent as they are now disconnected from the manifold.
			oneless(adjacent);
			oneless(current);
			
			current = prev;  // Take a step back.
			break;
		}
	}
	
	// The capoff has changed, so better re-find the capoff_tetra.
	for (capoff_tetra = first_tetra; capoff_tetra != NULL; capoff_tetra = capoff_tetra->get_next())
		if (capoff_tetra->get_category() == cap)
			break;
	
	return capoff_tetra;
}

// Duplicate layer #0 and insert the copy (layer #L) between layer #0 and layer #(L-1).
// Order: Bottom -> 0, L, (L-1), (L-2), ..., 2, 1 <- Top.
void manifold::insert_layer()
{
	output_debugging("insert_layer", 2);
	
	if (!last_tetra)  // Empty manifold?
		return;
	
	// Reset the layer flags on cubes and marked points.
	for (int i = 0; i < int(marked_points.size()); i++)
		marked_points[i] = standard;
	
	for (int i = 0; i < int(cubes.size()); i++)
		cubes[i]->set_status(true);
	
	num_layers++;
	
	// Duplicate layer #0.
	for (int i = 0; i < int(cubes.size()); i++)
	{
		cube *C = cubes[i];
		cube *C2 = new cube(*this, num_layers, C);
		
		C->get_tetra(top, left)->subbedby(C2->get_tetra(top, left), 1, perm(0,1,2,3));
		C->get_tetra(top, right)->subbedby(C2->get_tetra(top, right), 2, perm(0,1,2,3));
		
		C->get_tetra(top, left)->gluesym(C2->get_tetra(low, left), 1, perm(0,3,2,1));
		C->get_tetra(top, right)->gluesym(C2->get_tetra(low, right), 2, perm(0,1,3,2));
		
		delete C2;
	}
	
	// Glue the tetra in layer #L together and subbedby the top most faces in.
	for (tetra *current = last_tetra; current->get_layer() == num_layers; current = current->get_prev())
		for (int i = 0; i < 4; i++)
			if (current->get_gluedto(i) == NULL)
			{
				tetra *myparent = current->get_parent();
				tetra *myuncle = myparent->get_gluedto(i);
				
				if ((myuncle != NULL) && (myuncle->get_layer() == 0))
				{
					tetra *mycousin = myuncle->get_child();
					if (mycousin != NULL) 
						current->gluesym(mycousin, i, myparent->get_gluing(i));
				}
			}
	
	return;
}

void manifold::tidy_boundary()
{
	if (GLOBAL_calculate_peripheral_curves) insert_layer();  // We'll need an extra layer to find the longitudes in.
	
	tetra *capoff_tetra = capoff();  // Deal with boundary faces with a mandatory capoff().
	identify_cusps(capoff_tetra);  // Compute the cusps and their number.
	if (GLOBAL_calculate_peripheral_curves) canonical_peripheral_curves(capoff_tetra);
	if (GLOBAL_optimise) capoff_tetra = foldoff(capoff_tetra);  // If allowed to optimise, reduce the number of tetrahedra by folding off.
	return;
}

std::string manifold::to_string()
{
	std::stringstream ss;
	snap_print(ss);
	return ss.str();
}

///////////// Printing code //////////////////

std::ostream &operator<<(std::ostream &o, const perm &to_print)
{
	for (int i = 0; i < 4; i++)
		o << to_print.image[i];
	
	return o;
}

void manifold::snap_print(std::ostream &o)
{
	output_debugging("print", 2);
	
	// Watch out - blank names will kill snappea.
	if (name == "")
		name = "Blank name disallowed";
	
	// Print basic information.
	o << "% Triangulation" << std::endl;
	o << name << std::endl;
	o << "not_attempted  0.00000000" << std::endl;
	o << "oriented_manifold" << std::endl;  // All manifolds we build are orientable, as all face gluings are odd. 
	o << "CS_unknown" << std::endl;
	o << std::endl;
	
	// As our manifold is orientable we must only have torus cusps, so we write them down.
	o << num_cusps << " 0" << std::endl;
	for (int i = 0; i < num_cusps; i++)
		o << "    torus  0.00000000  0.00000000" << std::endl;
	
	o << std::endl;
	
	// First assign each tetra a snap_index
	int i = 0; 
	for (tetra *current = first_tetra; current != NULL; current = current->get_next())
		current->snap_index = i++;
	
	// Print the number of tetrahedra.
	o << i << std::endl;
	
	// Finally, print the face gluings for each tetrahedron and the cusp data, if any.
	for (tetra *current = first_tetra; current != NULL; current = current->get_next())
		current->print_wrt(o);
	
	o << std::endl;
	return;
}

void tetra::print_wrt(std::ostream &o)
{
	for (int i = 0; i < 4; i++)
		o << "   " << gluedto[i]->snap_index;
	
	o << std::endl;
	
	for (int i = 0; i < 4; i++)
		o << " " << gluing[i];
	
	o << std::endl;
	
	o << "  -1   -1   -1   " << cusp_number << std::endl;  // Material vertices 0 through 2 and possibly a ideal vertex at 3.
	o << "  0  0  0  0  0  0  0  0  0  0  0  0  " << peripheral_curves[meridian][0] << "  " << peripheral_curves[meridian][1] << "  " << peripheral_curves[meridian][2] << "  " << peripheral_curves[meridian][3] << std::endl;  //  meridian (right sheet)
	o << "  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << std::endl;  //  meridian (left sheet)
	o << "  0  0  0  0  0  0  0  0  0  0  0  0  " << peripheral_curves[longitude][0] << "  " << peripheral_curves[longitude][1] << "  " << peripheral_curves[longitude][2] << "  " << peripheral_curves[longitude][3] << std::endl;  //  longitude (right sheet)
	o << "  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0" << std::endl;  // longitude (left sheet)
	o << "  0.000000000000   0.000000000000" << std::endl;
	o << std::endl;
}


