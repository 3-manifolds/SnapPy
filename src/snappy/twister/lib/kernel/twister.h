#ifndef TWISTER_HEADER_H
#define TWISTER_HEADER_H

#include <iostream>
#include <vector>
#include "global.h"

// The different types of manifolds that we can build.
enum Manifold_type
{
	splitting, bundle
};

// Provides an argument for annulus::twist() and rectangle::half_twist().
enum Sign
{
	plus = 1,
	minus = -1
};

// Tells the tetrahedron where in the manifold it lives.
enum Category
{
	cubie, speedbump, handle, cap, gadget_left, gadget_right, gadget_pup, filler, marker, none
};

enum Position
{
	top, mid, low, other
};

enum Side
{
	left, right, neither
};

enum Curves
{
	longitude, meridian
};

enum Marked_status
{
	standard, drilled, half_twisted
};

// Permutations (low-level; used by other classes).
class perm
{
	int image[4];
	
public:
	perm();
	
	perm(const perm &to_copy);
	
	// Image of 0-3 under the permutation.
	perm(const int* of_in);
	
	// Constructor:  perm(0,1,2,3) is the identity.
	perm(int im0, int im1, int im2, int im3)
	{
		image[0] = im0;
		image[1] = im1;
		image[2] = im2;
		image[3] = im3;
		
		// Make sure the thing that we've just made is a permutation.
		// As users have no direct control over permutations, this is to check hard coded permutations are valid.
		for (int i=0; i<4; i++)
			for (int j=i+1; j<4; j++)
				if (image[i] == image[j]) output_error("Invalid permutation.");
	}
	
	// Inverse.
	perm inverse() const;
	
	// Overloaded subscript operator.
	int operator[](int input) const;
	
	// Composition of permutations.
	perm of(const perm &other) const;
	
	// Allows printing.
	friend std::ostream &operator<<(std::ostream &o, const perm &to_print);
};

// Tetrahedron  (low-level; used by other classes)
class tetra
{
	tetra *next;  // For the doubly linked list
	tetra *prev;
	
	tetra *gluedto[4];  // Gluing information
	perm gluing[4];
	
	Category category;  // Location
	Position position;
	Side side;
	
	tetra *parent;  // For inserting layers
	tetra *child;
	
	int layer;  // Which layer it belongs to
	
	int peripheral_curves[2][4];
	int cusp_number;
	int temp;
	
public:
	tetra(class manifold *M, Category mycategory = none, Position myposition = other, Side myside = neither, int mylayer = -1, tetra *my_parent = NULL);
	~tetra();
	
	int snap_index;
	
	void print_wrt(std::ostream &o);
	
	// These let us get / set private properties.
	tetra *get_next(){return next;}
	void set_next(tetra *target){next = target;}
	
	tetra *get_prev(){return prev;}
	void set_prev(tetra *target){prev = target;}
	
	tetra *get_child(){return child;}
	void set_child(tetra *target){child = target; target->parent = this;}
	
	tetra *get_parent(){return parent;}
	void set_parent(tetra *target){parent = target; target->child = this;}
	
	tetra *get_gluedto(int whichface){return gluedto[whichface];}
	void set_gluedto(int whichface, tetra *target){gluedto[whichface] = target;}
	
	perm get_gluing(int whichface){return gluing[whichface];}
	void set_gluing(int whichface, perm how)
	{
		// Check how is an odd permutation. This can only arrise from
		// programmer error as users have no direct control over permutation.
		bool c = true;  // c == 'how is even'.
		for (int i = 1; i < 4; i++)
			for (int j = 0; j < i; j++)
				if (how[j] > how[i])
					c ^= true;
		
		if (c) output_error("Even gluing detected.");
		
		gluing[whichface] = how;
	}
	
	Category get_category(){return category;}
	void set_category(Category new_category){category = new_category;}
	
	int get_peripheral_curves(Curves c,int whichface){return peripheral_curves[c][whichface];}
	void set_peripheral_curves(Curves c, int whichface, int target){peripheral_curves[c][whichface] = target;}
	
	int get_cusp_number(){return cusp_number;}
	void set_cusp_number(int target){cusp_number = target;}
	
	int get_temp(){return temp;}
	void set_temp(int target){temp = target;}
	
	// There are no set functions for these properties as they cannot be changed after creation.
	
	Position get_position(){return position;}
	Side get_side(){return side;}
	int get_layer(){return layer;}
	
	void walk_about(int fromface);
	
	// gluesym may only be used to glue a pair of free faces.
	// 'how' is, in our applications, always odd --- i.e. our
	// manifolds are oriented, with all tetrahedra right-handed. 
	// See Jeff's comments in kernel_typedefs.h line 154.
	void gluesym(tetra *whereglue, int whichface, const perm &how);
	void ungluesym(int whichface);
	
	// whichface must be glued already to some face. 
	// 'how' should be an even permutation.
	// See the function body for usage.
	void subbedby(tetra *whereglue, int whichface, const perm &how);
};

// An intersection point, thickened in both directions
// For examples of use, see the file sample.c
class cube 
{
	tetra *topleft, *midleft, *lowleft;
	tetra *topright, *midright, *lowright;
	
	bool is_glued[2];  // Has it been glued in each direction?
	bool status;
	
public:
	class manifold *home;
	// A cube needs to know what manifold it belongs to.
	cube(class manifold &home_in, const int layer_number = 0, cube *parent_cube = NULL);
	
	// We need to be able to request the specific tetra that make up a cube.
	tetra *get_tetra(int i)  // This for looping.
	{
		switch (i)
		{
			case 0:
				return topleft;
			case 1:
				return midleft;
			case 2:
				return lowleft;
			case 3:
				return topright;
			case 4:
				return midright;
			case 5:
				return lowright;
		}
		return NULL;
	}
	
	tetra *get_tetra(Position p, Side s)  // This is for requesting by name
	{
		if (s == left)
		{
			if (p == top)
				return topleft;
			else if (p == mid)
				return midleft;
			else if (p == low)
				return lowleft;
			else
				output_error("Invalid position requested.");
		}
		else if (s == right)
		{
			if (p == top)
				return topright;
			else if (p == mid)
				return midright;
			else if (p == low)
				return lowright;
			else
				output_error("Invalid position requested.");
		}
		else
			output_error("Invalid side requested.");
		
		return NULL;
	}
	
	bool get_status(){return status;};
	void set_status(bool target){status = target;};
	
	void closeup_cube()
	{
		// Glue the bottoms to the tops.
		lowleft->gluesym(topleft, 3, perm(0,3,2,1));
		lowright->gluesym(topright, 3, perm(0,1,3,2));
	};
	
	friend class annulus;
	friend class rectangle;
};

// A thickened closed curve, i.e. a thickened annulus.
// For examples of use, see the file sample.c
class annulus
{
	int length;
	cube **sq;
	bool *upright;  // true for each + intersection
	
public:
	annulus(const std::vector<cube*> &sq_in, const std::vector<bool> &upright_in);
	annulus(const annulus &a);
	
	~annulus();
	
	void drill();  // Removes a core of the annulus.
	
	void twohandle(bool is_above);
	
	// Plumbing on a curve on a surface in a manifold,
	// or cutting regluing after a Dehn twist
	// (positive/negative "speed bump" glued above)
	void twist(Sign whichway);
	
};

// A thickened arc
// See the file sample.c for an example of use
class rectangle 
{
	int length;
	cube **sq;
	bool *upright;  // true for each + intersection
	int front;
	int back;
	bool one_ended;  // True if and only if both ends of the rectangle lie in the same vertical boundary component of the thickened surface.
	
public:
	// Same construction options as for annulus
	rectangle(const std::vector<cube*> &sq_in, const std::vector<bool> &upright_in);
	rectangle(const rectangle &r);
	~rectangle();
	
	void drill();
	void half_twist(Sign whichway);
	
	int get_front(){return front;}
	void set_front(int target){front = target;}
	int get_back(){return back;}
	void set_back(int target){back = target;}
	bool get_one_ended(){return one_ended;}
	void set_one_ended(bool target){one_ended = target;}
};

// A three-manifold, including a list of its tetrahedra
class manifold
{
	// Pointers to the linked lists.
	// One for the tetrahedra:
	tetra *first_tetra;
	tetra *last_tetra;
	
	Manifold_type manifold_type;
	
	std::string name;
	int num_layers;
	int num_cusps;
	void onemore(tetra *newguy);
	void oneless(tetra *oldguy);
	
	tetra *capoff();
	void identify_cusps(tetra *capoff_tetra);
	void canonical_peripheral_curves(tetra *capoff_tetra);
	tetra *foldoff(tetra *capoff_tetra);
	
public:
	manifold(std::string name_in, Manifold_type mytype);
	~manifold();
	
	void insert_layer();
	void tidy_boundary();
	
	void snap_print(std::ostream &o);
	std::string to_string();
	
	int get_num_layers(){return num_layers;}
	
	Manifold_type get_manifold_type(){return manifold_type;}
	
	std::vector<cube*> cubes;
	std::vector<Marked_status> marked_points;
	
	friend tetra::tetra(class manifold *M, Category mycategory, Position myposition, Side myside, int mylayer, tetra *my_parent);
};

#endif
