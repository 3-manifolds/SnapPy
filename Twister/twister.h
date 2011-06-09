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

// Provides an argument for annulus::twohandle().
enum Direction
{
	above = true,
	below = false
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
	cube, speedbump, handle, cap, gadget_left, gadget_right, gadget_pup, filler, marker, none
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
	tetra(class manifold *M, Category mycategory = none, Position myposition = other, Side myside = neither, int mylayer = -1);
	~tetra();
	
	int snap_index;
	
	void print_wrt(std::ostream &o);
	
	// These let us get / set private properties.
	tetra *get_next()
	{
		return next;
	}
	
	void set_next(tetra *target)
	{
		next = target;
	}
	
	tetra *get_prev()
	{
		return prev;
	}
	
	void set_prev(tetra *target)
	{
		prev = target;
	}
	
	tetra *get_child()
	{
		return child;
	}
	
	void set_child(tetra *target)
	{
		child = target;
	}
	
	tetra *get_parent()
	{
		return parent;
	}
	
	void set_parent(tetra *target)
	{
		parent = target;
	}
	
	tetra *get_gluedto(int whichface)
	{
		return gluedto[whichface];
	}
	
	void set_gluedto(int whichface, tetra *target)
	{
		gluedto[whichface] = target;
	}
	
	perm get_gluing(int whichface)
	{
		return gluing[whichface];
	}
	
	void set_gluing(int whichface, perm how)
	{
		// Check how is a odd permutation
		
		bool c = true;  // c == 'how is even'.
		bool b[4] = {false, false, false, false};
		for (int i = 0; i < 4; i++)
		{
			if (b[i]) continue;
			c ^= true;
			for (int j = i; not b[j]; j = how[j])
			{
				b[j] = true;
				c ^= true;
			}
		}
		
		if (c) output_error("Even gluing detected.");
		
		gluing[whichface] = how;
	}
	
	Category get_category()
	{
		return category;
	}
	
	void set_category(Category new_category)
	{
		category = new_category;
	}
	
	int get_peripheral_curves(Curves c,int whichface)
	{
		return peripheral_curves[c][whichface];
	}
	
	void set_peripheral_curves(Curves c, int whichface, int target)
	{
		peripheral_curves[c][whichface] = target;
	}
	
	int get_cusp_number()
	{
		return cusp_number;
	}
	
	void set_cusp_number(int target)
	{
		cusp_number = target;
	}
	
	int get_temp()
	{
		return temp;
	}
	
	void set_temp(int target)
	{
		temp = target;
	}
	
	// There are no set functions for these properties as they cannot be changed after creation.
	
	Position get_position()
	{
		return position;
	}
	
	Side get_side()
	{
		return side;
	}
	
	int get_layer()
	{
		return layer;
	}
	
	void walk_about(int fromface);
	
	// gluesym may only be used to glue a pair of free faces.
	// 'how' is, in our applications, always odd --- i.e. our
	// manifolds are oriented, with all tetrahedra right-handed. 
	void gluesym(tetra *whereglue, int whichface, const perm &how);
	void ungluesym(int whichface);
	
	// whichface must be glued already to some face. 
	// 'how' should be an even permutation.
	// See the function body for usage.
	void subbedby(tetra *whereglue, int whichface, const perm &how);
};

// Oriented square
// Provides an input format for annulus and rectangle constructors
struct orisquare
{
	class square *sq;
	bool is_upright;
	orisquare *next;
	int depth;
	
	orisquare(Sign whichway, class square *mysquare)
	{
		sq = mysquare;
		is_upright = (whichway == plus);
		next = NULL;
		depth = 1;
	}
	
	orisquare &operator<<(orisquare &other)
	{
		orisquare *iter = this;
		
		while (iter->next)
		{
			iter->depth += other.depth;
			iter = iter->next;
		}
		iter->depth += other.depth;
		
		iter->next = &other;
		
		return *this;
	}
};

// An intersection point, thickened in both directions
// For examples of use, see the file sample.c
class square 
{
	tetra *topleft, *midleft, *lowleft;
	tetra *topright, *midright, *lowright;
	
	bool is_glued[2];  // Has it been glued in each direction?
	
public:
	class manifold *home;
	// A square needs to know what manifold it belongs to.
	square(class manifold &home_in);
	
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
	
	// The constructor for annuli and rectangles requires arguments
	// of the form +square, -square.
	orisquare &operator+();
	orisquare &operator-();
	
	friend class annulus;
	friend class rectangle;
};

// A thickened closed curve, i.e. an annulus.
// For examples of use, see the file sample.c
class annulus
{
	int length;
	square **sq;
	bool *upright;  // true for each + intersection
	
public:
	annulus(const std::vector<square*>&, const std::vector<bool>&);
	
	// The old way of constructing an annulus: give it an
	// array of pointers to squares, and an array of booleans
	// indicating orientation (true = '+').
	annulus(int length_in, square **sq_in, bool *upright_in);
	annulus(orisquare &ori);
	
	annulus(const annulus &a);
	
	~annulus();
	
	void drill();  // Removes a core of the annulus.
	
	void twohandle(Direction is_above);
	
	// For backwards compatibility
	void twohandleabove()
	{
		twohandle(above);
		return;
	}
	
	void twohandlebelow()
	{
		twohandle(below);
		return;
	}
	
	// Plumbing on a curve on a surface in a manifold,
	// or cutting regluing after a Dehn twist
	// (positive/negative "speed bump" glued above)
	void twist(Sign whichway);
	
	// alternate calling method for twist()
	void operator++();  
	void operator--();  
};

// A thickened arc
// See the file sample.c for an example of use
class rectangle 
{
	int length;
	square **sq;
	bool *upright;  // true for each + intersection
	
public:
	// Same construction options as for annulus
	rectangle(const std::vector<square*>&, const std::vector<bool>&);
	rectangle(orisquare &ori);
	rectangle(int length_in, square **sq_in, bool *upright_in);
	
	rectangle(const rectangle &r);
	~rectangle();
	
	void drill();
	void half_twist(Sign whichway);
};

// A three-manifold, including a list of its tetrahedra
class manifold
{
	// Pointers to the linked lists.
	// One for the tetrahedra:
	tetra *first_tetra;
	tetra *last_tetra;
	
	Manifold_type manifold_type;
	
	char *name;
	int num_layers;
	int num_cusps;
	void onemore(tetra *newguy);
	
	tetra *capoff();
	void identify_cusps(tetra *capoff_tetra);
	void canonical_peripheral_curves(tetra *capoff_tetra);
	tetra *foldoff(tetra *capoff_tetra);
	
public:
	manifold(char *name_in, Manifold_type mytype);
	~manifold();
	
	void oneless(tetra *oldguy);
	
	void insert_layer();
	void tidy_boundary();
	
	void snap_print(std::ostream &o);
	
	int get_num_layers()
	{
		return num_layers;
	}
	
	Manifold_type get_manifold_type()
	{
		return manifold_type;
	}
	
	friend tetra::tetra(class manifold *M, Category mycategory = none, Position myposition = other, Side myside = neither, int mylayer = -1);
};

#endif
