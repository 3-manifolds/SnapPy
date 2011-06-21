#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "global.h"
#include "parsing.h"
#include "twister.h"

// States for the FSM.
enum State
{
	twisting, drilling
};

int main(int argc, char **argv)
{
	// Set default values.
	Manifold_type manifold_type = bundle;
	std::string manifold_name = "", handles = "", gluing = "", surface_file = "", output_file = "";
	
	parse_input(argc, argv, surface_file, output_file, manifold_name, manifold_type, gluing, handles);
	
	// Build a (blank) surface with the correct name.
	manifold M(manifold_name, manifold_type);
	
	// These will hold all the pieces of information for the manifold.
	std::vector<square> squares;
	std::vector<annulus> annuli;
	std::vector<rectangle> rectangles;
	std::vector<std::string> annulus_name;
	std::vector<std::string> annulus_name_inverse;
	std::vector<std::string> rectangle_name;
	std::vector<std::string> rectangle_name_inverse;
	
	// Fill up the surface with interesting information from the surface_file.
	build_manifold(M, squares, annuli, rectangles, annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, surface_file);
	
	if (handles != "")  // Attach handles (if needed).
	{
		size_t marker = 0;
		
		do
		{
			std::string handle_location = find_next_substring(handles, seperator, marker);
			
			bool performed_action = false;  // Records if we actually performed an action.
			for (int i = 0; i < int(annulus_name.size()); i++)
			{
				if (annulus_name[i] == handle_location) 
				{
					annuli[i].twohandle(above);
					performed_action = true;
					break;
				}
				
				if (annulus_name_inverse[i] == handle_location)
				{
					annuli[i].twohandle(below);
					performed_action = true;
					break;
				}
			}
			if (!performed_action) output_warning("No annulus named '" + handle_location + "' found.");
		} while (marker != 0);
	}
	
	M.insert_layer();  // Insert layer #1 to prevent rooftop access.
	
	State current_state = twisting;  // Starting mode.
	
	// We parse the gluing w = w_1 ... w_n from left to right, so our twists are stacked
	// from top to bottom. So when a curve c at the bottom is pushed up 
	// the resulting curve at the top is w(c) = w_1( ... w_n(c) ... ), i.e. mapping
	// classes act on the left.
	if (gluing != "")  // Twist / drill (if needed).
	{
		size_t marker = 0;
		
		do
		{
			// Work out what state we will be working in.
			if (marker != 0)
			{
				if (gluing.substr(marker - 1, 1) == start_drill)
					current_state = drilling;
				else if (gluing.substr(marker - 1, 1) == start_twist)
					current_state = twisting;
			}
			
			std::string action_name = find_next_substring(gluing, product, marker);
			
			if (action_name == "") continue;  // Skip it if there is no action.
			
			bool performed_action = false;  // Records if we actually performed an action.
			for (int i = 0; i < int(annulus_name.size()); i++)
			{
				if ((annulus_name[i] != action_name) && (annulus_name_inverse[i] != action_name)) continue;
				if (current_state == twisting)
					annuli[i].twist((annulus_name[i] == action_name) ? plus : minus);
				else
					annuli[i].drill();
				
				performed_action = true;
				break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
			}
			if (performed_action) continue;
			
			for (int i = 0; i < int(rectangle_name.size()); i++)
			{
				if ((rectangle_name[i] != action_name) && (rectangle_name_inverse[i] != action_name)) continue;
				if (current_state == twisting)
					rectangles[i].half_twist((rectangle_name[i] == action_name) ? plus : minus);
				else
					rectangles[i].drill();
				
				performed_action = true;
				break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
			}
			if (performed_action) continue;
			
			// If our action_name didn't match any curves then throw a warning.
			output_warning("No annulus or rectangle named '" + action_name + "' found.");
		} while (marker != 0);
	}
	
	M.tidy_boundary();
	
	// Now output this information.
	if (output_file != "")
	{
		// If we were given a destination file, write it there.
		std::ofstream myfile(output_file.c_str());
		if (not myfile.is_open())
			output_error("Unable to write to " + output_file + ".");
		
		M.snap_print(myfile);
		myfile.close();
	}
	else
	{
		// Else write it to std::cout.
		M.snap_print(std::cout);
	}
	
	return 0;
}
