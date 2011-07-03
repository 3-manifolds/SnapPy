#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "global.h"
#include "parsing.h"
#include "twister.h"

int main(int argc, char **argv)
{
	// Set default values.
	set_globals_to_defaults();
	Manifold_type manifold_type = bundle;
	std::string manifold_name = "", handles = "", gluing = "", surface_file = "", output_file = "";
	
	// Parse the given command.
	parse_input(argc, argv, surface_file, output_file, manifold_name, manifold_type, gluing, handles);
	
	// Build a (blank) surface with the correct name.
	manifold M(manifold_name, manifold_type);
	
	// These will hold all the pieces of information for the manifold.
	std::vector<square> squares;
	std::vector<annulus> annuli;
	std::vector<rectangle> rectangles;
	std::vector<std::string> annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse;
	
	// Fill up the surface with interesting information from the surface_file.
	build_manifold(M, squares, annuli, rectangles, annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, surface_file);
	
	int num_annuli = int(annulus_name.size());
	int num_rectangles = int(rectangle_name.size());
	
	if (handles != "")  // Attach handles (if needed).
	{
		size_t marker = 0;
		
		do
		{
			std::string handle_location = find_next_substring(handles, seperator, marker);
			
			bool performed_action = false;  // Records if we actually performed an action.
			for (int i = 0; i < num_annuli; i++)
			{
				if ((annulus_name[i] == handle_location) || (annulus_name_inverse[i] == handle_location))
				{
					annuli[i].twohandle((annulus_name[i] == handle_location) ? above : below);
					performed_action = true;
					break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
				}
			}
			if (performed_action) continue;
			
			output_warning("No annulus named '" + handle_location + "' found to attach a handle to.");
		} while (marker != 0);
	}
	
	M.insert_layer();  // Insert layer #1 to prevent rooftop access.
	
	// We parse the gluing w = w_1 ... w_n from left to right, so our twists are 
	// stacked from top to bottom. So when a curve c at the bottom is pushed up 
	// the resulting curve at the top is w(c) = w_1( ... w_n(c) ... ), i.e. mapping
	// classes act on the left.
	if (gluing != "")  // Twist / drill (if needed).
	{
		size_t marker = 0;
		
		do
		{
			std::string action_name = find_next_substring(gluing, seperator, marker);
			
			if (action_name == "") continue;  // Skip it if there is no action.
			
			bool drilling = (action_name.substr(0, 1) == drill);  // Are we drilling?
			if (drilling) action_name = action_name.substr(1);  // If so, drop the drill symbol.
			
			bool performed_action = false;  // Records if we actually performed an action.
			for (int i = 0; i < num_annuli; i++)
			{
				if ((annulus_name[i] == action_name) || (annulus_name_inverse[i] == action_name))
				{
					if (drilling)
						annuli[i].drill();
					else
						annuli[i].twist((annulus_name[i] == action_name) ? plus : minus);
					
					performed_action = true;
					break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
				}
			}
			if (performed_action) continue;
			
			for (int i = 0; i < num_rectangles; i++)
			{
				if ((rectangle_name[i] == action_name) || (rectangle_name_inverse[i] == action_name))
				{
					if (drilling)
						rectangles[i].drill();
					else
						rectangles[i].half_twist((rectangle_name[i] == action_name) ? plus : minus);
					
					performed_action = true;
					break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
				}
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
		// If we were given an output file, write it there.
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
