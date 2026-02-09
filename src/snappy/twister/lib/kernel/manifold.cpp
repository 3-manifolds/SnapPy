#include <string>
#include <vector>
#include <cstdlib>
#include "global.h"
#include "parsing.h"
#include "twister.h"

void initialise_manifold(manifold &M, 
	std::vector<annulus> &annuli,
	std::vector<rectangle> &rectangles,
	std::vector<std::string> &annulus_name,
	std::vector<std::string> &annulus_name_inverse,
	std::vector<std::string> &rectangle_name,
	std::vector<std::string> &rectangle_name_inverse,
	std::vector<std::string> &macro_name,
	std::vector<std::string> &macro,
	std::string surface_file_contents)
{
	// Check that we actually got the contents of a surface file.
	if (surface_file_contents.substr(0,24) != "# A Twister surface file")
		output_error("A surface file must start with \"# A Twister surface file\".");
	
	std::vector<std::string> surface_description;
	
	{
		size_t marker = 0;
		std::string parsed_line = "";
		surface_file_contents += "\n";  // We add one extra new line, just to make sure / our lives easier.
		for (std::string line = find_next_substring(surface_file_contents, "\n", marker); marker != 0; line = find_next_substring(surface_file_contents, "\n", marker))
		{
			size_t marker2 = 0;
			parsed_line += find_next_substring(line, comment, marker2);  // Remove comments.
			if ((count_substring(line, comment) != 0) && (parsed_line != ""))
			{
				std::string cleaned_line = remove_whitespace(parsed_line);  // Clean up the information that we got.
				if (cleaned_line != "") surface_description.push_back(cleaned_line);
				output_debugging("Parsed line: " + cleaned_line + "\n", 1);
				parsed_line = "";
			}
		}
	}
	
	if (surface_description.size() == 0) output_error("Empty surface description.");
	
	// Determine the number of cubes needed.
	// int num_cubes = atoi((char *) surface_description[0].c_str());
	
	// Determine the number of cubes needed by finding the largest cube number referenced.
	int num_cubes = 0;
	int size = int(surface_description.size());
	for (int i = 0; i < size; i++)
	{
		size_t marker = 0;
		std::string line = surface_description[i];
		
		// Get information about this line.
		std::string type = find_next_substring(line, delimiter, marker);
		if ((type == "annulus") || (type == "rectangle"))
		{
			std::string curve_name = find_next_substring(line, delimiter, marker);
			std::string curve_name_inverse = find_next_substring(line, delimiter, marker);
			
			int sq_num = count_substring(line, delimiter) - 2;  // How many squares there are in this curve.
			if (sq_num < 1) output_error("Invalid curve description.");
			
			// Build the chain.
			for (int j = 0; j < sq_num; j++)
			{
				std::string cube_description = find_next_substring(line, delimiter, marker);
				
				bool orientation = true;
				int number = extract_info(cube_description, orientation);
				if (number+1 >  num_cubes)
					num_cubes = number + 1;
			}
		}
	}
	
	if (num_cubes == 0) output_error("No cubes requested.");
	
	for (int i = 0; i < num_cubes; i++)
		M.cubes.push_back(new cube(M));
	
	// If the manifold we're in is a surface bundle, close
	// up the cube (top to bottom) to get a solid torus. 
	if (M.get_manifold_type() == bundle)
		for (int i = 0; i < num_cubes; i++)
			M.cubes[i]->closeup_cube();
	
	// For storing the information needed to calculate the boundary component
	// indices of the rectangle names.
	std::vector<int> first_number, last_number;  // The first / last cube number of the thickened rectangles.
	std::vector<bool> first_orientation, last_orientation;  // The first / last cube orientations of the thickened rectangles.
	
	std::vector<bool> glued_edges (4*num_cubes, false);
	std::vector<int> vertex_indices;
	vertex_indices.reserve(4*num_cubes);
	for (int i = 0; i < 4*num_cubes; i++)
		vertex_indices.push_back(i);
	
	// Read in each line.
	for (int i = 0; i < size; i++)
	{
		size_t marker = 0;
		std::string line = surface_description[i];
		
		// Get information about this line.
		std::string type = find_next_substring(line, delimiter, marker);
		
		// We extract lines labelled "macro" - these are for parsing shortcuts.
		if (type == "macro")
		{
			macro_name.push_back(find_next_substring(line, delimiter, marker));  // curve_name
			macro.push_back(find_next_substring(line, delimiter, marker));  
		}
		else if ((type == "annulus") || (type == "rectangle"))
		{
			std::string curve_name = find_next_substring(line, delimiter, marker);
			std::string curve_name_inverse = find_next_substring(line, delimiter, marker);
			
			std::vector<int> chain_indices;
			std::vector<cube*> chain_cubes;
			std::vector<bool> chain_upright;
			
			int sq_num = count_substring(line, delimiter) - 2;  // How many squares there are in this curve.
			if (sq_num < 1) output_error("Invalid curve description.");
			chain_cubes.reserve(sq_num);
			chain_upright.reserve(sq_num);
			
			// Build the chain.
			for (int j = 0; j < sq_num; j++)
			{
				std::string cube_description = find_next_substring(line, delimiter, marker);
				
				bool orientation = true;
				int number = extract_info(cube_description, orientation);
				if (number >= num_cubes) output_error(curve_name + " uses a cube that does not exist.");
				chain_indices.push_back(number);
				chain_cubes.push_back(M.cubes[number]);
				chain_upright.push_back(orientation);
			}
			
			// Glue up the edges / vertices as needed.
			for (int j = 0; j < ((type == "annulus") ? sq_num : sq_num-1); j++)
			{
				glued_edges[4*chain_indices[j] + ((chain_upright[j]) ? 1 : 0)] = true;
				glued_edges[4*chain_indices[(j+1) % sq_num] + ((chain_upright[(j+1) % sq_num]) ? 3 : 2)] = true;
				
				int a = vertex_indices[4*chain_indices[j] + ((chain_upright[j]) ? 1 : 0)];
				int b = vertex_indices[4*chain_indices[(j+1) % sq_num] + ((chain_upright[(j+1) % sq_num]) ? 0 : 3)];
				
				// Replace max(a,b) with min(a,b).
				for (int k = 0; k < 4*num_cubes; k++)
					if (vertex_indices[k] == ((a < b) ? b : a)) vertex_indices[k] = ((a < b) ? a : b);
				
				a = vertex_indices[4*chain_indices[j] + ((chain_upright[j]) ? 2 : 1)];
				b = vertex_indices[4*chain_indices[(j+1) % sq_num] + ((chain_upright[(j+1) % sq_num]) ? 3 : 2)];
				
				// Replace max(a,b) with min(a,b).
				for (int k = 0; k < 4*num_cubes; k++)
					if (vertex_indices[k] == ((a < b) ? b : a)) vertex_indices[k] = ((a < b) ? a : b);
			}
			
			// Store it according to its type.
			if (type == "annulus")
			{
				annulus_name.push_back(curve_name);
				annulus_name_inverse.push_back(curve_name_inverse);
				annuli.push_back(annulus(chain_cubes, chain_upright));
			}
			else if (type == "rectangle")
			{
				rectangle_name.push_back(curve_name);
				rectangle_name_inverse.push_back(curve_name_inverse);
				rectangles.push_back(rectangle(chain_cubes, chain_upright));
				
				// Push the first / last cube info onto the relevant lists.
				first_orientation.push_back(chain_upright[0]);
				first_number.push_back(chain_indices[0]);
				last_orientation.push_back(chain_upright[sq_num-1]);
				last_number.push_back(chain_indices[sq_num-1]);
			}
		}
		else
		{
			output_error("Invalid curve or macro description.");
		}
	}
	
	// Now we collate together all of the gluing information to determine the
	// boundary component index of each unglued edge.
	std::vector<int> boundary_index (4*num_cubes, -1);
	int boundary_components = 0;
	for (int i = 0; i < 4*num_cubes; i++)
	{
		if (glued_edges[i]) continue;
		if (boundary_index[i] != -1) continue;
		
		int start_index = vertex_indices[i];
		int current_index = start_index;
		do
		{
			for (int j = 0; j < 4*num_cubes; j++)
			{
				if (glued_edges[j]) continue;
				if (boundary_index[j] != -1) continue;
				if (vertex_indices[j] != current_index) continue;
				
				current_index = vertex_indices[j + ((j % 4 == 3) ? -3 : 1)];  // Move on to the next vertex
				boundary_index[j] = boundary_components;  // Remember to write down the boundary index for that component.
				break;
			}
		} while (current_index != start_index);
		boundary_components++;
	}
	
	// Create the needed marked points.
	for (int i = 0; i < boundary_components; i++)
		M.marked_points.push_back(half_twisted);
	
	// Set the front and back boundary component indices of each rectangle.
	for (int i = 0; i < int(rectangles.size()); i++)
	{
		rectangles[i].set_front(boundary_index[4*first_number[i] + ((first_orientation[i]) ? 3 : 2)]);
		rectangles[i].set_back(boundary_index[4*last_number[i] + ((last_orientation[i]) ? 1 : 0)]);
		rectangles[i].set_one_ended(rectangles[i].get_front() == rectangles[i].get_back());
	}
	
	// Collate all names together.
	std::vector<std::string> arc_names;
	arc_names.insert(arc_names.end(), annulus_name.begin(), annulus_name.end());
	arc_names.insert(arc_names.end(), annulus_name_inverse.begin(), annulus_name_inverse.end());
	arc_names.insert(arc_names.end(), rectangle_name.begin(), rectangle_name.end());
	arc_names.insert(arc_names.end(), rectangle_name_inverse.begin(), rectangle_name_inverse.end());
	arc_names.insert(arc_names.end(), macro_name.begin(), macro_name.end());
	
	// And check that they are a valid collection of names.
	check_valid_names(arc_names);
	
	return;
}

void construct_manifold(manifold &M, std::string surface_file_contents, std::string gluing, std::string handles)
{
	// These will hold all the pieces of information for the manifold.
	std::vector<annulus> annuli;
	std::vector<rectangle> rectangles;
	std::vector<std::string> annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, macro_name, macro;
	
	// Fill up the manifold with interesting information from the surface_file_contents.
	initialise_manifold(M, annuli, rectangles, annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, macro_name, macro, surface_file_contents);
	
	// Format the commands.
	format_command(handles, annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, macro_name, macro);
	format_command(gluing, annulus_name, annulus_name_inverse, rectangle_name, rectangle_name_inverse, macro_name, macro);
	
	output_debugging("\nParsed handles: " + handles + "\nParsed gluing: " + gluing + "\n", 1);
	
	int num_annuli = int(annulus_name.size());
	int num_rectangles = int(rectangle_name.size());
	
	if (handles != "")  // Attach handles (if needed).
	{
		size_t marker = 0;
		
		do
		{
			std::string handle_location = find_next_substring(handles, separator, marker);
			
			if (handle_location == "") continue;  // Skip it if there is no action.
			
			bool performed_action = false;  // Records if we actually performed an action.
			for (int i = 0; i < num_annuli; i++)
			{
				if ((annulus_name[i] == handle_location) || (annulus_name_inverse[i] == handle_location))
				{
					annuli[i].twohandle((annulus_name[i] == handle_location));
					performed_action = true;
					break;  // Annuli / rectangle names are unique, so we can give up if we get a hit.
				}
			}
			if (performed_action) continue;
			
			output_warning("No annulus named '" + handle_location + "' found to attach a handle to.");
		} while (marker != 0);
	}
	
	M.insert_layer();  // Insert layer #1 to prevent rooftop access.
	
	// Reset the layer flags on cubes and marked points.
	for (int i = 0; i < int(M.marked_points.size()); i++)
		M.marked_points[i] = half_twisted;
	
	for (int i = 0; i < int(M.cubes.size()); i++)
		M.cubes[i]->set_status(false);
	
	// We parse the gluing w = w_1 ... w_n from left to right, so our twists are 
	// stacked from top to bottom. So when a curve c at the bottom is pushed up 
	// the resulting curve at the top is w(c) = w_1( ... w_n(c) ... ), i.e. mapping
	// classes act on the left.
	if (gluing != "")  // Twist / drill (if needed).
	{
		size_t marker = 0;
		
		do
		{
			std::string action_name = find_next_substring(gluing, separator, marker);
			
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
			output_warning("No annulus or rectangle named '" + action_name + "' found to twist or half-twist.");
		} while (marker != 0);
	}
	
	M.tidy_boundary();
	
	return;
}
