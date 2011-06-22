#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include "twister.h"
#include "global.h"

// Please make sure that this help message and the function
// parse_input remain in sync.  (It would be nice if the various
// clauses of print_help messages were actually stored in the various
// clauses of parse_input...)

void print_help(std::ostream &o)
{
	o << std::endl;
	o << "Twister:  A program for generating SnapPea triangulation files from" << std::endl;
	o << "mapping class group data." << std::endl;
	o << std::endl;
	// In different environments the commands will look different.
	#if OS_ENVIRONMENT == 0  // Windows.
		o << "Usage: Twister.exe [options]" << std::endl;
	#elif OS_ENVIRONMENT == 1  // UNIX.
		o << "Usage: Twister.out [options]" << std::endl;
	#endif
	o << std::endl;
	o << "Options:" << std::endl;
	o << " --help   Display this message." << std::endl;
	o << " --version   Display the current version number." << std::endl;
	o << std::endl;
	o << " -f, --surface 'surface'   Use the specified surface file."  << std::endl;
	o << " -o, --output 'file'   Write results to the specified file."  << std::endl;
	o << " -ow, --output-warnings 'file' Write warnings and errors to the specified file." << std::endl;
	o << " -b, --bundle 'monodromy'   Create a surface bundle." << std::endl;
	o << " -s, --splitting 'gluing'   Create a Heegaard splitting." << std::endl;
	o << " -h, --handles 'handles'   Attach 2-handles." << std::endl;
	o << " -n, --name 'name'   Name of the 3-manifold." << std::endl;
	o << std::endl;
	o << " -ml   Turn off longitude and meridian calculations." << std::endl;
	o << " -w, --warnings   Turn off warnings." << std::endl;
	o << " -O, --optimisations   Turn off optimisations." << std::endl;
	o << " -d, --debugging   Turn on debugging mode." << std::endl;
	o << std::endl;
	o << "'surface' is the path to a .sur file and must be provided." << std::endl;
	o << std::endl;
	o << "'monodromy' and 'gluing' are words of annulus and rectangle names (or" << std::endl;
	o << "inverses). These are read from left to right and determine a sequence of" << std::endl;
	o << "(half) Dehn twists and drillings.  For example, 'a*B*a*B*A*A*[a*b]' will" << std::endl;
	o << "perform 6 twists and then drill twice." << std::endl;
	o << std::endl;
	o << "'handles' is a word of annulus names (or inverses).  For example, 'a*c*A'" << std::endl;
	o << " means attach three 2-handles, two above and one below." << std::endl;
	o << std::endl;
	o << "Examples:" << std::endl;
	// In different environments the examples will look different.
	#if OS_ENVIRONMENT == 0  // Windows.
		o << "   Twister.exe -f \".\\Surfaces\\S_1_1.sur\" -b \"a*B\" " << std::endl;
		o << "produces a triangulation of the figure eight knot complement." << std::endl;
		o << "   Twister.exe -f \".\\Surfaces\\S_2.sur\" -s \"\" -h \"a*B*c\" " << std::endl;
		o << "produces a triangulation of the genus two splitting of the solid torus." << std::endl;
	#elif OS_ENVIRONMENT == 1  // UNIX.
		o << "   Twister.out -f \"Surfaces/S_1_1.sur\" -b \"a*B\" " << std::endl;
		o << "produces a triangulation of the figure eight knot complement." << std::endl;
		o << "   Twister.out -f \"Surfaces/S_2.sur\" -s \"\" -h \"a*B*c\" " << std::endl;
		o << "produces a triangulation of the genus two splitting of the solid torus." << std::endl;
	#endif
	
	return;
}

void parse_input(int argc, char **argv, std::string &surface_file, std::string &output_file, std::string &name, Manifold_type &manifold_type, std::string &gluing, std::string &handles)
{
	// Set default values.
	surface_file = "";
	output_file = "";
	GLOBAL_message_stream = "";
	name = "";
	manifold_type = bundle;
	gluing = "";
	handles = "";
	
	if (argc == 1)
	{
		print_help(std::cout);
		exit(0);
	}
	
	if (argc >= 2)
	{
		for (int i = 1; i < argc; i++)
		{
			if (strcmp(argv[i], "--help") == 0)
			{
				print_help(std::cout);
				exit(0);
			}
				else if (strcmp(argv[i], "--version") == 0)
			{
				std::cout << "Twister 2.3.0" << std::endl;
				exit(0);
			}
			else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--surface") == 0)
			{
				// Really should check that the next entry is a word, _not_ a flag.
				if (++i == argc) 
					output_error("-f requires a surface file.");
				
				surface_file = argv[i];
			}
			else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0)
			{
				// Really should check that the next entry is a word, _not_ a flag.
				if (++i == argc) 
					output_error("-o requires an output file.");
				
				output_file = argv[i];
			}
			else if (strcmp(argv[i], "-ow") == 0 || strcmp(argv[i], "--output-warning") == 0)
			{
				// Really should check that the next entry is a word, _not_ a flag.
				if (++i == argc) 
					output_error("-ow requires an output warning file.");
				
				GLOBAL_message_stream = argv[i];
			}
			else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--splitting") == 0)
			{
				// Really should check that the next entry is a gluing, _not_ a flag.
				if (++i == argc)
					output_error("-s requires a gluing.");
				
				manifold_type = splitting;
				gluing = argv[i];
			}
			else if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bundle") == 0)
			{
				// Really should check that the next entry is a monodromy, _not_ a flag.
				if (++i == argc)
					output_error("-b requires a monodromy.");
				
				manifold_type = bundle;
				gluing = argv[i];
			}
			else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--handles") == 0)
			{
				// Really should check that the next entry is a handles_list, _not_ a flag.
				if (++i == argc)
					output_error("-h requires a collecton of handles to attach.");
				
				handles = argv[i];
				manifold_type = splitting;
			}
			else if (strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--name") == 0)
			{
				// Really should check that the next entry is a name, _not_ a flag.
				if (++i == argc)
					output_error("-n requires a manifold name.");
				
				name = argv[i];
			}
			else if (strcmp(argv[i], "-ml") == 0)
			{
				GLOBAL_calculate_peripheral_curves = false;
			}
			else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--warnings") == 0)
			{
				GLOBAL_warnings = false;
			}
			else if (strcmp(argv[i], "-O") == 0 || strcmp(argv[i], "--optimisations") == 0)
			{
				GLOBAL_optimise = false;  // Turns off foldoff really.
			}
			else if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--debugging") == 0)
			{
				GLOBAL_debugging = true;
			}
			else
			{
				output_warning("Unexpected input at the command line.");
			}
		}
	}
	
	if (surface_file == "") output_error("No surface file specified.");
	
	// If no name was specified we rebuild the command given (in a nice order).
	if (name == "")
	{
		if (manifold_type == bundle) name = "-b " + gluing;
		else if (manifold_type == splitting) name = "-s " + gluing + " -h " + handles;
		else name = "";
		
		if (GLOBAL_calculate_peripheral_curves == false) name += " -ml ";
		if (GLOBAL_warnings == false) name += " -w ";
		if (GLOBAL_optimise == false) name += " -O ";
		if (GLOBAL_debugging == true) name += " -d ";
		name += " -f " + surface_file;
	}
	
	// SnapPea can't handle long names though.
	if (int(name.size()) > max_name_length)
	{
		output_warning("Name truncated. SnapPea can't handle triangulations with > " + max_name_length_string + " character names.");
		name = name.substr(0, max_name_length);
	}
	
	return;
}

int count_substring(const std::string inpt, const std::string search_for)
{
	int count = 0;
	for (size_t pos = inpt.find(search_for); pos != std::string::npos; count++)
		pos = inpt.find(search_for, pos + 1);
	
	return count;
}

std::string find_next_substring(const std::string inpt, const std::string search_for, size_t &start_point)
{
	size_t start = start_point;
	start_point = inpt.find_first_of(search_for, start_point) + 1;
	return inpt.substr(start, start_point - start - 1);
}

int extract_info(const std::string inpt, bool &orientation)
{
	if (inpt.substr(0,1) == "+")
	{
		orientation = true;
		return atoi((char *) (inpt.substr(1)).c_str());
	}
	else if (inpt.substr(0,1) == "-")
	{
		orientation = false;
		return atoi((char *) (inpt.substr(1)).c_str());
	}
	else
	{
		orientation = true;
		return atoi((char *) inpt.c_str());
	}
}

std::string remove_whitespace(const std::string inpt)
{
	std::string output = "";
	output.reserve(inpt.size());
	size_t marker = 0, next_marker = 0;
	
	do
	{
		next_marker = inpt.find_first_of(" \t", marker);  // Could filter out more characters here.
		output += inpt.substr(marker, next_marker - marker);
		marker = next_marker + 1;  // Skip over the bad character.
	} while (marker != 0);  // When next_marker == std::string::npos, marker == std::string::npos + 1 == 0.
	
	return output;
}

void check_valid_names(const std::vector<std::string> &arc_names)
{
	int size = int(arc_names.size());
	
	// Check for blank names.
	for (int i = 0; i < size; i++)
		if (arc_names[i] == "")
			output_error("Empty curve name.");
	
	// Check for names that dont contain 
	for (int i = 0; i < size; i++)
		if (arc_names[i].find_first_not_of(valid_arc_name_characters) != std::string::npos)
			output_error("Invalid curve name.");
	
	for (int i = 0; i < size; i++)
		for (int j = i + 1; j < size; j++)
			if (arc_names[i]== arc_names[j])
				output_error("Duplicated curve name.");
	
	return;
}

void build_manifold(manifold &M, 
	std::vector<square> &squares,
	std::vector<annulus> &annuli,
	std::vector<rectangle> &rectangles,
	std::vector<std::string> &annulus_name,
	std::vector<std::string> &annulus_name_inverse,
	std::vector<std::string> &rectangle_name,
	std::vector<std::string> &rectangle_name_inverse,
	std::string surface_file)
{
	std::vector<std::string> surface_description;
	
	std::ifstream myfile(surface_file.c_str());
	if (not myfile.is_open())
		output_error("Unknown surface requested.");
	
	std::string line = "";
	std::string parsed_line = "";
	while (myfile.good())
	{
		getline (myfile, line);
		size_t marker = 0;
		parsed_line += find_next_substring(line, comment, marker);  // Remove comments.
		if ((count_substring(line, comment) != 0) && (parsed_line != ""))
		{
			std::string cleaned_line = remove_whitespace(parsed_line);  // Clean up the information that we got.
			if (cleaned_line != "") surface_description.push_back(cleaned_line);
			output_debugging("Parsed line: " + cleaned_line);
			parsed_line = "";
		}
	}
	myfile.close();
	
	// Determine the number of squares in the surface and build them.
	int num_squares = atoi((char *) surface_description[0].c_str());
	for (int i = 0; i < num_squares; i++)
		squares.push_back(square(M));
	
	// Read in each line.
	int size = int(surface_description.size());
	for (int i = 1; i < size; i++)
	{
		size_t marker = 0;
		line = surface_description[i];
		
		// Get information about this line.
		int sq_num = count_substring(line, delimiter) - 2;  // How many squares there are in this curve.
		std::string type = find_next_substring(line, delimiter, marker);
		if (type == "macro") continue;  // We ignore lines labelled "macro" - these are for parsing shortcuts.
		std::string curve_name = find_next_substring(line, delimiter, marker);
		std::string curve_name_inverse = find_next_substring(line, delimiter, marker);
		
		std::vector<square*> chain_squares;
		std::vector<bool> chain_upright;
		
		chain_squares.reserve(sq_num);
		chain_upright.reserve(sq_num);
		
		// Build the chain.
		for (int j = 0; j < sq_num; j++)
		{
			std::string square_description = find_next_substring(line, delimiter, marker);
			
			bool orientation = true;
			int number = extract_info(square_description, orientation);
			if (number >= num_squares) output_error("Insufficient squares built.");
			chain_squares.push_back(&squares[number]);
			chain_upright.push_back(orientation);
		}
		
		// Store it according to its type.
		if (type == "annulus")
		{
			annulus_name.push_back(curve_name);
			annulus_name_inverse.push_back(curve_name_inverse);
			annuli.push_back(annulus(chain_squares, chain_upright));
		}
		else if (type == "rectangle")
		{
			rectangle_name.push_back(curve_name);
			rectangle_name_inverse.push_back(curve_name_inverse);
			rectangles.push_back(rectangle(chain_squares, chain_upright));
		}
		else
		{
			output_error("Invalid arc type.");
		}
	}
	
	// Collate all names together.
	std::vector<std::string> arc_names;
	size = int(annulus_name.size());
	for (int i = 0; i < size; i++)
	{
		arc_names.push_back(annulus_name[i]);
		arc_names.push_back(annulus_name_inverse[i]);
	}
	
	size = int(rectangle_name.size());
	for (int i = 0; i < size; i++)
	{
		arc_names.push_back(rectangle_name[i]);
		arc_names.push_back(rectangle_name_inverse[i]);
	}
	
	// And check that they are a valid collection of names.
	check_valid_names(arc_names);
	
	return;
}
