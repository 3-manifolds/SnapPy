#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include "global.h"
#include "parsing.h"
#include "manifold.h"
#include "twister.h"

// Please make sure that this help message and the function
// parse_input remain in sync.  (It would be nice if the various
// clauses of print_help messages were actually stored in the various
// clauses of parse_input...)

void print_help(std::ostream &o)
{
	o << std::endl;
	o << "Twister:  A program for generating SnapPea triangulation files from mapping" << std::endl;
	o << "class group data. By Mark Bell, Tracy Hall and Saul Schleimer." << std::endl;
	o << std::endl;
	// In different environments the commands will look different.
	#if OS_ENVIRONMENT == 0  // Windows.
		o << "Usage: Twister.exe [options]" << std::endl;
	#elif OS_ENVIRONMENT == 1  // UNIX.
		o << "Usage: Twister.out [options]" << std::endl;
	#endif
	o << std::endl;
	o << "Options:" << std::endl;
	o << " --help   Display this information." << std::endl;
	o << " --version   Display the current version number." << std::endl;
	o << std::endl;
	o << " -f, --surface \"surface\"   Use the specified surface file."  << std::endl;
	o << " -o, --output \"file\"   Write results to the specified file."  << std::endl;
	o << " -ow, --output-warnings \"file\" Write warnings and errors to the specified file." << std::endl;
	o << " -b, --bundle \"monodromy\"   Create a surface bundle." << std::endl;
	o << " -s, --splitting \"gluing\"   Create a Heegaard splitting." << std::endl;
	o << " -h, --handles \"handles\"   Attach 2-handles." << std::endl;
	o << " -n, --name \"name\"   Name of the 3-manifold." << std::endl;
	o << std::endl;
	o << " -p   Turn off peripheral curve calculations." << std::endl;
	o << " -w, --warnings   Turn off warnings." << std::endl;
	o << " -O, --optimisations   Turn off optimisations." << std::endl;
	o << " -d, --debugging \"level\"  Set debugging level." << std::endl;
	o << std::endl;
	o << "\"surface\" is the path to a .sur file and must be provided." << std::endl;
	o << std::endl;
	o << "\"monodromy\" and \"gluing\" are words of annulus and rectangle names (or " << std::endl;
	o << "inverses) separated by '*'. These are read from left to right and determine a" << std::endl;
	o << "sequence of (half) Dehn twists and, when prefixed with an \"!\", drillings. " << std::endl;
	o << "For example, \"a*B*a*B*A*A*!a*!b\" will perform 6 twists and then drill twice." << std::endl;
	o << std::endl;
	o << "\"handles\" is a word of annulus names (or inverses) separated by '*'. For " << std::endl;
	o << "example, 'a*c*A' means attach three 2-handles, two above and one below." << std::endl;
	o << std::endl;
	o << "monodromy, gluing and handles all support macros, specified in the surface" << std::endl;
	o << "file, as well as nested powers that are indicated by the format:" << std::endl;
	o << "\t(command)^power" << std::endl;
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

void parse_input(int argc, char **argv, std::string &surface_file, std::string &output_file, std::string &warning_file, std::string &name, Manifold_type &manifold_type, std::string &gluing, std::string &handles)
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
				std::cout << "Twister " + version << std::endl;
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
			else if (strcmp(argv[i], "-ow") == 0 || strcmp(argv[i], "--output-warnings") == 0)
			{
				// Really should check that the next entry is a word, _not_ a flag.
				if (++i == argc) 
					output_error("-ow requires an output warning file.");
				
				warning_file = argv[i];
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
			else if (strcmp(argv[i], "-p") == 0)
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
				if (++i == argc)
					output_error("-d requires a debugging level.");
				if (! isdigit(argv[i][0]))
					output_error("debugging level must be an integer.");
				
				GLOBAL_debugging_level = atoi(argv[i]);
			}
			else
			{
				output_warning("Unexpected input at the command line.");
			}
		}
	}
	
	if (surface_file == "")
		output_error("No surface file specified.");
	
	// If no name was specified we rebuild the command given (in a nice order).
	if (name == "")
	{
		name = "-f \"" + surface_file + "\"";
		if (manifold_type == bundle) name += " -b \"" + gluing + "\"";
		else if (manifold_type == splitting) name += " -s \"" + gluing + "\" -h \"" + handles + "\"";
		
		if (GLOBAL_calculate_peripheral_curves == false) name += " -p ";
		if (GLOBAL_optimise == false) name += " -O ";
	}
	
	// SnapPea can't handle long names though.
	if (int(name.size()) > max_name_length)
	{
		output_warning("Name truncated. SnapPea can't handle triangulations with more that " + max_name_length_string + " character names.");
		name = name.substr(0, max_name_length);
	}
	
	return;
}

int main(int argc, char **argv)
{
	// Set default values.
	set_globals_to_defaults();
	Manifold_type manifold_type = bundle;
	std::string manifold_name = "", handles = "", gluing = "", surface_file = "", surface_file_contents = "", output_file = "", warning_file = "";
	std::string manifold_contents = "";
	bool error_occured = false;
	
	try
	{
		// Parse the given command, extracting out all the interesting bits.
		parse_input(argc, argv, surface_file, output_file, warning_file, manifold_name, manifold_type, gluing, handles);
		surface_file_contents = load_file_contents(surface_file);
		
		// Build the manifold and get its information.
		manifold M(manifold_name, manifold_type);
		construct_manifold(M, surface_file_contents, gluing, handles);
		manifold_contents = M.to_string();
	}
	catch (...)
	{
		error_occured = true;
	}
	
	// Write out any messages stored in the GLOBAL_message_stream.
	try
	{
		if (GLOBAL_message_stream != "")
		{
			if (warning_file != "")
			{
				std::ofstream myfile(warning_file.c_str());
				if (not myfile.is_open())
					output_error("Unable to write to " + warning_file + ".");
				
				myfile << GLOBAL_message_stream;
				myfile.close();
			}
			else  // Else write it to std::cerr.
			{
				std::cerr << GLOBAL_message_stream << std::endl;
			}
		}
	}
	catch (...)
	{
		error_occured = true;
	}
	
	if (error_occured) exit(1);
	
	try
	{
		if (output_file != "")  // If we were given an output file, write it there.
		{
			std::ofstream myfile(output_file.c_str());
			if (not myfile.is_open())
				output_error("Unable to write to " + output_file + ".");
			
			myfile << manifold_contents;
			myfile.close();
		}
		else  // Else write it to std::cout.
		{
			std::cout << manifold_contents;
		}
	}
	catch (...)
	{
		error_occured = true;
	}
	
	return 0;
}
