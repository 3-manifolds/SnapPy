
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "global.h"

// Global options (and their defaults).
bool GLOBAL_debugging = false;
bool GLOBAL_warnings = true;
bool GLOBAL_optimise = true;
bool GLOBAL_calculate_peripheral_curves = true;
std::string GLOBAL_message_stream = "";

void set_globals_to_defaults()
{
	// Reset the global options to their defaults.
	GLOBAL_debugging = false;
	GLOBAL_warnings = true;
	GLOBAL_optimise = true;
	GLOBAL_calculate_peripheral_curves = true;
	GLOBAL_message_stream = "";
}


// Global output management.
void output_debugging(std::string strn)
{
	if (GLOBAL_debugging)
	{
		if (GLOBAL_message_stream != "")
		{
			std::ofstream myfile(GLOBAL_message_stream.c_str(), std::ios::app);
			if (not myfile.is_open()) exit(1);
			myfile << strn << " " << std::flush;
			myfile.close();
		}
		else
			std::cerr << strn << " " << std::flush;
	}
	
	return;
}

void output_warning(std::string strn)
{
	if (GLOBAL_warnings)
	{
		if (GLOBAL_message_stream != "")
		{
			std::ofstream myfile(GLOBAL_message_stream.c_str(), std::ios::app);
			if (not myfile.is_open()) exit(1);
			myfile << "Warning: " << strn << std::endl;
			myfile.close();
		}
		else
			std::cerr << "Warning: " << strn << std::endl;
	}
	return;
}

void output_error(std::string strn)
{
	// Errors are never optional.
	if (GLOBAL_message_stream != "")
	{
		std::ofstream myfile(GLOBAL_message_stream.c_str(), std::ios::app);
		if (not myfile.is_open()) exit(1);
		myfile << "Error: " << strn << std::endl;
		myfile.close();
	}
	else
		std::cerr << "Error: " << strn << std::endl;
	
	exit(1);
}

