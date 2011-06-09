
#include <string>
#include <iostream>
#include <cstdlib>

// Global options (and their defaults).
bool GLOBAL_debugging = false;
bool GLOBAL_warnings = true;
bool GLOBAL_optimise = true;
bool GLOBAL_calculate_peripheral_curves = true;

// Global output management.
void output_debugging(std::string strn)
{
	if (GLOBAL_debugging) std::cerr << strn << " " << std::flush;
	return;
}

void output_warning(std::string strn)
{
	if (GLOBAL_warnings) std::cerr << "Warning: " << strn << std::endl;
	return;
}

void output_error(std::string strn)
{
	// Errors are never optional.
	std::cerr << "Error: " << strn << std::endl;
	exit(1);
}

