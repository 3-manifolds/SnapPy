
#include <string>
#include <iostream>
#include <fstream>

#include "global.h"

// Global options (and their defaults).
int GLOBAL_debugging_level = 0;
bool GLOBAL_warnings = true;
bool GLOBAL_optimise = true;
bool GLOBAL_calculate_peripheral_curves = true;
std::string GLOBAL_message_stream = "";

void set_globals_to_defaults()
{
	// Reset the global options to their defaults.
	GLOBAL_debugging_level = 0;
	GLOBAL_warnings = true;
	GLOBAL_optimise = true;
	GLOBAL_calculate_peripheral_curves = true;
	GLOBAL_message_stream = "";
}

// Global output management.
void output_debugging(const std::string strn, const int level)
{
	if (GLOBAL_debugging_level >= level)
		GLOBAL_message_stream += " " + strn;
	
	return;
}

void output_warning(const std::string strn)
{
	if (GLOBAL_warnings)
		GLOBAL_message_stream += "Warning: " + strn + "\n";
	
	return;
}

void output_error(const std::string strn)
{
	if (GLOBAL_warnings)
		GLOBAL_message_stream += "Error: " + strn + "\n";
	
	throw -1;
	return;
}

