#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>

// Global options.
extern bool GLOBAL_debugging;
extern bool GLOBAL_warnings;
extern bool GLOBAL_optimise;
extern bool GLOBAL_calculate_peripheral_curves;

// Global compile time constants.
const int max_name_length = 99;  // This is the maximum length name that SnapPea can handle. See unix_file_io.c (L:208).
const std::string max_name_length_string = "99";

// Globals for reading from a file.
const std::string delimiter = ",";
const std::string comment = "#";
const std::string valid_arc_name_characters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_";

// Globals for parsing word.
// These MUST be single characters.
const std::string seperator = "*";
const std::string start_drill = "[";
const std::string start_twist = "]";
const std::string product = seperator + start_drill + start_twist;
const std::string power = "^";

void output_debugging(std::string strn);
void output_warning(std::string strn);
void output_error(std::string strn);

#endif
