#ifndef PARSING_HEADER_H
#define PARSING_HEADER_H

#include <string>
#include <vector>
#include "twister.h"

std::string load_file_contents(const std::string file);

std::string find_next_substring(const std::string inpt, const std::string search_for, size_t &start_point);
int count_substring(const std::string inpt, const std::string search_for);
int extract_info(const std::string inpt, bool &orientation);
std::string remove_whitespace(const std::string inpt);
void check_valid_names(const std::vector<std::string> &arc_names);

void format_command(std::string &command,
	const std::vector<std::string> &annulus_name,
	const std::vector<std::string> &annulus_name_inverse,
	const std::vector<std::string> &rectangle_name,
	const std::vector<std::string> &rectangle_name_inverse,
	const std::vector<std::string> &macro_name,
	const std::vector<std::string> &macro);

#endif
