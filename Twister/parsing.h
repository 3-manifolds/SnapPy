#ifndef PARSING_HEADER_H
#define PARSING_HEADER_H

#include <string>
#include <vector>
#include "twister.h"

void parse_input(int argc, char **argv, std::string &surface_file, std::string &output_file, std::string &name, Manifold_type &manifold_type, std::string &gluing, std::string &handles);
std::string find_next_substring(std::string inpt, std::string search_for, size_t &start_point);

void build_manifold(manifold &M,
	std::vector<square> &squares,
	std::vector<annulus> &annuli,
	std::vector<rectangle> &rectangles,
	std::vector<std::string> &annulus_name,
	std::vector<std::string> &annulus_name_inverse,
	std::vector<std::string> &rectangle_name,
	std::vector<std::string> &rectangle_name_inverse,
	std::string surface_file);

#endif
