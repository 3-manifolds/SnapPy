#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
#include "twister.h"
#include "global.h"

std::string load_file_contents(const std::string file)
{
	std::ifstream myfile(file.c_str());
	if (! myfile.is_open())
		output_error("Unknown file requested.");
	
	std::string line = "";
	std::string contents = "";
	while (myfile.good())
	{
		getline (myfile, line);
		contents += line + "\n";
	}
	myfile.close();
	
	return contents;
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
	orientation = (inpt.substr(0,1) == "-") ? false : true;
	
	return abs(atoi((char *) inpt.c_str()));
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
			output_error("Empty curve or macro name.");
	
	// Check for names that contain a non-valid character or start with a digit.
	for (int i = 0; i < size; i++)
		if ((arc_names[i].find_first_not_of(valid_arc_name_characters) != std::string::npos) || isdigit(arc_names[i][0]) || arc_names[i][0] == '_')
			output_error("Invalid curve or macro name.");
	
	for (int i = 0; i < size; i++)
		for (int j = i+1; j < size; j++)
			if (arc_names[i]== arc_names[j])
				output_error("Duplicated curve or macro name.");
	
	return;
}

void find_and_replace(std::string &inpt, const std::string &search, const std::string &replace)
{
	for(size_t next = inpt.find(search); next != std::string::npos; next = inpt.find(search, next))
	{
		inpt.replace(next, search.length(), replace);   // Do the replacement.
		next += replace.length();                       // Move to just after the replace This is the point were we start the next search from. 
	}
}

void recursive_find_and_replace(std::string &command, const std::string &search, const std::string &replace)
{
	while (command.find(search) != std::string::npos)
		find_and_replace(command, search, replace);
	
	return;
}

int digit_length(int n)
{
	int c = 0;
	if (n <= 0) c++;
	while (n != 0)
	{
	   n  /= 10;
	   c++;
	}
	return c;
}

void nicen_command(std::string &command)
{
	command = separator + remove_whitespace(command) + separator;
	find_and_replace(command, "(", separator + "(" + separator);
	find_and_replace(command, ")", separator + ")" + separator);
	find_and_replace(command, ")" + separator + "^", ")^");
	find_and_replace(command, drill, separator + drill);
	
	size_t marker = 0;
	for (find_next_substring(command, ")", marker); marker != 0; find_next_substring(command, ")", marker))
		if (command.substr(marker,1) != "^")
			command.replace(marker, 0, "^1");
	
	for (find_next_substring(command, "^", marker); marker != 0; find_next_substring(command, "^", marker))
	{
		if ((! isdigit(command[marker])) && (command[marker] != '-'))
			output_error("In commands ^ must be followed by an integer power.");
		if (command.substr(marker-2,1) != ")")
			output_error("Substrings of a command must be enclosed in ( ) to raise them to a power.");
	}
	
	for (find_next_substring(command, drill, marker); marker != 0; find_next_substring(command, drill, marker))
		if (valid_arc_name_characters.find(command[marker]) == std::string::npos)
			output_error("In commands ! must be followed by a curve to drill.");
	
	int c = 0;
	for (find_next_substring(command, "()", marker); marker != 0; find_next_substring(command, "()", marker))
	{
		if (command.substr(marker-1,1) == "(") c++;
		if (command.substr(marker-1,1) == ")") c--;
		if (c < 0) output_error("Unbalanced braces.");
	}
	if (c != 0) output_error("Unbalanced braces.");
	
	return;
}

void tokenise_command(std::string &command, const std::set<std::string> &valid_names)
{
	std::string new_command = "";
	std::string new_subcommand = "";
	size_t marker = 0;
	
	do
	{
		if (command.substr(marker,1) == drill)
		{
			new_command += separator + drill;
			marker++;
		}
		if (command.substr(marker,1) == "(")
		{
			new_command += separator + "(";
			marker++;
		}
		if (command.substr(marker,1) == ")")
		{
			int power = atoi((char *) command.substr(marker).c_str());
			int power_length = digit_length(power);
			
			new_command += separator + command.substr(marker, 2+power_length);
			marker += 2+power_length;
		}
		
		size_t old_marker = marker;
		new_subcommand = find_next_substring(command, separator, marker);
		
		for (int j = new_subcommand.size(); j >= 0; j--)
		{
			if (j > 0)
			{
				if (valid_names.find(new_subcommand.substr(0,j)) != valid_names.end())
				{
					new_command += new_subcommand.substr(0,j) + separator;
					marker = old_marker + j;
					break;
				}
			}
			else
				new_command += new_subcommand + separator;
		}
	} while (marker != 0);
	
	command = new_command;
	
	return;
}

void replace_macros(std::string &command, const std::vector<std::string> &macro_name, const std::vector<std::string> &macro)
{
	std::string old_command = "";
	int num_macros = int(macro_name.size());
	
	for (int i = 0; i < num_macros; i++)
		if (command.find(drill + macro_name[i]) != std::string::npos)
			output_error("Cannot drill a macro.");
	
	int c = 0;  // Count the number of replacements we make.
	while (command != old_command)
	{
		old_command = command;
		for (int i = 0; i < num_macros; i++)
			find_and_replace(command, separator + macro_name[i] + separator, separator + macro[i] + separator);
		
		if (c++ > num_macros)
			output_error("Circular references detected in macros in command.");
	}
	
	return;
}

void expand_brackets(std::string &command, std::vector<std::string> &names, std::vector<std::string> &inverse_names)
{
	int num_names = int(names.size());
	
	for(size_t next = command.find(")"); next != std::string::npos; next = command.find(")", next))
	{
		size_t previous = command.rfind("(", next);
		std::string replacement = "", subreplacement = "";
		int power = atoi((char *) command.substr(next+2, std::string::npos).c_str());
		int power_length = digit_length(power);
		
		subreplacement = command.substr(previous+1, next-previous-1);
		if (power < 0)
		{
			size_t marker = 0;
			std::string new_subreplacement = "";
			for (std::string s = find_next_substring(subreplacement, separator, marker); marker != 0; s = find_next_substring(subreplacement, separator, marker))
			{
				for (int i = 0; i < num_names; i++)
				{
					if (names[i] == s)
					{
						s = inverse_names[i];
						break;
					}
					if (drill + names[i] == s)
					{
						s = drill + inverse_names[i];
						break;
					}
				}
				new_subreplacement = s + separator + new_subreplacement; // Restack the results.
			}
			subreplacement = new_subreplacement;
		}
		
		power = abs(power);
		
		for (int i = 0; i < power; i++)
			replacement += subreplacement;
		
		command.replace(previous, next-previous+1+1+power_length, replacement);
		
		next++;
	}
	
	return;
}

void cleanup_command(std::string &command)
{
	recursive_find_and_replace(command, separator+separator, separator);
	
	return;
}

void format_command(std::string &command,
	const std::vector<std::string> &annulus_name,
	const std::vector<std::string> &annulus_name_inverse,
	const std::vector<std::string> &rectangle_name,
	const std::vector<std::string> &rectangle_name_inverse,
	const std::vector<std::string> &macro_name,
	const std::vector<std::string> &macro)
{
	// We now combine together all the curve names for formatting.
	std::set<std::string> valid_names;
	valid_names.insert(annulus_name.begin(), annulus_name.end());
	valid_names.insert(annulus_name_inverse.begin(), annulus_name_inverse.end());
	valid_names.insert(rectangle_name.begin(), rectangle_name.end());
	valid_names.insert(rectangle_name_inverse.begin(), rectangle_name_inverse.end());
	valid_names.insert(macro_name.begin(), macro_name.end());
	
	std::vector<std::string> names, inverse_names;
	names.insert(names.end(), annulus_name.begin(), annulus_name.end());
	names.insert(names.end(), annulus_name_inverse.begin(), annulus_name_inverse.end());
	names.insert(names.end(), rectangle_name.begin(), rectangle_name.end());
	names.insert(names.end(), rectangle_name_inverse.begin(), rectangle_name_inverse.end());
	inverse_names.insert(inverse_names.end(), annulus_name_inverse.begin(), annulus_name_inverse.end());
	inverse_names.insert(inverse_names.end(), annulus_name.begin(), annulus_name.end());
	inverse_names.insert(inverse_names.end(), rectangle_name_inverse.begin(), rectangle_name_inverse.end());
	inverse_names.insert(inverse_names.end(), rectangle_name.begin(), rectangle_name.end());
	
	nicen_command(command);
	tokenise_command(command, valid_names);
	replace_macros(command, macro_name, macro);
	expand_brackets(command, names, inverse_names);
	cleanup_command(command);
	
	return;
}
