/*
 *  unix_file_io.h
 *
 *  These three functions allow unix-style programs
 *  to read and save Triangulations.
 */

#include "SnapPea.h"

extern Triangulation    *read_triangulation(char *file_name);
extern Triangulation    *read_triangulation_from_string(char *file_data);
extern void             write_triangulation(Triangulation *manifold, char *file_name);
extern char             *string_triangulation(Triangulation *manifold);
