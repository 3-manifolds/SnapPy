/**
 * casson_io.h
 *
 * Functions for reading and writing a triangulation format similar to the one
 * used by Casson's geo.
 *
 * The format starts with the name or the triangulation,
 * an optional "vertices_known" and one to three blocks (depending
 * on the configuration flags) with each block having one line
 * per edge.
 *
 */

#ifndef _casson_io_
#define _casson_io_

#include "SnapPea.h"

typedef struct OStream OStream;

/*
 * Reads the triangulation from the string starting at *str and advances
 * *str to point to the end of the section corresponding to the triangulation.
 *
 * This can be used if additional information was simply concatenated to
 * the triangulation format.
 */
Triangulation *read_casson_format(char **str);

void write_casson_format_to_stream(
    OStream *stream,
    Triangulation *manifold,
    Boolean include_angular_error,
    Boolean include_geometric_structure_and_cusps,
    Boolean include_peripheral_curves);

char *write_casson_format_to_string(
    Triangulation *manifold,
    Boolean include_angular_error,
    Boolean include_geometric_structure_and_cusps,
    Boolean include_peripheral_curves);

#endif
