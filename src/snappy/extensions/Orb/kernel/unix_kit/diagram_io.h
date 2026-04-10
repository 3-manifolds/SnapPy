/**
 * diagram_io.h
 *
 * Functions for reading and writing a Diagram in a format
 * similar to the link projection format.
 *
 */

#ifndef _diagram_io_
#define _diagram_io_

typedef struct Diagram Diagram;
typedef struct OStream OStream;

Diagram * read_diagram_from_string(char *str);
char * write_diagram_to_string(Diagram *diagram);
void write_diagram_to_stream(OStream * stream, Diagram *diagram);

#endif
