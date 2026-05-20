/**
 * diagram_io.h
 *
 * Functions for reading and writing a Diagram in a format
 * similar to the link projection format.
 *
 */

#ifndef _orb_diagram_io_
#define _orb_diagram_io_

typedef struct OrbDiagram OrbDiagram;
typedef struct OStream OStream;

OrbDiagram * orb_read_diagram_from_string(char *str);
char * orb_write_diagram_to_string(OrbDiagram *diagram);
void orb_write_diagram_to_stream(OStream * stream, OrbDiagram *diagram);

#endif
