#include "orb_io.h"

#include "casson_io.h"
#include "diagram_io.h"
#include "ostream.h"

#include "kernel.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

/* What is happening in the original Orb

   loadOrbSlot (in gui/organizer.cpp):
       * reads the line "% orb"
       * reads the next line for the name (but ignores it)
       - readTriangulation
            - readCassonFormat
            - verifyCassonFormat
            - cassonToTriangulation
            - freeCassonFormat
       - readDiagram
            - DiagramCanvas::readDiagram (in interface.cpp)
                 * seems to be somewhat similar to the plink format (vertices, edges, crossings)
                 * EndData, Edge, Vertex, Crossing in diagram_canvas.h
                 * DiagramCanvas::outputTriangulation creates a Graph

   ManifoldInterface::saveSlot (in gui/manifold_interface.cpp)
       - Console::saveTriangulation (in gui/console.cpp)
       - DiagramCanvas::saveDiagram (in gui/interface.cpp)

*/

static Boolean is_eol_char(char c){
    return c == '\n' ||  c == '\r' || c == '\0';
}

void read_orb_from_string(
    char *str,
    Triangulation ** trig,
    Diagram ** diagram)
{
    char * p = str;

    /*
     * Read and ignore the header (% orb).
     */
    while (!is_eol_char(*p)) {
        p++;
    }

    /*
     * Find first non-empty line.
     */

    while (is_eol_char(*p) && *p != '\0') {
        p++;
    }

    char * name_start = p;
    while (!is_eol_char(*p)) {
        p++;
    }
    size_t name_length = p - name_start;

    *trig = read_casson_format(&p);

    if (*trig && name_length > 0)
    {
        (*trig)->name = (char*) my_malloc(name_length + 1);
        if ((*trig)->name == NULL)
        {
            /* uFatalError */
            return;
        }
        memcpy((*trig)->name, name_start, name_length);
        (*trig)->name[name_length] = '\0';
    }
    while (isspace(*p)) {
        p++;
    }
    if (*p == '\0') {
        return;
    }

    *diagram = read_diagram_from_string(p);
}

void read_orb(
    const char *file_name,
    Triangulation ** trig,
    Diagram ** diagram)
{
    /* Follows unit_kit/unix_file_io.c */

    FILE * fp = fopen(file_name, "rb");
    if (fp == NULL) {
        return;
    }

    long filesize;
    if ( fseek(fp, 0, SEEK_END) != 0 ||
         (filesize = ftell(fp) ) == -1 ||
         fseek(fp, 0, SEEK_SET) != 0) {
        return;
    }

    char * buffer = (char*) malloc(filesize + 1);
    buffer[filesize] = '\0';

    if ( fread(buffer, filesize, 1, fp) != 1) {
        free(buffer);
        return;
    }

    read_orb_from_string(buffer, trig, diagram);

    free(buffer);
}

static void write_orb_to_stream(
    OStream * stream,
    Triangulation *trig,
    Diagram * diagram)
{
    ostream_printf(stream, "%% orb\n");
    if (trig)
    {
        if (trig->name) {
            ostream_printf(stream, "%s\n", trig->name);
        } else {
            ostream_printf(stream, "untitled\n");
        }
        write_casson_format_to_stream(
            stream,
            trig,
            FALSE, TRUE, TRUE);
    }

    if (diagram)
    {
        if (trig) {
            ostream_printf(stream, "\n");
        }

        write_diagram_to_stream(
            stream,
            diagram);
    }
}

char * write_orb_to_string(
    Triangulation *trig,
    Diagram * diagram)
{
    OStream stream;
    string_stream_init(&stream);

    write_orb_to_stream(&stream, trig, diagram);

    return stream.buffer;
}
