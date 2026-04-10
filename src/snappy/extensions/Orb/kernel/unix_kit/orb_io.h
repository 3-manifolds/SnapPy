/**
 * orb_io.h
 *
 * Functions to read and write the orb file format used by
 * Damian Heard's original orb.
 *
 * It simply concatenates:
 * - The header: % orb
 * - The triangulation in a Casson-like format, see casson_io.h
 * - If given, the planar diagram of the knotted graph, see diagram_io.h
 */

#ifndef _orb_io_
#define _orb_io_

#include "diagram.h"

#include "SnapPea.h"

/* Corresponds to Organizer::loadOrbSlot in gui/organizer.cpp */
extern void read_orb_from_string(
    char *file_data,
    Triangulation ** trig,
    Diagram ** diagram);

extern void read_orb(
    const char *file_name,
    Triangulation ** trig,
    Diagram ** diagram);

/* Corresponds to ManifoldInterface::saveSlot in gui/manifold_interface.cpp */
char * write_orb_to_string(
    Triangulation *trig,
    Diagram * diagram);

#endif
