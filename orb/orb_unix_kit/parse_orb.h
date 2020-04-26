#ifndef _parse_orb_
#define _parse_orb_

#include "casson_typedefs.h"

extern void read_orb_from_string(
        char *file_data,
        char **name,
        CassonFormat ** cf,
        char **orb_link_projection_data);

extern void read_orb(
        const char *file_name,
        char **name,
        CassonFormat ** cf,
        char **orb_link_projection_data);

#endif

