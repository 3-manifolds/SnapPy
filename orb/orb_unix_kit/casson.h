#ifndef _casson_
#define _casson_

#include "casson_typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern Boolean verify_casson(
        CassonFormat *cf);

extern void free_casson(
        CassonFormat * cf);

extern Triangulation * casson_to_triangulation(
        CassonFormat * cf);

#ifdef __cplusplus
}
#endif

#endif
