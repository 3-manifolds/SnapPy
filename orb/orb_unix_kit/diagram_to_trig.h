#ifndef _diag_to_trig_h_
#define _diag_to_trig_h_

typedef struct Triangulation Triangulation;

#ifdef __cplusplus
extern "C" {
#endif

Triangulation * diagram_data_to_triangulation(const char *d);

#ifdef __cplusplus
}
#endif

#endif
