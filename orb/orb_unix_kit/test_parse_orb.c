#include "parse_orb.h"

#include "casson.h"

#include <stdio.h>

int main()
{
    char *name = NULL;
    CassonFormat * casson = NULL;
    char *orb_link_projection_data = NULL;
    Triangulation * trig = NULL;
    Boolean ok = FALSE;

    read_orb("example.orb", &name, &casson, &orb_link_projection_data);

    if (casson) {
        printf("Got it\n");
    }

    if (name) {
        printf("Name: %s\n", name);
    }

    if (orb_link_projection_data) {
        printf("plink: %s\n", orb_link_projection_data);
    }

    if (!verify_casson(casson)) {
        printf("verify_casson failed\n");
        return 1;
    }

    printf("verify_casson passed\n");

    trig = casson_to_triangulation(casson);

    printf("got trig\n");

    printf("Solution type: %d\n", find_structure(trig, FALSE));

    printf("Volume: %lf\n", my_volume(trig, &ok));
    printf("Vol ok: %d\n", ok);

    free_casson(casson);
    
    free_triangulation(trig);

    return 0;
}
