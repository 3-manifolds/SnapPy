#include "parse_util.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

int main()
{
    char *a = my_strdup("Hello World\n\rNext line\n  3  4\n4  5  \n");
    char *l, *t;

    while ((l = parse_line(&a))) {
        printf("Line: '%s'\n", l);
	
	while ((t = parse_token(&l))) {
  	    printf("           Token: '%s'\n", t);
	}
    }
  
    printf("End\n");

    return 0;
}
