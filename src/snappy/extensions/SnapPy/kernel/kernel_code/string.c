/*
 *  string.c
 *
 *  This file provides the function
 *
 *      char *my_strdup(CONST char *)
 *
 *  my_strdup is a version of the POSIX strdup using my_malloc.
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

char *my_strdup(const char *s1)
{
    if (s1 == NULL)
        return NULL;

    char *result = NEW_ARRAY(strlen(s1) + 1, char);
    strcpy(result, s1);
    return result;
}

SNAPPEA_NAMESPACE_END_SCOPE


