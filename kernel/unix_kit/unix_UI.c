/*

IMPORTANT NOTE:  No longer used SnapPy interface library.

*/

/*
 *  unix_UI.c
 *
 *  This file contains a quick and dirty implementation of some UI functions
 *  required for using the SnapPea kernel in a stdio.h environment.
 *  It's intended for use by mathematicians who want to call the SnapPea
 *  kernel functions from within their own C code.
 */

#include "SnapPea.h"
#include <stdio.h>
#include <stdlib.h>

#include "kernel_namespace.h"

void uAcknowledge(
    const char *message)
{
    fprintf(stderr, "%s\n", message);
}

void uFatalError(
    const char    *function,
    const char    *file)
{
    fprintf(
        stderr,
        "A fatal error has occurred in the function %s() in the file %s.c.\n",
        function,
        file);

    exit(1);
}

void uAbortMemoryFull(void)
{
    fprintf(stderr, "out of memory\n");
    exit(2);
}

int uQuery(
    const char  *message,
    const int   num_responses,
    const char  *responses[],
    const int   default_response)
{
    /*
     *  If desired you could write this function to obtain a response from the user,
     *  but for now it set up to return the default response, to facilitate batch
     *  computations.
     */
    fprintf(stderr, "Q: %s\nA:  %s\n", message, responses[default_response]);

    return default_response;
}

/*
 *  The "long computation" feature is unused, but we define its
 *  global variables and functions to avoid a link error.
 */
Boolean gLongComputationInProgress,
        gLongComputationCancelled;
void uLongComputationBegins(
    const char    *message,
    Boolean is_abortable)
{
}
FuncResult uLongComputationContinues()
{
    return func_OK;
}
void uLongComputationEnds()
{
}

#include "end_namespace.h"
