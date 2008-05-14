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

void uAcknowledge(
    const char *message)
{
    fprintf(stderr, "%s\n", message);
}

void uFatalError(
    char    *function,
    char    *file)
{
    fprintf(
        stderr,
        "A fatal error has occurred in the function %s() in the file %s.c.\n",
        function,
        file);

    //    exit(1);
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
 *  Callbacks can be used to set up a signal for user interrupts.
 */
void no_op() {
}

void (*begin_long_comp_callback)() = no_op;
void (*continue_long_comp_callback)() = no_op;
void (*end_long_comp_callback)() = no_op;

void register_callbacks(void(*begin_callback)(),
			void(*middle_callback)(),
			void(*end_callback)()){
  begin_long_comp_callback = begin_callback;
  continue_long_comp_callback = middle_callback;
  end_long_comp_callback = end_callback;
}
  
Boolean gLongComputationInProgress,
        gLongComputationCancelled;

void cancel_computation(){
  gLongComputationCancelled = 1;
}

void uLongComputationBegins(
    char    *message,
    Boolean is_abortable)
{
  gLongComputationCancelled = 0;
  gLongComputationInProgress = 1;
  begin_long_comp_callback();
}

FuncResult uLongComputationContinues()
{
  continue_long_comp_callback();
  if ( gLongComputationCancelled ){
    return func_cancelled;
  }
  return func_OK;
}

void uLongComputationEnds()
{
  end_long_comp_callback();
  gLongComputationCancelled = 0;
  gLongComputationInProgress = 0;
}

