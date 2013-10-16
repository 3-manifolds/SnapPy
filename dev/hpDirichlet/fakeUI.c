/*
 * For testing purposes we provide definitions of functions that the
 * SnapPea kernel expects its UI to provide.
 */
#include <stdio.h>
#include "SnapPea.h"
#include "Dirichlet.h"

  void uFatalError(const char *function, const char *file) {
    printf("Fatal error in function %s in file %s\n", function, file);
  }

  void uAcknowledge(const char *message){
    printf("%s\n", message);
  }
  
  int uQuery(const char *message, const int num_responses,
	     const char *responses[], const int default_response){
    return default_response;
  }

  void uLongComputationBegins(char *message, Boolean is_abortable){
    return;
  }

  FuncResult uLongComputationContinues(void){
    return func_OK;
  }

  void uLongComputationEnds(void){
    return;
  }

  void uAbortMemoryFull(void){
    return;
  }

  /* void precise_generators(MatrixPairList *gen_list){ */
  /*     return; */
  /* } */

  /* void precise_o31_product(O31Matrix a, O31Matrix b, O31Matrix product){ */
  /*     return; */
  /* } */
