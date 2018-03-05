#include "c_library.h" 

void multiply_GL2R(GL2RMatrix* A, GL2RMatrix* B, GL2RMatrix* C){
    int i, j, k;
    for (i=0; i < 2; i++){
        for (j=0; j < 2; j++){
            for (k=0; k < 2; k++){
                C->entries[i][j] +=
		    (A->entries[i][k])*(B->entries[k][j]) +
		    (A->entries[i][k])*(B->entries[k][j]);
            }
        }
    }
}
