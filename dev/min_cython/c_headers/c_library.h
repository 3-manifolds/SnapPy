typedef struct {
    int  entries[2][2]; 
}  GL2RMatrix;

extern int timestwo(int x);
extern int square(int x);

extern void multiply_GL2R(GL2RMatrix* A, GL2RMatrix* B, GL2RMatrix* C); 
