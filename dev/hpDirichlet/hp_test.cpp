#include <iostream>
#include "SnapPea.h"
#include "hp_Dirichlet.h"

using std::cout;
using std::endl;

#include "examples/K13a1779.h"

void qdprint(qd_real x);
void load_moebius(COMPLEX *coeffs, hp_MoebiusTransformation *M);
void print_moebius(hp_MoebiusTransformation *M);
void print_o31(hp_O31Matrix O);

extern hp_WEPolyhedron  *hp_Dirichlet_from_generators(
    hp_O31Matrix            generators[],
    int                     num_generators,
    REAL                    vertex_epsilon,
    DirichletInteractivity  interactivity,
    Boolean                 maximize_injectivity_radius);

int main() {
  char output[200];
  unsigned int oldcw;
  hp_MoebiusTransformation M[2];
  hp_O31Matrix generators[2];
  hp_WEPolyhedron* P;
  COMPLEX det;
  Boolean check;
  int i, j;

  fpu_fix_start(&oldcw);

  load_moebius(A_coeffs, &M[0]);
  load_moebius(B_coeffs, &M[1]);
  hp_Moebius_array_to_O31_array(M, generators, 2);

  /***
  cout << "Moebius transformation A:" << endl;
  print_moebius(&M[0]);
  det = hp_sl2c_determinant(M[0].matrix);
  cout << "determinant: " << det.real() << "+ i*" << det.imag() << endl;
  //  hp_Moebius_to_O31(&M[0], generators[0]);
  print_o31(generators[0]);

  cout << "Moebius transformation B:" << endl;
  print_moebius(&M[1]);
  det = hp_sl2c_determinant(M[1].matrix);
  cout << "determinant: " << det.real() << "+ i*" << det.imag() << endl;
  print_o31(generators[1]);

  check = hp_O31_determinants_OK(generators, 2, (REAL)"1.0e-55");
  if (check)
    cout << "determinants look OK" << endl;
  else
    cout << "determinants look bad" << endl;
  cout << "gen 0: ";
  qdprint(hp_gl4R_determinant(generators[0]));
  cout << endl << "gen 1: ";
  qdprint(hp_gl4R_determinant(generators[1]));
  cout << endl;
  **/

  P = hp_Dirichlet_from_generators(generators, 2, (REAL)"1.0e-24", 0, 1);
  cout << P->num_vertices << " vertices" << endl;
  cout << P->num_edges << " edges" << endl;
  cout << P->num_faces << " faces" << endl;
 
  fpu_fix_end(&oldcw);
  return 0;
}

void qdprint(qd_real x) {
  char buffer[100];
  x.write(buffer, 100, 68);
  cout << buffer << endl;
}

void load_moebius(COMPLEX *coeffs, hp_MoebiusTransformation *M) {
  M->matrix[0][0] = coeffs[0];
  M->matrix[0][1] = coeffs[1];
  M->matrix[1][0] = coeffs[2];
  M->matrix[1][1] = coeffs[3];
  M->parity = orientation_preserving;
}

void print_moebius(hp_MoebiusTransformation *M) {
  int i, j;
  for (i=0; i<2; i++)
    for (j=0;j<2;j++){
      qdprint(M->matrix[i][j].real());
      cout << " + i ";
      qdprint(M->matrix[i][j].imag());
      cout << endl;
    }
}

void print_o31(hp_O31Matrix O) {
  int i, j;
  for (i=0; i<4; i++)
    for (j=0;j<4;j++){
      qdprint(O[i][j]);
    }
}
