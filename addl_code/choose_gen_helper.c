#include "kernel.h"
#include "kernel_namespace.h"

void choose_gen_tetrahedron_info(Triangulation    *manifold, 
				 int tet_index, 
				 int *generator_path,
				 int *face0_gen,
				 int *face1_gen,
				 int *face2_gen,
				 int *face3_gen,
				 Complex *corner0, 
				 Complex *corner1, 
				 Complex *corner2, 
				 Complex *corner3){
  Tetrahedron   *tet;
  FaceIndex f;
  int *gens[4];

  // choose_generators(manifold, TRUE, FALSE);
  for (tet = manifold->tet_list_begin.next; tet->index != tet_index; tet = tet->next){}
  
  *corner0 = tet->corner[0];
  *corner1 = tet->corner[1];
  *corner2 = tet->corner[2];
  *corner3 = tet->corner[3];
  *generator_path = tet->generator_path;
  
  gens[0] = face0_gen;
  gens[1] = face1_gen;
  gens[2] = face2_gen;
  gens[3] = face3_gen;

  for (f = 0; f < 4; f++){
    if (tet->generator_status[f] == outbound_generator)
      *gens[f] = tet->generator_index[f] + 1;
    if (tet->generator_status[f] == inbound_generator)
      *gens[f] = -(tet->generator_index[f] + 1);
    if (tet->generator_status[f] == not_a_generator)
      *gens[f] = 0;
  }
}
  
#include "end_namespace.h"
