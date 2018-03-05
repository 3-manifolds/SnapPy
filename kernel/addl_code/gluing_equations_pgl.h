/* Prototypes for functions defined in gluing_equations.c */

#ifndef __GLUING_EQUATIONS_PGL_H__
#define __GLUING_EQUATIONS_PGL_H__

#include "ptolemy_types.h"

/*

        Returns a matrix with rows of the form

                  a b c  d e f  ...

        which means

	    z0^a * z0'^b * z0''^c  *  z1^d * z1'^e * z1''^f  * ... = 1

        for an edge equation, and (same) for a cusp equation.

   	          z0' = 1 / (1-z0), z0'' = (z0 - 1) / z0

        In terms of the tetrahedra, z0 is the invariant of the edge
        (2,3), z0' the invariant of the edge (0,2) and z0'' is the
        invariant of the edge (1,2).  See kernel_code/edge_classes.c
        for a detailed account of the convention.

	See
	Garoufalidis, Zickert, Goerner: 
	Gluing Equations for PGL(n,C)-Representations of 3-Manifolds 
	http://arxiv.org/abs/1207.6711
	for detailed description.
*/

void get_edge_gluing_equations_pgl(
    Triangulation *manifold, Integer_matrix_with_explanations *m, int N);

void get_face_gluing_equations_pgl(
    Triangulation *manifold, Integer_matrix_with_explanations *m, int N);

void get_internal_gluing_equations_pgl(
    Triangulation *manifold, Integer_matrix_with_explanations *m, int N);

void get_cusp_equations_pgl(
    Triangulation *manifold, Integer_matrix_with_explanations *m, int N,
    int cusp_num, int meridians, int longitudes);

#endif
