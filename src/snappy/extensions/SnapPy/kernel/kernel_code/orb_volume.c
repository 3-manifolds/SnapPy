/**
 *  @file orb_volume.c
 *
 *  Ported from snappea/code/my_volume.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_volume.c
 */

#include "kernel.h"

#include "dilog.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static Complex U(Complex z, Real *angles);
static Real    tetrahedron_volume(Real *angles);

/*
 *  Ported from my_volume in snappea/code/my_volume.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_volume.c#L9-L41
 */

Real orb_volume(
    Triangulation *manifold)
{
    Real volume = 0.0;

    for (Tetrahedron *tet = manifold->tet_list_begin.next;
         tet != &manifold->tet_list_end;
         tet = tet->next)
    {
        Real tet_vol = tetrahedron_volume(
            tet->orb_tet_shape->dihedral_angle[ultimate]);

        if (tet->orb_tet_shape->orientation_parameter[ultimate] > 0)
            volume += tet_vol;
        else
            volume -= tet_vol;
    }

    return volume;
}

/*  Volume computed using formula in "A volume forumla for generalized hyperbolic tetrahedra" by Ushijima
 *
 *  Ported from tetrahedron_volume in snappea/code/my_volume.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_volume.c#L43-L119
 */
static Real tetrahedron_volume(
    Real angles[6])
{
    GL4RMatrix G;
    Complex w1, w2, w, z1, z2, bottom;
    Real sqrt_det, real_top;

    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            G[i][j] = (i == j) ? (Real)1.0 : -cos(angles[edge_between_faces[i][j]]);

    /* calculate the top of the complex numbers z1 and z2 */
    real_top = 0.0;
    for (int i = 0; i < 3; i++)
        real_top -= 2.0 * sin(angles[i]) * sin(angles[5 - i]);

    sqrt_det = sqrt(ABS(gl4R_determinant(G)));

    z1.real = real_top;
    z1.imag = 2.0 * sqrt_det;

    z2.real = real_top;
    z2.imag = -2.0 * sqrt_det;

    /* now for the bottom */
    bottom = Zero;

    for (int i = 0; i < 3; i++)
    {
        w1.real = cos(angles[i]);
        w1.imag = sin(angles[i]);
        w2.real = cos(angles[5 - i]);
        w2.imag = sin(angles[5 - i]);
        w = complex_mult(w1, w2);
        bottom = complex_plus(bottom, w);
    }

    for (int i = 0; i < 4; i++)
    {
        w = One;
        for (int j = 0; j < 4; j++)
            if (i != j)
            {
                w1.real = cos(angles[edge_between_faces[i][j]]);
                w1.imag = sin(angles[edge_between_faces[i][j]]);
                w = complex_mult(w, w1);
            }
        bottom = complex_plus(bottom, w);
    }

    w = One;
    for (int i = 0; i < 6; i++)
    {
        w1.real = cos(angles[i]);
        w1.imag = sin(angles[i]);
        w = complex_mult(w, w1);
    }

    bottom = complex_plus(bottom, w);

    z1 = complex_div(z1, bottom);
    z2 = complex_div(z2, bottom);

    return complex_minus(U(z1, angles), U(z2, angles)).imag / 2;
}

/*
 *  Ported from U in snappea/code/my_volume.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/my_volume.c#L122-L171
 */
static Complex U(
    Complex  z,
    Real    *angles)
{
    Complex result = complex_volume_dilog(z), w, w1, w2, dilogw;

    for (int i = 0; i < 3; i++)
    {
        w = One;

        for (int j = 0; j < 3; j++)
            if (i != j)
            {
                w1.real = cos(angles[j]);
                w1.imag = sin(angles[j]);
                w2.real = cos(angles[5 - j]);
                w2.imag = sin(angles[5 - j]);
                w = complex_mult(w, w1);
                w = complex_mult(w, w2);
            }

        w = complex_mult(w, z);
        dilogw = complex_volume_dilog(w);
        result = complex_plus(result, dilogw);
    }

    for (int i = 0; i < 4; i++)
    {
        w = MinusOne;

        for (int j = 0; j < 4; j++)
            if (i != j)
            {
                w1.real = cos(angles[edge_between_vertices[i][j]]);
                w1.imag = sin(angles[edge_between_vertices[i][j]]);
                w = complex_mult(w, w1);
            }

        w = complex_mult(w, z);
        dilogw = complex_volume_dilog(w);
        result = complex_minus(result, dilogw);
    }

    return complex_real_mult(0.5, result);
}

SNAPPEA_NAMESPACE_END_SCOPE
