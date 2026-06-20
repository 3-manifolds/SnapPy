/**
 *  @file orb_tetrahedron_development.c
 *
 *  Ported from:
 *  snappea/code/tetrahedron_realization.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/tetrahedron_realization.c
 *  snappea/code/new_choose_generators.c
 *  https://github.com/DamianHeard/orb/blob/f1bbe9a2170b172278c6fa43bd8039dfd6a66276/snappea/code/new_choose_generators.c
 *
 *  The geometric data lives in Tetrahedron::orb_tet_shape,
 *  EdgeClass::orb_edge_shape and Cusp::orb_cusp_shape.
 */

#include "kernel.h"

SNAPPEA_NAMESPACE_BEGIN_SCOPE

static void initialize_matrix_from_angles(Tetrahedron *tet, GL4RMatrix Gram_matrix);
static void transpose(GL4RMatrix matrix);
static Boolean eigsrt(Real *d, GL4RMatrix v);
static void normalize(Real *d, GL4RMatrix v);
static void tred2(GL4RMatrix a, int n, Real d[], Real e[]);
static void tqli(Real d[], Real e[], int n, GL4RMatrix z);
static Boolean realize_tetrahedron_from_angles(Tetrahedron *tet);
static void compute_lorentz_transformation(
    GL4RMatrix basis,
    GL4RMatrix dual_basis,
    GL4RMatrix image_basis,
    GL4RMatrix image_dual_basis,
    FaceIndex face,
    Permutation gluing,
    GL4RMatrix transformation);

#define SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))
#define EULER_EPSILON 0.001

static Boolean realize_tetrahedron_from_angles(
    Tetrahedron *tet)
{
    int         i,
                j;
    GL4RMatrix  g,
                g_inverse;
    Real        length1,
                length2,
                d[4],
                e[4];

    initialize_matrix_from_angles(tet, g);

    gl4R_invert(g, g_inverse);

    /* check topology */

    for (i = 0; i < 4; i++)
    {
        Real orbifold_euler_characteristic =
            orb_compute_orbifold_cusp_euler_characteristic(tet->cusp[i]);

        if (ABS(orbifold_euler_characteristic) < EULER_EPSILON)
        {
            if (ABS(g_inverse[i][i]) > 0.0001)
                return FALSE;
        }
        else if (g_inverse[i][i] * orbifold_euler_characteristic > 0)
            return FALSE;
    }

    tred2(g, 4, d, e);
    tqli(d, e, 4, g);
    eigsrt(d, g);
    normalize(d, g);
    transpose(g);

    gl4R_invert(g, tet->orb_tet_shape->dual_basis);
    transpose(g);

    for (i = 0; i < 4; i++)
    {
        g[i][0] = -g[i][0];
        /* these are normalised onto the unit sphere */

        length1 = 0.0;
        length2 = 0.0;

        for (j = 0; j < 4; j++)
        {
            length1 += g[i][j] * g[i][j];
            length2 += tet->orb_tet_shape->dual_basis[i][j]
                     * tet->orb_tet_shape->dual_basis[i][j];
        }

        length1 = sqrt(length1);
        length2 = sqrt(length2);

        if (length1 < 0.0001 || length2 < 0.0001)
            uFatalError("realize_tetrahedron_from_angles",
                        "orb_tetrahedron_development");

        for (j = 0; j < 4; j++)
        {
            tet->orb_tet_shape->basis[i][j] = g[i][j] / length1;
            tet->orb_tet_shape->dual_basis[i][j] /=
                length2;
        }
    }

    for (i = 0; i < 4; i++)
    {
        if (tet->orb_tet_shape->basis[i][0] < -1 / sqrt(2) + 0.001)
            for (j = 0; j < 4; j++)
                tet->orb_tet_shape->basis[i][j] =
                    -tet->orb_tet_shape->basis[i][j];

        if (tet->orb_tet_shape->basis[i][0] < -1 / sqrt(2) + 0.001)
            for (j = 0; j < 4; j++)
                tet->orb_tet_shape->dual_basis[i][j] =
                    -tet->orb_tet_shape->dual_basis[i][j];
    }

    return TRUE;
}


Boolean orb_realize_tetrahedron_from_Gram_matrix(
    Tetrahedron *tet)
{
    int         i,
                j,
                k,
                l;
    GL4RMatrix  g,
                minor_matrix;
    Real        d[4],
                sqrt_d[4],
                e[4];
    EdgeClass   *edge;
    Cusp        *cusp;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            if (i != j)
            {
                edge = tet->edge_class[edge_between_vertices[i][j]];
                g[i][j] = edge->orb_edge_shape->inner_product[ultimate];
            }
            else
            {
                cusp = tet->cusp[i];
                g[i][i] = cusp->orb_cusp_shape->inner_product[ultimate];
            }

    tred2(g, 4, d, e);
    tqli(d, e, 4, g);
    eigsrt(d, g);

    for (i = 0; i < 4; i++)
        tet->orb_tet_shape->eigenvalue[i] = d[i];

    sqrt_d[0] = sqrt(fabs(d[0]));

    for (i = 0; i < 3; i++)
        if (tet->orb_tet_shape->dihedral_angle[ultimate][i] > PI)
            sqrt_d[i + 1] = -sqrt(fabs(d[i + 1]));
        else
            sqrt_d[i + 1] =  sqrt(fabs(d[i + 1]));

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            tet->orb_tet_shape->basis[i][j] = g[i][j] * sqrt_d[j];

    /* now to find the dual basis */

    for (i = 0; i < 4; i++)
    {
        for (k = 0; k < 4; k++)
            minor_matrix[0][k] = 0.0;

        for (l = 0, j = 1; l < 4 && j < 4; l++)
            if (l != i)
            {
                for (k = 0; k < 4; k++)
                    minor_matrix[j][k] =
                        tet->orb_tet_shape->basis[l][k];
                j++;
            }

        for (j = 0; j < 4; j++)
            tet->orb_tet_shape->dual_basis[i][j] =
                (j == 0) ? -orb_minor1(minor_matrix, 0, j)
                         :  orb_minor1(minor_matrix, 0, j);

        if (o31_inner_product(tet->orb_tet_shape->dual_basis[i],
                              tet->orb_tet_shape->basis[i])
                * tet->orb_tet_shape->orientation_parameter[ultimate] > 0)
            for (j = 0; j < 4; j++)
                tet->orb_tet_shape->dual_basis[i][j] *= -1;
    }

    return TRUE;
}


static void initialize_matrix_from_angles(
    Tetrahedron *tet,
    GL4RMatrix  Gram_matrix)
{
    int i,
        j;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            Gram_matrix[i][j] = (i == j) ? (Real)1.0 :
                ((tet->orb_tet_shape->dihedral_angle[ultimate]
                    [edge_between_faces[i][j]] < 0) ?
                    -cosh(tet->orb_tet_shape->dihedral_angle[ultimate]
                        [edge_between_faces[i][j]]) :
                    -cos(tet->orb_tet_shape->dihedral_angle[ultimate]
                        [edge_between_faces[i][j]]));
}


static void transpose(
    GL4RMatrix matrix)
{
    O31Matrix temp;
    int       i,
              j;

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            temp[i][j] = matrix[j][i];

    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            matrix[i][j] = temp[i][j];
}

/* moves columns of v so that the negative eigenvalue is the first element of d */
static Boolean eigsrt(
    Real        *d,
    GL4RMatrix  v)
{
    int     neg,
            index,
            i;
    Real    temp;

    neg = 0;
    index = 0;

    for (i = 0; i < 4; i++)
        if (d[i] < 0 && ABS(d[i]) < 0.00001)
            d[i] = 0.0;

    for (i = 0; i < 4; i++)
        if (d[i] <= d[index])
        {
            index = i;
            if (d[i] < 0)
                neg++;
        }

    if (neg != 1 && d[index] == 0)
        uFatalError("eigsrt", "orb_tetrahedron_development");

    for (i = 0; i < 4; i++)
    {
        temp = v[i][index];
        v[i][index] = v[i][0];
        v[i][0] = temp;
    }

    temp = d[index];
    d[index] = d[0];
    d[0] = temp;

    return TRUE;
}


static void normalize(
    Real        *d,
    GL4RMatrix  v)
{
    int     i,
            j;
    Real    x;

    for (i = 0; i < 4; i++)
    {
        x = sqrt(fabs(d[i]));

        for (j = 0; j < 4; j++)
            v[j][i] = v[j][i] / x;
    }
}


static void tred2(
    GL4RMatrix  a,
    int         n,
    Real        d[],
    Real        e[])
{
    int     l,
            k,
            j,
            i;
    Real    scale,
            hh,
            h,
            g,
            f;

    for (i = n - 1; i >= 1; i--)
    {
        l = i - 1;
        h = scale = 0.0;

        if (l > 0)
        {
            for (k = 0; k <= l; k++)
                scale += fabs(a[i][k]);

            if (scale == 0.0)
                e[i] = a[i][l];
            else
            {
                for (k = 0; k <= l; k++)
                {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }

                f = a[i][l];
                g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                f = 0.0;

                for (j = 0; j <= l; j++)
                {
                    a[j][i] = a[i][j] / h;
                    g = 0.0;

                    for (k = 0; k <= j; k++)
                        g += a[j][k] * a[i][k];

                    for (k = j + 1; k <= l; k++)
                        g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }

                hh = f / (h + h);

                for (j = 0; j <= l; j++)
                {
                    f = a[i][j];
                    e[j] = g = e[j] - hh * f;

                    for (k = 0; k <= j; k++)
                        a[j][k] -= (f * e[k] + g * a[i][k]);
                }
            }
        }
        else
            e[i] = a[i][l];
        d[i] = h;
    }

    d[0] = 0.0;
    e[0] = 0.0;

    for (i = 0; i < n; i++)
    {
        l = i - 1;

        if (d[i] != 0.0)
        {
            for (j = 0; j <= l; j++)
            {
                g = 0.0;

                for (k = 0; k <= l; k++)
                    g += a[i][k] * a[k][j];

                for (k = 0; k <= l; k++)
                    a[k][j] -= g * a[k][i];
            }
        }

        d[i] = a[i][i];
        a[i][i] = 1.0;

        for (j = 0; j <= l; j++)
            a[j][i] = a[i][j] = 0.0;
    }
}


static void tqli(
    Real        d[],
    Real        e[],
    int         n,
    GL4RMatrix  z)
{
    int     m,
            l,
            iter,
            i,
            k;
    Real    s,
            r,
            p,
            g,
            f,
            dd,
            c,
            b;

    for (i = 1; i < n; i++)
        e[i - 1] = e[i];
    e[n - 1] = 0.0;

    for (l = 0; l < n; l++)
    {
        iter = 0;

        do
        {
            for (m = l; m < n - 1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m + 1]);

                if (fabs(e[m]) + dd == dd) /* float */
                    break;
            }

            if (m != l)
            {
                if (iter++ == 30)
                    uFatalError("tqli", "orb_tetrahedron_development");

                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt(g * g + 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--) /* not sure about this */
                {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = (r = sqrt(f * f + g * g));

                    if (r == 0.0)
                    {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;

                    for (k = 0; k < n; k++)
                    {
                        f = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }

                if (r == 0.0 && i >= 0) /* 1 */
                    continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }
        while (m != l);
    }
}


static void compute_lorentz_transformation(
    GL4RMatrix  basis,
    GL4RMatrix  dual_basis,
    GL4RMatrix  image_basis,
    GL4RMatrix  image_dual_basis,
    FaceIndex   face,
    Permutation gluing,
    GL4RMatrix  transformation)
{
    GL4RMatrix  m1,
                m2,
                m1_inverse;
    int         i,
                j,
                sign;
    Real        length1 = 1.0,
                length2 = 1.0,
                length3 = 1.0,
                length4 = 1.0;

    sign = (parity[gluing] == orientation_preserving) ? -1 : 1;

    for (i = 0; i < 4; i++)
    {
        if (i != face)
        {
            length1 = sqrt(ABS(o31_inner_product(basis[i], basis[i])));
            length2 = sqrt(ABS(o31_inner_product(
                image_basis[EVALUATE(gluing, i)],
                image_basis[EVALUATE(gluing, i)])));
        }
        else
        {
            length3 = sqrt(ABS(o31_inner_product(
                dual_basis[i], dual_basis[i])));
            length4 = sqrt(ABS(o31_inner_product(
                image_dual_basis[EVALUATE(gluing, i)],
                image_dual_basis[EVALUATE(gluing, i)])));
        }

        if (length1 < 0.0001)
            length1 = 1.0;
        if (length2 < 0.0001)
            length2 = 1.0;
        if (length3 < 0.0001)
            length3 = 1.0;
        if (length4 < 0.0001)
            length4 = 1.0;

        for (j = 0; j < 4; j++)
        {
            m1[i][j] = (i == face) ?
                dual_basis[i][j] / length3 :
                basis[i][j] / length1;
            m2[i][j] = (i == face) ?
                sign * image_dual_basis[EVALUATE(gluing, i)][j] / length4 :
                image_basis[EVALUATE(gluing, i)][j] / length2;
        }
    }

    gl4R_invert(m1, m1_inverse);
    o31_product(m1_inverse, m2, transformation);
}


void orb_compute_corners_of_neighbor(
    Tetrahedron *tet,
    FaceIndex   face)
{
    Tetrahedron *nbr_tet = tet->neighbor[face];
    Permutation gluing = tet->gluing[face];
    FaceIndex   nbr_face = EVALUATE(gluing, face);
    GL4RMatrix  transformation;

    if (orb_realize_tetrahedron_from_Gram_matrix(nbr_tet) == FALSE)
        uFatalError("orb_compute_corners_of_neighbor",
                    "orb_tetrahedron_development.c");

    compute_lorentz_transformation(
        nbr_tet->orb_tet_shape->basis,
        nbr_tet->orb_tet_shape->dual_basis,
        tet->orb_tet_shape->basis,
        tet->orb_tet_shape->dual_basis,
        nbr_face,
        inverse_permutation[gluing],
        transformation);

    o31_product(
        nbr_tet->orb_tet_shape->basis,
        transformation,
        nbr_tet->orb_tet_shape->basis);
    o31_product(
        nbr_tet->orb_tet_shape->dual_basis,
        transformation,
        nbr_tet->orb_tet_shape->dual_basis);
}

SNAPPEA_NAMESPACE_END_SCOPE
