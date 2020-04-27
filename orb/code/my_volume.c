#include "kernel.h"
#if 1
#include "dilog.h"
#else
#include "gsl_sf_dilog.h"
#endif

static Complex Li_2( Complex z, Boolean *ok );
static Complex U( Complex z, double *angles, Boolean *ok );

/* Volume computed using formula in "A volume forumla for generalized hyperbolic tetrahedra" by Ushijima */

extern double my_volume( Triangulation *manifold, Boolean *ok )
{
  Tetrahedron *tet;
  double       volume=0,angles[6], tet_vol=0;
  int i;
  Boolean ok1;

 *ok = TRUE;

  for( tet = manifold->tet_list_begin.next;
       tet!=&manifold->tet_list_end;
       tet = tet->next )
  {
       for(i=0;i<6;i++) angles[i] = tet->dihedral_angle[ultimate][i];

	tet_vol = tetrahedron_volume( angles , &ok1 );

	if (!ok1)
	{
		if (!flat_tet(tet))
		{	
			*ok = FALSE;
			uFatalError("my_volume", "my_volume");
		}
	}

       if (tet->orientation_parameter[ultimate] > 0 ) 
		volume += tet_vol;
       else	volume -= tet_vol;
  }

  return volume;
}

extern double tetrahedron_volume( double *angles, Boolean *ok )
{
int i,j;
 static const int opposite[]={5,4,3};
Complex w1,w2,w,z1,z2,bottom;
GL4RMatrix G;
double sqrt_det,real_top;

*ok = TRUE;

for(i=0;i<4;i++)
   for(j=0;j<4;j++)
      G[i][j] = (i==j) ? 1 : -cos(angles[edge_between_faces[i][j]]); 

/* calculate the top of the complex numbers z1 and z2 */

real_top = 0;

for(i=0;i<3;i++)
   real_top -= 2 * sin(angles[i]) * sin(angles[opposite[i]]);

sqrt_det = sqrt( ABS( gl4R_determinant( G )));

z1.real = real_top;
z1.imag =  2*sqrt_det;

z2.real = real_top;
z2.imag = -2*sqrt_det;

/* now for the bottom */

bottom = Zero;

for(i=0;i<3;i++)
{
   w1.real = cos(angles[i]);
   w1.imag = sin(angles[i]);
   w2.real = cos(angles[opposite[i]]);
   w2.imag = sin(angles[opposite[i]]);

   w  = complex_mult( w1, w2 );

   bottom = complex_plus( bottom, w); 
}

for(i=0;i<4;i++)
{
   w = One;

   for(j=0;j<4;j++)
   if (i!=j)
   {
       w1.real = cos(angles[edge_between_faces[i][j]]);
       w1.imag = sin(angles[edge_between_faces[i][j]]);
       w = complex_mult( w, w1);
   }

   bottom = complex_plus( bottom , w);
}

w = One;

for(i=0;i<6;i++)
{
   w1.real = cos(angles[i]);
   w1.imag = sin(angles[i]);
   w = complex_mult( w, w1);
}

bottom = complex_plus(bottom, w);

z1 = complex_div( z1, bottom );
z2 = complex_div( z2, bottom );

return complex_minus( U(z1,angles, ok), U(z2,angles, ok)).imag / 2;

}


static Complex U( Complex z, double *angles, Boolean *ok)
{
  int i,j;
static const int opposite[]={5,4,3};
  Complex result = Li_2(z, ok), w, w1,w2, dilogw;

  for(i=0;i<3;i++)
  {
      w = One;

      for(j=0;j<3;j++)
      if (i!=j)
      {
           w1.real = cos(angles[j]);
           w1.imag = sin(angles[j]);
	   w2.real = cos(angles[opposite[j]]);
	   w2.imag = sin(angles[opposite[j]]);

           w  = complex_mult( w, w1 );
	   w = complex_mult( w, w2 );
      }

      w = complex_mult( w, z );
      dilogw = Li_2( w, ok );

      result = complex_plus( result, dilogw );
  }

  for(i=0;i<4;i++)
  {
      w = MinusOne;

      for(j=0;j<4;j++)
      if (i!=j)
      {
          w1.real = cos(angles[edge_between_vertices[i][j]]);
          w1.imag = sin(angles[edge_between_vertices[i][j]]);
          w  = complex_mult( w, w1 );
      }

      w = complex_mult( w, z );
      dilogw = Li_2( w, ok );

      result = complex_minus( result, dilogw );
  }

  result = complex_real_mult( .5, result);

  return result;
}

#if 1
static Complex Li_2( Complex z, Boolean *ok )
{
    *ok = TRUE;

    return complex_volume_dilog(z);
}

#else

static Complex Li_2( Complex z, Boolean *ok )
{

  Complex w;
  gsl_sf_result re, im; 

  gsl_sf_complex_dilog_e( complex_modulus(z), atan2(z.imag,  z.real), &re, &im);

  if (re.val != re.val )
  {
	w.real = 0;
	*ok = FALSE;
  }
  else  w.imag = re.val;   	/* ok yes. this looks suspect.  the gsl dilog functions is returning */
				/* nan in some answers so this is a way around in the mean time */
  if (im.val != im.val )	/* from what i've seen the nan are occuring when the answer should be zero */
  {
	w.imag = 0;
	*ok = FALSE;
  }
  else  w.imag = im.val;

  return w;
}

#endif



