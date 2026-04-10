#include "kernel.h"

static void	initialize_matrix_from_angles( Tetrahedron *, GL4RMatrix );
static void	transpose( GL4RMatrix );
static Boolean	eigsrt( double *, GL4RMatrix );
static void	normalize( double *, GL4RMatrix);
static void	tred2(GL4RMatrix a, int n, double d[], double e[]);
static void	tqli(double d[], double e[], int n, GL4RMatrix z);

#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
#define EULER_EPSILON 0.001 



extern Boolean realize_tetrahedron_from_angles( Tetrahedron *tet )
{
 int		i,
		j;
 GL4RMatrix	g,g1;
 double		length1,
		length2,
		d[4],
		e[4];

 initialize_matrix_from_angles( tet, g );

 gl4R_invert(g,g1);

 /* check topology */

 for(i=0;i<4;i++)
 if (ABS(tet->cusp[i]->orbifold_euler_characteristic)< EULER_EPSILON)
 {
        if (ABS(g1[i][i]) > 0.0001)
              return FALSE;
 }
 else if (g1[i][i]*tet->cusp[i]->orbifold_euler_characteristic > 0)
        return FALSE;

 tred2(g,4,d,e);

 tqli( d, e, 4, g);

 eigsrt(d,g);

 normalize(d,g);

 transpose(g);

 gl4R_invert( g, tet->dual_basis );

 transpose(g);

 for(i=0;i<4;i++)
 {
	g[i][0]= -g[i][0];
	/* these are normalised onto the unit sphere */

	length1 = 0;
	length2 = 0;

	for(j = 0; j < 4; j++)
	{
		length1 += g[i][j]*g[i][j];
		length2 += tet->dual_basis[i][j]*tet->dual_basis[i][j];
	}

	length1 = sqrt(length1);
	length2 = sqrt(length2);

	if (length1< 0.0001 || length2 < 0.0001 )
	{
		uFatalError("realize_tetrahedron","tetrahedra_realization");
	}

	for(j = 0; j < 4; j++)
	{
		tet->basis[i][j] = g[i][j] / length1;
		tet->dual_basis[i][j] = tet->dual_basis[i][j] /length2;
	}

 }

 for(i=0;i<4;i++)
 {
   if (tet->basis[i][0] < -1 / sqrt(2) + 0.001)
   for(j=0;j<4;j++) tet->basis[i][j] = -tet->basis[i][j];

   if (tet->basis[i][0] < -1 / sqrt(2) + 0.001)
   for(j=0;j<4;j++) tet->dual_basis[i][j] = -tet->dual_basis[i][j];
 }

     return TRUE;
}


extern Boolean realize_tetrahedron_from_Gram_matrix( Tetrahedron *tet )
{
 int            i,
                j,
                k,
                l;
 GL4RMatrix     g1,g2;
 double         d1[4],
                sqrt_d1[4],
                e1[4];
 EdgeClass *edge;
 Cusp      *cusp;

 for(i=0;i<4;i++)
   for(j=0;j<4;j++)
   if (i!=j)
   {
          edge = tet->edge_class[edge_between_vertices[i][j]];
          g1[i][j] = edge->inner_product[ultimate];
   }
   else
   {
          cusp = tet->cusp[i];
          g1[i][i] = cusp->inner_product[ultimate];
   }

 tred2(g1,4,d1,e1);
 tqli( d1, e1, 4, g1);

 eigsrt(d1,g1);

 for(i=0;i<4;i++)
      tet->eigenvalue[i] = d1[i];

 sqrt_d1[0] = sqrt(fabs(d1[0]));

 for(i=0;i<3;i++)
 if (tet->dihedral_angle[ultimate][i] > PI)
       sqrt_d1[i+1] = -sqrt(fabs(d1[i+1]));
 else
       sqrt_d1[i+1] =  sqrt(fabs(d1[i+1]));

 for(i=0;i<4;i++)
     for(j=0;j<4;j++)
             tet->basis[i][j] = g1[i][j]*sqrt_d1[j];

 /* now to find the dual basis */

 for(i=0;i<4;i++)
 {
       l = -1;

       for(k=0;k<4;k++) g2[0][k] = 0;

       for(l=0,j=1;l<4 && j<4;l++)
       if (l==i)
            continue;
       else
       {
            for(k=0;k<4;k++)
                g2[j][k] = tet->basis[l][k];
            j++;
       }

       for(j=0;j<4;j++)
                 tet->dual_basis[i][j] =
                     (j==0) ? -minor1(g2,0,j)
                            :  minor1(g2,0,j);

       if (o31_inner_product(tet->dual_basis[i],tet->basis[i]) *
                             tet->orientation_parameter[ultimate] > 0)
           for(j=0;j<4;j++)
                 tet->dual_basis[i][j] *= -1;
 }

  return TRUE;
}



static void initialize_matrix_from_angles( Tetrahedron *tet, GL4RMatrix Gram_matrix )
{

 int 	i,
	j;


 for( i=0; i<4; i++)
	for(j=0;j<4;j++)
 		Gram_matrix[i][j] =  (i==j)  ? 1:
 			(( tet->dihedral_angle[ultimate][ edge_between_faces[i][j] ] < 0 )?
 			-cosh( tet->dihedral_angle[ultimate][ edge_between_faces[i][j] ])
 			:-cos( tet->dihedral_angle[ultimate][ edge_between_faces[i][j] ] ));

}

static void transpose( GL4RMatrix a)
{
 O31Matrix	temp;
 int		i,
		j;


 for( i=0; i<4; i++)
 	for( j=0; j<4; j++)
		temp[i][j] = a[j][i];

 for( i=0; i<4; i++)
	for( j=0; j<4; j++)
		a[i][j] = temp[i][j];


}



/* moves columns of v so that the negative eigenvalue is the first element of d */
static Boolean	eigsrt( double *d, GL4RMatrix v)
{
	int	neg,
		index,
		i;
	double	temp;

	neg = 0;
	index = 0;

        for(i=0;i<4;i++)
        if (d[i] < 0 && ABS(d[i]) < 0.00001 )
            d[i] = 0;

	for(i=0;i<4;i++)
	if (d[i] <= d[index])
	{
		index = i;
		if (d[i]<0) neg++;
	}

	if (neg != 1 && d[index]==0)
		uFatalError("eigsrt","tetrahedron_realization");;

	for(i=0;i<4;i++)
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


static void	normalize( double *d, GL4RMatrix v)
{
	int	i,
		j;
	double	x;

	for(i=0;i<4;i++)
	{
		x = sqrt(fabs(d[i]));

		for(j=0;j<4;j++)
			v[j][i] = v[j][i] / x;
	}
}



static void	tred2(GL4RMatrix a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for (i=n-1;i>=1;i--) {
	    
		l=i-1;
		h=scale=0.0;
		
		if (l > 0) {
		    
			for (k=0;k<=l;k++)
				scale += fabs(a[i][k]);
			
			if (scale == 0.0) 
				e[i]=a[i][l];
			else {
			    
				for (k=0;k<=l;k++) {
					a[i][k] /= scale;
					h += a[i][k]*a[i][k];
				}
				
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; 
				a[i][l]=f-g;
				f=0.0;
				
				for (j=0;j<=l;j++) {
					a[j][i]=a[i][j]/h; 
					g=0.0;
				
					for (k=0;k<=j;k++)
						g += a[j][k]*a[i][k];
					
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h;
					f += e[j]*a[i][j];
				}

				hh=f/(h+h);

				for (j=0;j<=l;j++) {
				    
				    f=a[i][j];
				    e[j]=g=e[j]-hh*f;

				    for (k=0;k<=j;k++)
					a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}

	d[0]=0.0;
	e[0]=0.0;

	for (i=0;i<n;i++) {
		l=i-1;

		if (d[i]) {

		    	for (j=0;j<=l;j++) {
				g=0.0;

				for (k=0;k<=l;k++)
					g += a[i][k]*a[k][j];

				for (k=0;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		
		d[i]=a[i][i];
		a[i][i]=1.0;

		for (j=0;j<=l;j++)
		    a[j][i]=a[i][j]=0.0;
	}
}


static void tqli(double d[], double e[], int n, GL4RMatrix z)
{

	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	
	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;

    	for (l=0;l<n;l++) {
		iter=0;
		
		do {
			for (m=l;m<n-1;m++) {

				dd=fabs(d[m])+fabs(d[m+1]);

				if (fabs(e[m])+dd == dd) break; /* float */
			}

		    	if (m != l) {
				if (iter++ == 30) uFatalError("tqli", "tetrahedra_realization");

				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=sqrt(g*g + 1);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;

				for (i=m-1;i>=l;i--) { /* not sure about this */

					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=sqrt(f*f + g*g));

					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}

					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;

								  
					for (k=0;k<n;k++) {
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}

				if (r == 0.0 && i >= 0) continue; /* 1 */
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}


