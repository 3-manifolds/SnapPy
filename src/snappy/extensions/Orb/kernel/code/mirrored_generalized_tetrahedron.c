#include "kernel.h"


extern int mirrored_generalized_tetrahedron( O31Matrix *gen, GL4RMatrix basis, GL4RMatrix dual_basis )
{
        int             num_gens,
                        index;

        num_gens = 0;

        for( index=0;index<4;index++)
                compute_reflection( index,gen[num_gens++], dual_basis);

        for( index=0;index<4;index++)
                if (ABS(basis[index][0]) < (1 / sqrt(2) - 0.001))
                        compute_reflection( index,gen[num_gens++], basis);

        return num_gens;
}


extern void compute_reflection( int index, O31Matrix gen, GL4RMatrix basis )
{

        int     i,
                j;
        double  length;

        length = o31_inner_product( basis[index], basis[index] );

        for(i=0;i<4;i++) for(j=0;j<4;j++)
		gen[i][j] = -2 * basis[index][i] * basis[index][j] / length;

	for(i=0;i<4;i++)
		gen[i][0] = -gen[i][0];

	for(i=0;i<4;i++)
		gen[i][i] = 1 + gen[i][i]; 

}
