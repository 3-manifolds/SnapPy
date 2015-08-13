/*
  Test the change that John Berge requested so that
  relator from peripheral fillings comes from a 
  geometric curve on the cusp.  

  Since the elements in the fundamental group coming
  from the meridian and longitude commute, the key
  thing is just to make sure we're adding the correct
  number of them. 

  This version deals with the case (m,l) != 1 as well.
 */

#include "stdio.h"
#define ABS(x)  (((x) >= 0) ? (x) : -(x))
#define TRUE 1

void append_copies(int *curr, int change, int ignored);
void test_new_code(int m, int l);
long int gcd(long int a, long int b);

/* from kernel_code/gcd.c */
long int gcd(
    long int    a,
    long int    b)
{
    a = ABS(a);
    b = ABS(b);
    
    if (a == 0)
    {
        if (b == 0)
	  printf("Error");
        else
            return b;
    }

    while (TRUE)
    {
        if ((b = b%a) == 0)
            return a;
        if ((a = a%b) == 0)
            return b;
    }
}

void append_copies(int *curr, int change, int ignored){
    *curr = *curr + change;
}

void test_new_code(int m, int l){
    int meridian_value = 0;
    int longitude_value = 0;
    int new_word = 0;
    int *meridian = &meridian_value;
    int *longitude = &longitude_value;

    /* code from fundamental_group.c */
    int M=ABS(m), L=ABS(l);
    int m_sign=(M == m ? 1 : -1), l_sign=(L == l ? 1 : -1);
    int i=0, d = gcd(M,L);
    if (d > 1) {
      M /= d;
      L /= d;
    }
    while (d > 0) {
      do {
	if (i < M) append_copies(meridian, m_sign, new_word);
	else append_copies(longitude, l_sign, new_word);
	i += M;
	if (i >= M+L ) i -= (M+L);
      } while (i != 0);
      d -= 1;
    }
    if ( m != meridian_value || l != longitude_value)
	printf("(%d %d) -> (%d %d)\n", m, l, meridian_value, longitude_value);
}

int main(int argv, char* argc){
    int m, l, N = 100;
    for (m = -N; m <= N; m++){
	for (l = -N; l <= N; l++){
	  if (m != 0 || l != 0) 
	    test_new_code(m, l);
	}
    }
}
