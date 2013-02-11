/*
  Test the change that John Berge requested so that
  relator from peripheral fillings comes from a 
  geometric curve on the cusp.  

  Since the elements in the fundamental group coming
  from the meridian and longitude commute, the key
  thing is just to make sure we're adding the correct
  number of them. 

 */

#include "stdio.h"

void append_copies(int *curr, int change, int ignored);
void test_new_code(int m, int l);


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

    int M=(m<0 ? -m : m), L=(l<0 ? -l : l);
    int m_sign=(M == m ? 1 : -1), l_sign=(L == l ? 1 : -1);
    int i=0;
    do {
	if (i < M) append_copies(meridian, m_sign, new_word);
	else append_copies(longitude, l_sign, new_word);
	i += M;
	if (i >= M+L ) i -= (M+L);
    } while (i != 0);

    if ( m != meridian_value || l != longitude_value)
	printf("(%d %d) -> (%d %d)\n", m, l, meridian_value, longitude_value);
}

int main(int argv, char* argc){
    int m, l, N = 5;
    for (m = -N; m <= N; m++){
	for (l = -N; l <= N; l++){
	    test_new_code(m, l);
	}
    }
}
