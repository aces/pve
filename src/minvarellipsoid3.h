/* Header for minimum variance ellipse estimator in 3 D. */
/* Jussi Tohka, jussi.tohka@cs.tut.fi 12th April 2002.  */
#include "matrix3.h"

typedef double elem_type ;
#define CHI2INV_3 2.36597388437534 
#define VERY_LARGE 1E16

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type kth_smallest(elem_type a[], long int n, long int k);

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

int minvarellipsoid3(elem_type data[], long int n, int trials, 
                     pVector mean , pMatrix  variance);
 
