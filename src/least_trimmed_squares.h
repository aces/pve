/* Header for MCD estimator in 1 D. */
/* Jussi Tohka, jussi.tohka@cs.tut.fi 29th April 2002.  */

typedef double elem_type ;
#define CHI2INV_1 0.45493642311957
#define VERY_LARGE 1E16

int least_trimmed_squares(elem_type data[], long int n, 
                     elem_type* mean , elem_type* variance);
int compare_reals (const elem_type *a, const elem_type *b);

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

elem_type kth_smallest(elem_type a[], long int n, long int k);

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

