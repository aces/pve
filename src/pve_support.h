/* PVE_SUPPORT.H
 */
#ifndef _PVE_SUPPORT_H_
#define _PVE_SUPPORT_H_ 1

typedef double elem_type;

#define CHI2INV_1 0.45493642311957
#define CHI2INV_3 2.36597388437534 

#define VERY_LARGE 1E16

#define ELEM_SWAP(a,b) { register elem_type t=(a);(a)=(b);(b)=t; }

#define median(a,n) kth_smallest(a,n,(((n)&1)?((n)/2):(((n)/2)-1)))

extern elem_type kth_smallest(elem_type a[], long int n, long int k);

#endif /* _PVE_SUPPORT_H_ */
