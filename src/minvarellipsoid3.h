/* Header for minimum variance ellipse estimator in 3 D. */
/* Jussi Tohka, jussi.tohka@cs.tut.fi 12th April 2002.  */
#include "pve_support.h"

extern int minvarellipsoid3(elem_type data[], long int n, int trials, 
			    pVector mean , pMatrix  variance);
 
