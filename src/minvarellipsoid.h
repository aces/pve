/* Header for minimum variance ellipse estimator in 1 D. */
/* Jussi Tohka, jussi.tohka@cs.tut.fi 12th April 2002.  */

#include "pve_support.h"

extern int minvarellipsoid(elem_type data[], long int n, int trials, 
			   elem_type* mean , elem_type* variance);
 
