/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* Header for MCD estimator in 1 D. */
/* Jussi Tohka, jussi.tohka@cs.tut.fi 29th April 2002.  */

#include "pve_support.h"
#include "pve_common.h"

extern int least_trimmed_squares(elem_type data[], long int n, int stencil,
				 elem_type* mean , elem_type* variance);
extern int compare_reals (const void *a, const void *b);

