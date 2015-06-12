/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "minvarellipsoid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Minimum volume ellipsoid estimator for the mean and variance.
   Jussi Tohka jussi.tohka@cs.tut.fi April 12th 2002
   Butchered by Claude Lepage for one-sided approximations. July 2014.

Input   : data is a vector containing n samples. Trials is the number 
          of trials that the algorithm uses, more trials more exact 
          the result is but also more time is spent. Data samples are 
          expected to be truly real valued (i.e too many samples having
          the same value might lead to problems. 
Output  : mean and variance of the data.
Algorithm (approximative) is described from page 259 onwards in
P.J. Rousseeuw and A.M. Leroy: Robust Regression and Outlier Detection
John Wiley & Sons 1987.

*/

int minvarellipsoid(elem_type data[], long int n, int trials, int stencil,
                    elem_type* mean , elem_type* variance) {

  long int i, j;
  double sample_mean, median_scale1, median_scale2;
  double best_median_scale1, best_median_scale2;
  double test_variance, minimum_variance;

  elem_type min_value = data[0];
  elem_type max_value = data[0];
  for( i = 1; i < n; i++ ) {
    if( data[i] < min_value ) min_value = data[i];
    if( data[i] > max_value ) max_value = data[i];
  }
  elem_type delta = ( max_value - min_value ) / (double)trials;

  minimum_variance = VERY_LARGE;
  for(i = 0;i < trials;i++) {
   
    sample_mean = min_value + (double)i * delta;

    median_scale1 = 0.0;
    median_scale2 = 0.0;
    if( stencil == EST_CENTERED ) {
      for(j = 0; j < n;j++) {   
        median_scale1 += (data[j] - sample_mean) * (data[j] - sample_mean);
      }
      median_scale1 = sqrt( median_scale1 / n );
      test_variance = median_scale1;
    } else {

      // Notes on Split Normal Distribution:
      // https://en.wikipedia.org/wiki/Split_normal_distribution
      // global mean = mu + sqrt(2/pi)*(scale2-scale1)
      // global mode = mu
      // global variance = (1-2/pi)*(scale2-scale1)**2 + scale1*scale2
      // There is another way proposed in the above link to compute
      // the one-sided variances median_scale1 and median_scale2
      // based on maximum likelihood. 

      for(j = 0; j < n; j++) {   
        if( data[j] < sample_mean ) {
          median_scale1 += (data[j] - sample_mean) * (data[j] - sample_mean);
        }
        if( data[j] > sample_mean ) {
          median_scale2 += (data[j] - sample_mean) * (data[j] - sample_mean);
        }
      }

      double term1 = pow( median_scale1, 1.0/3.0 );
      double term2 = pow( median_scale2, 1.0/3.0 );
      double term3 = sqrt( ( term1 + term2 ) / (double)n );

      test_variance = term1 + term2;

      median_scale1 = term1 * term3;
      median_scale2 = term2 * term3;

    }

    if( test_variance < minimum_variance ) {
      *mean = sample_mean;
      minimum_variance = test_variance;
      best_median_scale1 = median_scale1;
      best_median_scale2 = median_scale2;
    }
  }

  if( stencil == EST_CENTERED ) {
    *variance = best_median_scale1;
  } else if( stencil == EST_BELOW ) {
    *variance = best_median_scale1;
  } else if( stencil == EST_ABOVE ) {
    *variance = best_median_scale2;
  }
  *variance = (*variance)*(*variance);

  return(0);
}

