/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include "minvarellipsoid.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Minimum variance ellipsoid estimator for the mean and variance.
   Jussi Tohka jussi.tohka@cs.tut.fi April 12th 2002

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

int minvarellipsoid(elem_type data[], long int n, int trials, 
                     elem_type* mean , elem_type* variance)
{
  int i,j,random_index[2];
  double score,best_score,sample_mean,sample_var,median_scale;
  double best_median_scale;
  elem_type* scaled_data;

 
  scaled_data = malloc(n * sizeof(elem_type));
  if(scaled_data == NULL) return(1);   
  srand(time(0));

  best_score = VERY_LARGE;
  for(i = 0;i < trials;i++) {
    random_index[0] = rint(rand()*((double) (n - 1) /RAND_MAX));
    random_index[1] = rint(rand()*((double) (n - 1) /RAND_MAX));
    
    sample_mean = 0.5*(data[random_index[0]] + data[random_index[1]]);
    sample_var = pow(data[random_index[0]] - sample_mean,2) +
               pow(data[random_index[1]] - sample_mean,2);
  
    if(fabs(sample_var) > 1E-16) {
      for(j = 0; j < n;j++) {   
        scaled_data[j] = (data[j] - sample_mean)*(1 / sample_var) *
                        (data[j] - sample_mean);
      }
      median_scale = median(scaled_data, n);     
      score = sqrt(median_scale*sample_var);
    }
    else {
      score = VERY_LARGE;
    }
    if(score < best_score) {
      *mean = sample_mean;
      *variance = sample_var;
      best_median_scale = median_scale;
      best_score = score;

    }
  }
  *variance = best_median_scale*(*variance)*(1/CHI2INV_1);
  free(scaled_data);      
  return(0);
}

