#include "matrix3.h"
#include "minvarellipsoid3.h"
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

int minvarellipsoid3(double data[], long int n, int trials, 
                     pVector mean , pMatrix variance)
{
  int i,j,k,random_index[4];
  double score,best_score;
  Vector3D sample_mean, standardized_data;
  Matrix3D sample_var,inv_sample_var,tmpmatrix;
  double median_scale;
  double best_median_scale;
  double* scaled_data;
  double d;
 
  scaled_data = malloc(n * sizeof(double)); 
  if(scaled_data == NULL) return(1); 
  srand(time(0));

  best_score = VERY_LARGE;
  for(i = 0;i < trials;i++) {
    for(k = 0;k < 4;k++) {
      random_index[k] = rint(rand()*((double) (n - 1) /RAND_MAX));
    }  
    for(k = 0;k < 3;k++) {
      sample_mean[k] = 0.25*(data[3*random_index[0] + k] + data[3*random_index[1] + k]
                     + data[3*random_index[2] + k] + data[3*random_index[3] + k]);
    }
    SetZero(sample_var); 
    for(k = 0;k < 4;k++) { 
      for(j = 0;j < 3;j++) {
        standardized_data[j] = data[3*random_index[k] + j] - sample_mean[j];
      }
      VectorProduct(standardized_data,standardized_data,tmpmatrix);
      AddMatrices(sample_var,tmpmatrix,sample_var);
    }
    ScalarMultiply(sample_var,(double) 1/3,sample_var);
    d = Determinant(sample_var);
    if(fabs(d) > 1E-16) {
      Invert(sample_var,inv_sample_var);
      for(j = 0; j < n;j++) { 
        standardized_data[0] = data[3*j]     - sample_mean[0];
        standardized_data[1] = data[3*j + 1] - sample_mean[1];
        standardized_data[2] = data[3*j + 2] - sample_mean[2];
        scaled_data[j] = QuadraticForm(inv_sample_var,standardized_data);
      }
      median_scale = median(scaled_data, n);     
      score = sqrt(d)*(pow(sqrt(median_scale),3));
    }
    else {
      score = VERY_LARGE;
    }
    if(score < best_score) {
      CopyVector(sample_mean,mean);
      CopyMatrix(sample_var,variance);
      best_median_scale = median_scale;
      best_score = score;

    }
  }
  ScalarMultiply(variance,best_median_scale/CHI2INV_3,variance);
  free(scaled_data);
  return(0);    
}

