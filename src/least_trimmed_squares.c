#include "least_trimmed_squares.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Least trimmed squares estimates for univariate location and variance.
   Jussi Tohka jussi.tohka@cs.tut.fi April 29th 2002

Input   : data is a vector containing n samples. Data samples are 
          expected to be truly real valued (i.e too many samples having
          the same value might lead to problems. 
Output  : mean and variance of the data.
The algorithm (exact) is described in
P.J. Rousseeuw and A.M. Leroy: Robust Regression and Outlier Detection
John Wiley & Sons 1987.

*/

int least_trimmed_squares(elem_type data[], long int n, 
                     elem_type* mean , elem_type* variance)
{
  int i,j,h,h2;
  double score,best_score,loc,best_loc,old_sum,new_sum,medd;
  double old_power_sum,new_power_sum;
  elem_type* scaled_data;

  h = n - n/2;
  h2 = n/2;

  qsort(data, n, sizeof(elem_type),compare_reals); 

  
  old_sum = 0;
  old_power_sum = 0.0;
  for(i = 0;i < h;i++) {
    old_sum = old_sum + data[i];
    old_power_sum = old_power_sum + data[i]*data[i];
  }

  loc = old_sum/h;  
  /* For better understanding of the algorithm: 
    O(N^2) implementation of the algorithm would compute score as:
    score = 0.0;
    for(i = 0;i < h;i++) {
      score = score + (data[i] - loc)*(data[i] - loc);
    } 
    But there is a faster way to this: */
  
  score = old_power_sum - old_sum*loc; 

  best_score = score;
  best_loc = loc;
  
  for(j = 1;j < h2 + 1;j++) {
    new_sum = old_sum - data[j - 1] + data[h - 1 + j];
    old_sum = new_sum;
    loc = old_sum/h;
    new_power_sum = old_power_sum - data[j - 1]*data[j - 1] 
                  + data[h - 1 + j]*data[h - 1 + j];
    old_power_sum = new_power_sum;
    score = old_power_sum - old_sum*loc; 

    if(score < best_score) {
      best_score = score;
      best_loc = loc;
    }
  }  
  *mean = best_loc;

  /* For the variance, it is needed to calculate the ellipsoid covering one half of samples. 
     This is not implemented optimally here because data has already been sorted. */

  scaled_data = malloc(n*sizeof(elem_type));
  if(scaled_data == NULL) return(1);
  for(i = 0; i < n ;i++) {
    scaled_data[i] = (data[i] - best_loc)*((h - 1)/best_score)*(data[i] - best_loc);
  }
  medd = median(scaled_data,n);
  free(scaled_data);
  *variance = (best_score/(h - 1))*(medd/CHI2INV_1);
  return(0);
}

int compare_reals (const elem_type *a, const elem_type *b)
{
 if( *a > *b)
    return(1);
  else
   return(-1);
}


/*
 * Algorithm from N. Wirth's book, implementation by N. Devillard.
 * This code in public domain.
 */

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array
   Notice   :   use the median() macro defined below to get the median. 

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/


elem_type kth_smallest(elem_type a[], long int n, long int k)
{
    long int i,j,l,m ;
    elem_type x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[k] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<k) l=i ;
        if (k<i) m=j ;
    }
    return a[k] ;
}



