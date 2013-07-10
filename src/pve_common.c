/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* pve_common.c: auxilary functions for pve.c and pve3.c. */

#include "pve_common.h"

/* Downsample the samples. */

double * Downsample_values( double * samples, long int count, long int * new_count ) {

  int i;

  *new_count = 0;
  double * sample_vector = malloc(MAXSAMPLES * sizeof(double));
  if( sample_vector == NULL ) return(NULL);
  for(i =0;i < MAXSAMPLES;i++) {
    int sampleno = rint(rand() *( (double) (count - 1) / RAND_MAX));
    sample_vector[i] = samples[sampleno];
  }
  
  *new_count = MAXSAMPLES;
  return( sample_vector );
}

/* Add a bit of uniform random noise to the samples. */

void * AddNoise_values( double * samples, long int count, double voxel_value_interval ) {

  int i;

  for(i = 0; i < count; i++) {
    samples[i] = samples[i] + (rand() - RAND_MAX / 2)*  
                              voxel_value_interval/RAND_MAX;
  } 

}

/* -------------------------------------------------------------------
  Functions that do tiny tasks for other functions from this file.
  They are not called from main
  -------------------------------------------------------------------- */
/* No explanation needed ... */

/* Limits x to the range from 0 to 1. */ 

void Limit_0_1(double* x) 
{
  if( *x < 0) *x = 0.0;
  else if( *x > 1) *x = 1.0 ;

}

/* ---------------------------------------------------------------------------
   And finally functions that are called from main to do the actual computations.
  ---------------------------------------------------------------------------- */

/* Normalizes n probalities to sum to one */

int Normalize(double* pval, char n)

{ double sca = 0.0;
  char i;

  for(i = 0;i < n;i++) {
    sca = sca + pval[i];
  }
  if(fabs(sca) >  VERY_SMALL) {       /* To avoid divisions by zero */
    for(i = 0;i < n;i++) {
      pval[i] = pval[i]/sca;
    }
    return 0;
  } else {
    return 1;
  }
}


/* Finds maximum argument out of the n possibilities */

char Maxarg(double* pval,char n)
{
  double maximum,index;
  char i;
  
  maximum = pval[0];
  index = 1;
  for(i = 1;i < n;i++) {
    if(pval[i] > maximum) {
      index = i + 1;
      maximum = pval[i];
    }
  }
  return(index);
}


/* Computes prior probalibility of a voxel have a certain label given the neighbourhood labels.
   Takes classified volume and three coordinates.
    The first argument is of course the label under investigation.
   prior is the prior probability of the label to occur.
   Integers same , similar and different and double beta specify the MRF.
   No need to check for on_the_border as this is checked in the likelyhood loop.
 */

void Compute_mrf_probability(double* mrf_probability, Volume* pvolume, int x, int y , int z, 
                             double* width_stencil, double beta, double same, double similar, 
                             double different, double prior, BOOLEAN sc_region ) {

  int i,j,k,ii;
  char label, label2;  
  double similarity_value;

  for( label = 0; label < CLASSES; label++ ) {
    mrf_probability[label] = 0;
  }

  ii = 0; 
  for(i = -1; i < 2; i++) {
    for(j = -1; j < 2; j++) {
      for(k = -1; k < 2; k++) {
        if( i == 0 && j == 0 && k == 0 ) {
          for( label = 0; label < CLASSES; label++ ) {
            mrf_probability[label] += prior;
          }
        } else {
          label2 = get_volume_real_value(*pvolume, x + i, y + j , z + k,0,0);
          for( label = 0; label < CLASSES; label++ ) {
            if(Are_same(label+1,label2)) {
              similarity_value = same;
            } else if(Are_similar(label+1,label2)) {
              if (((label+1) == GMCSFLABEL)||((label+1) == CSFLABEL)) {
                similarity_value = similar;
              } else {
                similarity_value = -1;
              }
            } else {
              similarity_value = different;
            }
            mrf_probability[label] += similarity_value * width_stencil[ii];
          }
        }
        ii++;
      }
    }
  } 

  for( label = 0; label < CLASSES; label++ ) {

    if (sc_region) {
      if( label+1 == WMGMLABEL ) {
        mrf_probability[WMGMLABEL-1] = 0;
      } else {
        mrf_probability[label] = exp( -beta * mrf_probability[label] );
      }
    } else {
      if( label+1 == WMSCLABEL || label+1 == SCGMLABEL || label+1 == SCLABEL ) {
        mrf_probability[label] = 0;
      } else {
        mrf_probability[label] = exp( -beta * mrf_probability[label] );
      }
    }

  }
} 

void initialize_Potts_table() {

  int  i, j;

  for( i = 0; i <= CLASSES; i++ ) {
    for( j = 0; j <= CLASSES; j++ ) {
      POTTS_LOOKUP_TABLE[i][j] = 0;
    }
  }
  POTTS_LOOKUP_TABLE[BGLABEL][CSFBGLABEL] = 1;

  POTTS_LOOKUP_TABLE[CSFLABEL][GMCSFLABEL] = 1;
  POTTS_LOOKUP_TABLE[CSFLABEL][CSFBGLABEL] = 1;

  POTTS_LOOKUP_TABLE[GMLABEL][GMCSFLABEL] = 1;
  POTTS_LOOKUP_TABLE[GMLABEL][WMGMLABEL] = 1;
  POTTS_LOOKUP_TABLE[GMLABEL][SCGMLABEL] = 1;

  POTTS_LOOKUP_TABLE[WMLABEL][WMGMLABEL] = 1;
  POTTS_LOOKUP_TABLE[WMLABEL][WMSCLABEL] = 1;

  POTTS_LOOKUP_TABLE[SCLABEL][WMSCLABEL] = 1;
  POTTS_LOOKUP_TABLE[SCLABEL][SCGMLABEL] = 1;

  POTTS_LOOKUP_TABLE[GMCSFLABEL][CSFLABEL] = 1;
  POTTS_LOOKUP_TABLE[GMCSFLABEL][GMLABEL] = 1;

  POTTS_LOOKUP_TABLE[WMGMLABEL][GMLABEL] = 1;
  POTTS_LOOKUP_TABLE[WMGMLABEL][WMLABEL] = 1;

  POTTS_LOOKUP_TABLE[CSFBGLABEL][CSFLABEL] = 1;
  POTTS_LOOKUP_TABLE[CSFBGLABEL][BGLABEL] = 1;

  POTTS_LOOKUP_TABLE[WMSCLABEL][SCLABEL] = 1;
  POTTS_LOOKUP_TABLE[WMSCLABEL][WMLABEL] = 1;

  POTTS_LOOKUP_TABLE[SCGMLABEL][SCLABEL] = 1;
  POTTS_LOOKUP_TABLE[SCGMLABEL][GMLABEL] = 1;

  for( i = 0; i <= CLASSES; i++ ) {
    for( j = 0; j <= CLASSES; j++ ) {
      if( POTTS_LOOKUP_TABLE[i][j] > 0 ) {
        POTTS_LOOKUP_TABLE[i][j] = POTTS_LOOKUP_TABLE[j][i];
      }
    }
  }
}


