/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* pve_common.c: auxiliary functions for pve.c and pve3.c. */

#include "pve_common.h"

/* Function that defines the sampling voxels based on the label
   and the masked classification for a given neighbourhood. These 
   labels must be identical in the multi-spectral case. */

long int Collect_neighbourhood( Volume volume_ngh, Volume volume_restrict,
                                Volume volume_seg, char ref_label,
                                int neighbourhood ) {

  int i, j, k, i1, j1, k1;
  int sizes[MAX_DIMENSIONS];
  long int count = 0;
  BOOLEAN on_the_border;

  get_volume_sizes( volume_seg, sizes );

  for(i = 2; i < sizes[0] - 2; i++) {
    for(j = 2;j < sizes[1] - 2; j++) {
      for( k = 2; k < sizes[2] - 2;k++) {
        set_volume_real_value(volume_ngh,i,j,k,0,0,0);
        short mask_val = 1;
        if( volume_restrict ) {
          mask_val = floor( get_volume_real_value(volume_restrict,i,j,k,0,0) + 0.50 );
        }
        if( mask_val > MASK_TR ) {     /* Voxel in the brain */
          if(rint(get_volume_real_value(volume_seg,i,j,k,0,0)) == ref_label) {  /* And of right type */
            /* Check neigbourhood: Choices for neighbourhood are 6 and 26 neignbourhood or
               no neighbourhood checking.
               The choice is defined by the constant NEIGHBOURHOOD  */
            on_the_border = FALSE;

            if(neighbourhood == 125) { /* 125 - neighbourhood -- good for 0.5mm voxels */
              for(i1 = -2;i1 <= 2;i1++) {
                for(j1 = -2;j1 <= 2;j1++) {
                  for(k1 = -2;k1 <= 2;k1++) {
                    if( volume_restrict ) {
                      if( rint( get_volume_real_value(volume_restrict,i+i1,j+j1,k+k1,0,0) ) < MASK_TR ) break;
                    }
                    if(rint(get_volume_real_value(volume_seg,i + i1,j + j1,k + k1,0,0))
                        != ref_label) {
                      on_the_border = TRUE;
                      break;
                    }
                  }
                  if( on_the_border ) break;
                }
                if( on_the_border ) break;
              }
            } else if(neighbourhood == 26) { /* 26 - neighbourhood */
              for(i1 = -1;i1 < 2;i1++) {
                for(j1 = -1;j1 < 2;j1++) {
                  for(k1 = -1;k1 < 2;k1++) {
                    if( volume_restrict ) {
                      if( rint( get_volume_real_value(volume_restrict,i+i1,j+j1,k+k1,0,0) ) < MASK_TR ) break;
                    }
                    if(rint(get_volume_real_value(volume_seg,i + i1,j + j1,k + k1,0,0))
                        != ref_label) {
                      on_the_border = TRUE;
                      break;
                    }
                  }
                  if( on_the_border ) break;
                }
                if( on_the_border ) break;
              }
            } else if(neighbourhood == 19) { /* 19 - neighbourhood */
              for(i1 = -1;i1 < 2;i1++) {
                for(j1 = -1;j1 < 2;j1++) {
                  for(k1 = -1;k1 < 2;k1++) {
                    if( volume_restrict ) {
                      if( rint( get_volume_real_value(volume_restrict,i+i1,j+j1,k+k1,0,0) ) < MASK_TR ) break;
                    }
                    if( abs(i1)+abs(j1)+abs(k1) <= 2 ) {
                      if(rint(get_volume_real_value(volume_seg,i + i1,j + j1,k + k1,0,0))
                          != ref_label) {
                        on_the_border = TRUE;
                        break;
                      }
                    }
                  }
                  if( on_the_border ) break;
                }
                if( on_the_border ) break;
              }

            } else if(neighbourhood == 6) { /* 6-neighbourhood */
              for(i1 = -1;i1 < 2;i1++) {
                if( volume_restrict ) {
                  if( rint( get_volume_real_value(volume_restrict,i+i1,j,k,0,0) ) < MASK_TR ) break;
                }
                if(rint(get_volume_real_value(volume_seg,i + i1,j,k,0,0))
                        != ref_label) {
                  on_the_border = TRUE;
                  break;
                }
              }
              for(i1 = -1;i1 < 2;i1++) {
                if( volume_restrict ) {
                  if( rint( get_volume_real_value(volume_restrict,i,j+i1,k,0,0) ) < MASK_TR ) break;
                }
                if(rint(get_volume_real_value(volume_seg,i,j + i1,k,0,0))
                        != ref_label) {
                  on_the_border = TRUE;
                  break;
                }
              }
              for(i1 = -1;i1 < 2;i1++) {
                if( volume_restrict ) {
                  if( rint( get_volume_real_value(volume_restrict,i,j,k+i1,0,0) ) < MASK_TR ) break;
                }
                if(rint(get_volume_real_value(volume_seg,i,j,k + i1,0,0))
                        != ref_label) {
                  on_the_border = TRUE;
                  break;
                }
              }
            }
            if(!on_the_border) {
              count = count + 1;
              set_volume_real_value(volume_ngh,i,j,k,0,0,1);
            }
          }
        }
      }
    }
  }

  return( count );
}

/* Function that defines the sampling voxels for sub-cortical 
   tissue type. We ignore the neighbourhood stencil for now. */

long int Collect_neighbourhood_subcortical( Volume volume_ngh, 
                                            Volume volume_restrict,
                                            Volume volume_subcort,
                                            Volume volume_seg, 
                                            int neighbourhood ) {

  int i, j, k;
  int sizes[MAX_DIMENSIONS];

  long int count = 0;
  get_volume_sizes( volume_subcort, sizes );

  // We should consider the "neighbourhood" here below.

  for(i = 1; i < sizes[0] - 1;i++) {
    for(j = 1;j < sizes[1] - 1;j++) {
      for( k = 1; k < sizes[2] - 1;k++) {
        set_volume_real_value(volume_ngh,i,j,k,0,0,0);
        short mask_val = 1;
        if( volume_restrict ) {
          mask_val = floor( get_volume_real_value(volume_restrict,i,j,k,0,0) + 0.50 );
        }
        if( mask_val > MASK_TR ) {     /* Voxel in the brain */
          char c = rint(get_volume_real_value(volume_seg,i,j,k,0,0));
          if( ( c == GMLABEL || c == SCLABEL ) &&
              get_volume_real_value(volume_subcort,i,j,k,0,0) > SC_TR ) {
            count++;
            set_volume_real_value(volume_ngh,i,j,k,0,0,1);
          }
        }
      }
    }
  }
  return( count );
}


/* Function that collects intensity values for parameter estimation 
   into one vector. Returns the pointer to the vector where the samples 
   are collected. This function supports also MVE estimator */

double * Collect_values( Volume volume_in, Volume volume_ngh,
                         long int * pcount, char * estimator ) {

  int sizes[MAX_DIMENSIONS];
  int i, j,k;
  long int count = 0;
  long int sampleno;
  double *sample_vector_tmp;
  double *sample_vector;

  get_volume_sizes( volume_in, sizes );

  count = *pcount;
  if(count > 0) {

    sample_vector = malloc(count * sizeof (double));
    if(sample_vector == NULL) return(NULL);
    sampleno = 0;
    for(i = 2; i < sizes[0] - 2; i++) {
      for(j = 2; j < sizes[1] - 2; j++) {
        for( k = 2; k < sizes[2] - 2; k++) {
          if(get_volume_real_value(volume_ngh,i,j,k,0,0) > 0.5) {
            sample_vector[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
            sampleno++;
          }
        }
      }
    }
    count = sampleno;

    /* Computers today are fast enough to handle all samples, so
       keep them all. This will have an effect if volume is a
       1mm or 0.5mm. CL. */

    /* if we are using MVE or MCD estimator so
     1) we must get rid of some samples;
     2) we must add little bit of noise to the other samples; */
    if( strcmp(estimator,"MVE")==0 || strcmp(estimator,"MCD")==0 ) {

      if( count > MAXSAMPLES ) {
        // srand(time(0)); /* Set seed */
        srand(123456); /* Set seed */
        sample_vector_tmp = Downsample_values( sample_vector, count, &sampleno );
        free( sample_vector );
        sample_vector = sample_vector_tmp;
        count = sampleno;
      }

      // We don't want noise for real data, unless data come from a simulator. CL.
      // double voxel_value_interval = ( get_volume_real_max(volume_in) - 
      //                                 get_volume_real_min(volume_in) ) /
      //                               DATATYPE_SIZE;
      // AddNoise_values( sample_vector, count, voxel_value_interval );

    }
  } else {
    sample_vector = NULL;
  }
  *pcount = count;
  return(sample_vector);
}

/* Function that collects intensity values within the sub-cortical
   mask for parameter estimation. Returns the pointer to the vector 
   where the samples are collected.
   This function supports also MVE estimator */

double * Collect_values_subcortical( Volume volume_in, Volume volume_subcort,
                                     Volume volume_seg, long int * pcount,
                                     char * estimator ) {

  int sizes[MAX_DIMENSIONS];
  int i, j, k;
  long int count = 0;
  long int sampleno;
  double *sample_vector_tmp;
  double *sample_vector;

  get_volume_sizes(volume_subcort,sizes);


  if(count > 0) {

    sample_vector = malloc(count * sizeof (double));
    if(sample_vector == NULL) return(NULL);
    sampleno = 0;
    for(i = 1; i < sizes[0] - 1;i++) {
      for(j = 1;j < sizes[1] - 1;j++) {
        for( k = 1; k < sizes[2] - 1;k++) {
          char c = rint(get_volume_real_value(volume_seg,i,j,k,0,0));
          if( ( c == GMLABEL || c == SCLABEL ) &&
              get_volume_real_value(volume_subcort,i,j,k,0,0) > SC_TR ) {
            sample_vector[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
            sampleno = sampleno + 1;
          }
        }
      }
    }
    count = sampleno;

    /* if we are using MVE or MCD estimator so
     1) we must get rid of some samples;
     2) we must add little bit of noise to the other samples; */

    if( strcmp(estimator,"MVE")==0 || strcmp(estimator,"MCD")==0 ) {

      if( count > MAXSAMPLES ) {
        // srand(time(0)); /* Set seed */
        srand(123456); /* Set seed */
        sample_vector_tmp = Downsample_values( sample_vector, count, &sampleno );
        free( sample_vector );
        sample_vector = sample_vector_tmp;
        count = sampleno;
      }

      // double voxel_value_interval = ( get_volume_real_max(volume_in) - 
      //                                 get_volume_real_min(volume_in) ) /
      //                               DATATYPE_SIZE;
      // AddNoise_values( sample_vector, count, voxel_value_interval );

    }
  } else {
    sample_vector = NULL;
  }

  *pcount = count;
  return(sample_vector);
}

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

void Compute_mrf_probability(double* mrf_probability, Volume pvolume, int x, int y , int z, 
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
          label2 = get_volume_real_value(pvolume, x + i, y + j , z + k,0,0);
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


