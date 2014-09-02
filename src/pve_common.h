/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

#include <volume_io.h>
#include <stdlib.h>
#include <stdio.h>

#define VERY_SMALL 1E-16
#define DATATYPE_SIZE 4096 /* i.e 12 bits */
#define MASK_TR  0.5           /* Values below this constant in the volume_mask
                                 indicate that corresponding voxels are of non brain
                                 tissue. */
#define MASK_CHANGED 1         /* Values equal to this constant in the volume_mask
                                 indicate that corresponding voxels has changed label. */
#define MASK_CHANGED_TR 1.5
#define MASK_NOTCHANGED 2      /* Values equal to this constant in the volume_mask
                                 indicate that corresponding voxels has not changed label. */

#define SC_TR  1.5             /* Values below this constant in the volume_subcort
                                 indicate that corresponding voxels are not subcortical. */

#define MAX_ITERATIONS 20      /* Maximum number of iterations of ICM */
#define LIKELIHOOD_RANGE_MIN 0.0
#define LIKELIHOOD_RANGE_MAX 1.0

#define MAXSAMPLES 200000
#define MAXTRIALS 1000
#define NEIGHBOURHOOD 6     /* These are for image based parameter estimation */

#define MEASUREMENT_FRACTION 0.5 /* if the model with both measurent noise and 
                                    physiological noise is used, this defines how big
                                    part the measurement noise is from the less varying pure 
                                    tissue type */

#define PURE_CLASSES 4   /* The number of pure tissue classes */
#define MIXED_CLASSES 5  /* The number of mixed tissue classes */
#define CLASSES 9        /* The number of mixed and pure tissue classes */

/* Labels: It is assumed that:   BGLABEL = 0;
   pure tissue labels are from 1 to 3;
   mixed tissue classes have labels 4 to 6;
   otherwise you should have no problems in renumbering classes.
   However, you must update also the Potts look up table and the partial 
   volume class look up table if you change the labels.
  */

#define WMLABEL 3
#define GMLABEL 2
#define CSFLABEL 1
#define BGLABEL 0
#define GMCSFLABEL 5
#define WMGMLABEL 6
#define CSFBGLABEL 7
#define SCLABEL 4
#define WMSCLABEL 8
#define SCGMLABEL 9

#define NOML 0
#define MLONLY 1
#define ML 2

#define BETA 0
#define SAME 1
#define SMLR 2
#define DIFF 3

#define EST_BELOW 1
#define EST_ABOVE 2
#define EST_CENTERED 4

/* Look up table for Potts model: Tells which classes are similar. Look at the macro 
   Are_similar to figure out how this works. */
extern char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1];
void initialize_Potts_table();

/* Some macros for convenience */

#define Are_same(label1,label2) label1 == label2
#define Are_similar(label1,label2) POTTS_LOOKUP_TABLE[label1][label2]
#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)

void Limit_0_1(double* x);
int Normalize(double* pval, char n);
char Maxarg(double* pval, char n);

double * Downsample_values( double *, long int, long int * );
void * AddNoise_values( double *, long int, double );
void Compute_mrf_probability(double * mrf_prob, Volume* pvolume, int x, int y , int z,
                             double* width_stencil, double beta, double same, double similar,
                             double different, double prior, BOOLEAN sc_region );
