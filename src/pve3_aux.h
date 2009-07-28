/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* HEADER FOR THE AUXILIARY FUNCTIONS NEEDED BY pve3.c 
   Take a look at .c file for comments on functions. */


#include <volume_io.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix3.h"
#include "minvarellipsoid3.h"

#define DATATYPE_SIZE 4096 /* i.e 12 bits */
#define MASK_TR  0.5           /* Values below this constant in the volume_mask
                                 indicate that corresponding voxels are of non brain
                                 tissue. */
#define SC_TR  1.5             /* Values below this constant in the volume_subcort
                                 indicate that corresponding voxels are not subcortical. */

#define MAX_ITERATIONS 20      /* Maximum number of iterations of ICM */
#define LIKELIHOOD_RANGE_MIN 0.0
#define LIKELIHOOD_RANGE_MAX 1.0

#define MAXSAMPLES 50000 
#define MAXTRIALS 10000
#define DEFAULT_TAG_FILENAME "/usr/local/mni/data/classify/ntags_1000_prob_90_nobg.tag"
#define NEIGHBOURHOOD 6     /* These are for image based parameter estimation */
#define ESTIMATOR "ML"   /* MVE and ML supported. The usage of MVE estimator in the 
                             multispectral case is , however, not encouraged. 
                             This is since the implementation is slow, the results from 
                             the randomized algorithm vary too much and anyway the results 
                             obtained with ML estimator are as good. I suggest the usage 
                             of Matlab to obtain MCD estimates, or if that is not convinient 
                             then ML estimator is probably better choice than MVE. Note that MVE 
                             in the single spectral is still a good choice. */

#define MEASUREMENT_FRACTION 0.5 /* if the model with both measurent noise and 
                                    physiological noise is used, this defines how big
                                    part the measurement noise is from the less varying pure 
                                    tissue type */


#define PURE_CLASSES 4   /* The number of pure tissue classes */
#define MIXED_CLASSES 5  /* The number of mixed tissue classes */
#define CLASSES 9        /* The number of mixed and pure tissue classes */

/* Labels: It is assumed that:
   BGLABEL = 0;
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
/* Error messages */

#define ERROR_TOO_FEW    "Too few input argumants. Try pve3 -help \n \n"
#define ERROR_PARAMS_FILE "Problem with parameters file." 
#define ERROR_INPUT_FILE "Unable to open the input file or the brainmask file."
#define ERROR_PARAMS_IMAGE "Unable to estimate parameters based on the segmented image."
#define ERROR_PARAMS_TAG "Unable to estimate parameters based on tag points."
#define ERROR_SECOND_ARG "First argument must -help, -file or -cl. \n \n"

#define NOML 0
#define MLONLY 1
#define ML 2

#define BETA 0
#define SAME 1
#define SMLR 2
#define DIFF 3

/* Look up table for Potts model: Tells which classes are similar. Look at the macro 
Are_similar to figure out how this works. */
extern const char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1];
#define Are_same(label1,label2) label1 == label2
#define Are_similar(label1,label2) POTTS_LOOKUP_TABLE[label1][label2]


void Display_help(); 
int Are_proper(pMatrix var_wm,pMatrix var_gm,pMatrix var_csf, 
               pMatrix var_bg,pMatrix var_measurement);
int Get_params_from_file3(char* fn, pVector means[PURE_CLASSES + 1],
                          pMatrix vars[PURE_CLASSES + 1], pMatrix pvmeasurement);
/* int Get_params_from_cl3(char** inputstr, double* pmwm, double* pmgm, double* pmcsf, 
                         double* pmbg, double* pvwm, double* pvgm, double* pvcsf , 
                         double* pvbg , double* pvmeasurement, double* pbeta); */
int Estimate_params_from_image3(Volume volume_inT1, Volume volume_inT2, Volume volume_inPD, 
				Volume volume_mask, Volume volume_subcort, 
				char* segmentation_filename, pVector means[PURE_CLASSES + 1],
				pMatrix vars[PURE_CLASSES + 1] , pMatrix pvmeasurement);

int Estimate_ml3(Volume volume_inT1, Volume volume_inT2,Volume volume_inPD,
		 Volume volume_mask,Volume volume_seg, Volume volume_subcort,
		 char ref_label,pVector mean,pMatrix var);

int Estimate_mve3(Volume volume_in, Volume volume_inT2,Volume volume_inPD,
                  Volume volume_mask,Volume volume_seg, Volume volume_subcort,
                  char ref_label, pVector mean,pMatrix var);

double* Collect_values3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
                        Volume volume_mask,Volume volume_seg,char ref_label,
                        long int* pcount);

double* Collect_values3_subcortical(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
				    Volume volume_subcort,char ref_label, long int* pcount);

int Estimate_params_from_tags3(char* tag_filename,Volume volume_inT1, Volume volume_inT2, Volume volume_inPD,
                                 pMatrix means[PURE_CLASSES + 1], pMatrix vars[PURE_CLASSES + 1]);  
int Open_images3(char* inT1_fn, char* inT2_fn, char* inPD_fn, 
                 char* mask_fn, Volume* pvolume_inT1, 
                 Volume* pvolume_inT2,Volume* pvolume_inPD,
                 Volume* pvolume_mask);


/* char Are_same(char label1,char label2);
   char Are_similar(char label1, char label2); */
void Limit_0_1(double* x);

int Normalize(double* pval,char  n);
char Maxarg(double* pval, char  n);

double Compute_Gaussian_likelihood3(pVector value, pVector mean , pMatrix var);
double Compute_Gaussian_likelihood3_fast(pVector value, pVector mean, pMatrix invvar,
                                         double sq_det );
double Compute_marginalized_likelihood3(pVector value, pVector mean1 , pVector mean2, 
                                       pMatrix var1, pMatrix var2, 
                                       pMatrix measurement_var);
void Compute_mrf_probability(double * mrf_prob, Volume* pvolume, int x, int y , int z, 
                             double* width_stencil, double beta, double same, double similar, 
                             double different, double prior);
void Parameter_estimation3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD, Volume volume_mask, 
                         Volume probabilities[CLASSES],pVector* mean, pMatrix* var,
                         pMatrix var_measurement);
int Parameter_estimation_classified3(Volume volume_inT1, Volume volume_inT2, Volume volume_inPD,
                                     Volume volume_mask, Volume volume_subcort, Volume classified,
                                     pVector mean[4], pMatrix var[4], pMatrix var_measurement);

int Compute_partial_volume_vectors3(Volume volume_inT1, Volume volume_inT2, 
                                    Volume volume_inPD, Volume volume_classified,
                                    Volume volume_pve[PURE_CLASSES],
                                    pVector means[PURE_CLASSES + 1], 
                                    pMatrix vars[PURE_CLASSES + 1], pMatrix var_measurement );
int Compute_partial_volume_vectors3_ml(Volume volume_inT1, Volume volume_inT2, 
                                    Volume volume_inPD, Volume volume_classified,
                                    Volume volume_pve[PURE_CLASSES],
                                    pVector means[PURE_CLASSES + 1], 
                                    pMatrix vars[PURE_CLASSES + 1], pMatrix var_measurement);

void Usage_info(char* pname);

