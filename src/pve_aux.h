/* HEADER FOR THE AUXILIARY FUNCTIONS NEEDED BY pve.c 
   Take a look at .c file for comments on functions. */

#include <volume_io.h>
#include <stdlib.h>
#include <stdio.h>


#define VERY_SMALL 1E-16
#define DATATYPE_SIZE 4096 /* i.e 12 bits */
#define MASK_TR  0.5           /* Values below this constant in the volume_mask
                                 indicate that corresponding voxels are of non brain
                                 tissue. */
#define INTERVALS 50           /* For numerical integration */
#define MAX_ITERATIONS 50      /* Maximum number of iterations of ICM */
#define LIKELIHOOD_RANGE_MIN 0.0
#define LIKELIHOOD_RANGE_MAX 1.0

#define MAXSAMPLES 50000 
#define MAXTRIALS 10000
#define DEFAULT_TAG_FILENAME "/usr/local/mni/data/classify/ntags_1000_prob_90_nobg.tag"
#define NEIGHBOURHOOD 6     /* These are for image based parameter estimation */
#define ESTIMATOR "MCD"   /* MVE, MCD and ML supported */

#define MEASUREMENT_FRACTION 0.5 /* if the model with both measurent noise and 
                                    psysiological noise is used, this defines how big
                                    part the measurement noise is from the less varying pure 
                                    tissue type */

#define PURE_CLASSES 3   /* The number of pure tissue classes */
#define MIXED_CLASSES 3  /* The number of mixed tissue classes */
#define CLASSES 6        /* The number of mixed and pure tissue classes */

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
#define GMCSFLABEL 4
#define WMGMLABEL 5
#define CSFBGLABEL 6

/* Error messages */

#define ERROR_TOO_FEW    "Too few input argumants. Try pve -help \n \n"
#define ERROR_PARAMS_FILE "Problem with parameters file." 
#define ERROR_INPUT_FILE "Unable to open the input file or the brainmask file."
#define ERROR_PARAMS_IMAGE "Unable to estimate parameters based on the segmented image."
#define ERROR_PARAMS_TAG "Unable to estimate parameters based on tag points."
#define ERROR_SECOND_ARG "First argument must -help, -file or -image. \n \n"

/* Look up table for Potts model: Tells which classes are similar. Look at the macro 
Are_similar to figure out how this works. */
const char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1] = {{0, 0, 0, 0, 0, 0, 1},
                                                           {0, 0, 0, 0, 1, 0, 1},
                                                           {0, 0, 0, 0, 1, 1, 0},
                                                           {0, 0, 0, 0, 0, 1, 0},
                                                           {0, 1, 1, 0, 0, 0, 0},
                                                           {0 ,0, 1, 1, 0, 0, 0},
                                                           {1, 1, 0, 0, 0, 0, 0}};

/* Some macros for convinience */

#define Are_same(label1,label2) label1 == label2
#define Are_similar(label1,label2) POTTS_LOOKUP_TABLE[label1][label2]


void Display_help(); 
int Get_params_from_file(char* fn, double* mean, double* var, double* pvmeasurement);
int Estimate_params_from_image(Volume volume_in, Volume volume_mask, 
                               char* segmentation_filename, double* mean, 
                               double* var , double* pvmeasurement);
int Estimate_ml(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var);
int Estimate_mve(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var);
int Estimate_mcd(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var);
double* Collect_values(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                long int* pcount);

/* double Estimate_mean(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label);
double Estimate_variance(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,double mean);
*/
int Estimate_params_from_tags(char* tag_filename,Volume volume_in,
                                 double* means, double* vars);
int Get_params_from_cl(char** inputstr, double* pmwm, double* pmgm, double* pmcsf, 
                         double* pmbg, double* pvwm, double* pvgm, double* pvcsf , 
                         double* pvbg , double* pvmeasurement); 
int Open_images(char* in_fn, char* mask_fn, Volume* pvolume_in, 
                Volume* pvolume_mask);

/* char Are_same(char label1,char label2) {return(label1 == label2); }
   char Are_similar(char label1, char label2) { return POTTS_LOOKUP_TABLE[label1][label2];} */
void Limit_0_1(double* x);

void Normalize(double* pval, char n);
char Maxarg(double* pval, char n);

double Compute_Gaussian_likelihood(double value, double mean , double var);
double Compute_marginalized_likelihood(double value, double mean1 , double mean2, 
                                       double var1, double var2, 
                                       double measurement_var,
                                       unsigned int nof_intervals);
double Compute_mrf_probability(char label, Volume* pvolume, int x, int y , int z, 
                               double* slice_width, double beta, int same, int similar, 
                               int different, double prior, int* sizes);
void Parameter_estimation(Volume volume_in, Volume volume_mask, 
                         Volume probabilities[CLASSES],double* mean, double* var,
                         double* var_measurement);
int Compute_partial_volume_vectors(Volume volume_in,Volume volume_classified,
                                   Volume volume_pve[PURE_CLASSES], double* mean);
int Compute_partial_volume_vectors_ml(Volume volume_in, Volume volume_classified,
                                      Volume volume_pve_ml[PURE_CLASSES],
                                      double* mean, double* var, double var_measurement);


