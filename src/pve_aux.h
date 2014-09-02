/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* HEADER FOR THE AUXILIARY FUNCTIONS NEEDED BY pve.c 
   Take a look at .c file for comments on functions. */

#include "pve_common.h"

// #define ESTIMATOR "MCD"   /* MVE, MCD and ML supported */
#define ESTIMATOR "MVE"   /* MVE, MCD and ML supported */
// #define ESTIMATOR "ML"   /* MVE, MCD and ML supported */

/* Error messages */

#define ERROR_TOO_FEW    "Too few input argumants. Try pve -help \n \n"
#define ERROR_PARAMS_FILE "Problem with parameters file." 
#define ERROR_INPUT_FILE "Unable to open the input file or the brainmask file."
#define ERROR_PARAMS_IMAGE "Unable to estimate parameters based on the segmented image."
#define ERROR_PARAMS_TAG "Unable to estimate parameters based on tag points."
#define ERROR_SECOND_ARG "First argument must -help, -file or -image. \n \n"

void Display_help(); 
int Get_params_from_file(char* fn, double* mean, double* var, double* pvmeasurement);
int Estimate_params_from_image(Volume volume_in, Volume volume_mask, Volume volume_subcort,
                               Volume volume_seg, double* mean, 
                               double* var , double* pvmeasurement);

int Estimate_ml( double * samples, long int nofsamples, int stencil,
                 double* mean, double* var );
int Estimate_mve( double * samples, long int nofsamples, int stencil,
                  double* mean, double* var );
int Estimate_mcd( double * samples, long int nofsamples, int stencil,
                  double* mean, double* var );

double* Collect_values(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
		       long int* pcount, int neighbourhood);
double* Collect_values_subcortical(Volume volume_in,Volume volume_subcort, long int* pcount);

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

double Compute_Gaussian_likelihood(double value, double mean , double var);
double Compute_Gaussian_log_likelihood(double value, double mean , double var);
double Compute_marginalized_likelihood(double value, double mean1 , double mean2, 
                                       double var1, double var2, 
                                       double measurement_var);
void Parameter_estimation(Volume volume_in, Volume volume_mask, 
                         Volume probabilities[CLASSES],double* mean, double* var,
                         double* var_measurement);
int Parameter_estimation_classified(Volume volume_in, Volume volume_mask, 
                                    Volume volume_subcort, Volume classified, 
                                    double* mean, double* var, double * var_measurement );
int Compute_final_classification(Volume volume_in, Volume volume_classified,
                                 Volume final_cls, double * mean, double * var);
int Compute_partial_volume_vectors(Volume volume_in,Volume volume_classified,
                                   Volume volume_pve[PURE_CLASSES], double* mean);
int Compute_partial_volume_vectors_ml(Volume volume_in, Volume volume_classified,
                                      Volume volume_pve_ml[PURE_CLASSES],
                                      double* mean, double* var, double var_measurement);
void Usage_info(char* pname);
