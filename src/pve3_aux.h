/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* HEADER FOR THE AUXILIARY FUNCTIONS NEEDED BY pve3.c 
   Take a look at .c file for comments on functions. */

#include "pve_common.h"
#include "matrix3.h"
#include "minvarellipsoid3.h"

#define ESTIMATOR "ML"   /* MVE and ML supported. The usage of MVE estimator in the 
                            multispectral case is , however, not encouraged. 
                            This is since the implementation is slow, the results from 
                            the randomized algorithm vary too much and anyway the results 
                            obtained with ML estimator are as good. I suggest the usage 
                            of Matlab to obtain MCD estimates, or if that is not convinient 
                            then ML estimator is probably better choice than MVE. Note that MVE 
                            in the single spectral is still a good choice. */


/* Error messages */

#define ERROR_TOO_FEW    "Too few input argumants. Try pve3 -help \n \n"
#define ERROR_PARAMS_FILE "Problem with parameters file." 
#define ERROR_INPUT_FILE "Unable to open the input file or the brainmask file."
#define ERROR_PARAMS_IMAGE "Unable to estimate parameters based on the segmented image."
#define ERROR_PARAMS_TAG "Unable to estimate parameters based on tag points."
#define ERROR_SECOND_ARG "First argument must -help, -file or -cl. \n \n"

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
				Volume volume_seg, pVector means[PURE_CLASSES + 1],
				pMatrix vars[PURE_CLASSES + 1] , pMatrix pvmeasurement);

int Estimate_ml3( double * samples, long int nofsamples,
                  pVector mean, pMatrix var);

int Estimate_mve3( double * samples, long int nofsamples,
                   pVector mean,pMatrix var);

double* Collect_values3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
                        Volume volume_mask,Volume volume_seg,char ref_label,
                        long int* pcount, int neighbourhood);

double* Collect_values3_subcortical(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
				    Volume volume_subcort, long int* pcount);

int Estimate_params_from_tags3(char* tag_filename,Volume volume_inT1, Volume volume_inT2, Volume volume_inPD,
                                 pMatrix means[PURE_CLASSES + 1], pMatrix vars[PURE_CLASSES + 1]);  
int Open_images3(char* inT1_fn, char* inT2_fn, char* inPD_fn, 
                 char* mask_fn, Volume* pvolume_inT1, 
                 Volume* pvolume_inT2,Volume* pvolume_inPD,
                 Volume* pvolume_mask);

double Compute_Gaussian_likelihood3(pVector value, pVector mean , pMatrix var);
double Compute_Gaussian_likelihood3_fast(pVector value, pVector mean, pMatrix invvar,
                                         double sq_det );
double Compute_marginalized_likelihood3(pVector value, pVector mean1 , pVector mean2, 
                                       pMatrix var1, pMatrix var2, 
                                       pMatrix measurement_var);
void Parameter_estimation3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD, Volume volume_mask, 
                         Volume probabilities[CLASSES],pVector* mean, pMatrix* var,
                         pMatrix var_measurement);
int Parameter_estimation_classified3(Volume volume_inT1, Volume volume_inT2, Volume volume_inPD,
                                     Volume volume_mask, Volume volume_subcort, Volume classified,
                                     pVector mean[4], pMatrix var[4], pMatrix var_measurement);

int Compute_final_classification3(Volume volume_inT1, Volume volume_inT2, 
                                  Volume volume_inPD, Volume volume_classified,
                                  Volume final_cls,
                                  pVector means[PURE_CLASSES + 1], 
                                  pMatrix vars[PURE_CLASSES + 1] );
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

