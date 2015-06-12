/* ----------------------------- MNI Header -----------------------------------
@NAME       :  pve
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Partial volume clssification of the single channel MRI.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jan 2002 (Jussi Tohka, jupeto@cs.tut.fi) 
@MODIFIED   : Last modification in Aug 2004 by Vivek Singh
---------------------------------------------------------------------------- */

/*
   Copyright Alan C. Evans
   Professor of Neurology
   McGill University
*/

/* A partial volume estimation algorithm for single channel data. 
   Look for pve_aux.c for auxiliary functions. This file contains only main.
   Given all the necessary inputs , this program generates 4 - 7 minc volumes, 
   1 containing discrete classification (including mixed classes) and 3 containing 
   true ML PV estimates for pure classes and/or 3 containing approximative PV 
   estimates for pure classes. At some point there should be some documentation/user's
   guide also.
   Use pve3 for multispectral data.   
*/

#include "pve_aux.h"
#include <time_stamp.h>
#include "ParseArgv.h"

/* Moved here to avoid multiple definitions */
char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1];

int main(int argc, char** argv)

{
 
  char *history;
  char *pname = argv[0];
  char outfilename[255]; 
  char* filename = outfilename;

  int i,j,k,ii;                                    /* Loop variables */
  int error_code;  
  int sizes[MAX_DIMENSIONS];
  int changed_num;
  int changed_num_last;
  double cnt_num;
  char c, new_tissue_class,current_tissue_class;
  int iteration;
  double smlr;
  double subcort_val;

  BOOLEAN changed;
  BOOLEAN em = FALSE;
  BOOLEAN use_curve = FALSE;
  BOOLEAN use_subcort = FALSE;
  BOOLEAN sc_region = FALSE;

  static double mrf_params[4] = { 0.1, -2, -1, 1};
  static double curve_params[4] = { -2.0, -25, 1, 1.0 };
  static char* param_file = NULL;
  static char* seg_image = NULL;
  static char* tag_file = NULL;
  static char* mask_image = NULL;
  static char* mask_restrict = NULL;
  static char* curve_image = NULL;
  static char* subcort_image = NULL;
  static int est_params = FALSE;
  static int em_t = FALSE;
  static int mlestimates_t = MLONLY;
  static int classify = FALSE;
  static int clobber = FALSE;
  static int num_iterations = 0;

  static ArgvInfo argTable[] = {
    {"-file <paramfile>", ARGV_STRING, (char *) 1, (char *) &param_file,
     "File to use containing nuisance parameters"},
    {"-image <seg.mnc>", ARGV_STRING, (char *) 1, (char *) &seg_image,
     "Segmented image to use for nuisance parameter estimation"},
    {"-tags <tagfile.tag>", ARGV_STRING, (char *) 1, (char *) &tag_file,
     "File containg tag points to be used in nuisance parameter estimation"},
    {"-mask <mask.mnc>", ARGV_STRING,(char *) 1, (char *) &mask_image,
     "Brain mask used to specify region over which pve will be performed"},
    {"-restrict <restrict_mask.mnc>", ARGV_STRING,(char *) 1, 
     (char *) &mask_restrict,
     "Mask used to specify region over which nuisance paramater estimation will be performed"},
    {"-subcortical <subcort_mask.mnc>", ARGV_STRING, (char *) 1, (char *) &subcort_image,
     "Mask used to specify where subcortical estimation should be done"},
    {"-est_params",ARGV_CONSTANT, (char *) TRUE, (char *) &est_params,
     "Reestimate nuisance parameters every iteration"},
    {"-em",ARGV_CONSTANT, (char *) TRUE, (char *) &em_t,
     "Use expectation maximization in parameter estimation"},
    {"-noem",ARGV_CONSTANT, (char *) FALSE, (char *) &em_t,
     "Do not use expectation maximization in parameter estimation (default)"},
    {"-noml",ARGV_CONSTANT, (char *) NOML, (char *) &mlestimates_t,
     "Provide partial volume estimates based on simplified model"},
    {"-mlonly",ARGV_CONSTANT, (char *) MLONLY, (char *) &mlestimates_t,
     "Provide partial volume estimates based on complete model (default)"},
    {"-ml", ARGV_CONSTANT, (char *) ML, (char *) &mlestimates_t, 
     "Provide both simplified and complete model estimates"},
    {"-classify", ARGV_CONSTANT, (char *) TRUE, (char *) &classify, 
     "Provide a final classified image based on probabilistic maps (BG,CSF,GM,WM,SC)"},
    {"-mrf",ARGV_FLOAT, (char *) 4, (char *) mrf_params,
     "Parameters to use for Markovian Random Field: beta same similar different"},
    {"-curve",ARGV_STRING,(char *) 1, (char *) &curve_image,
     "Use specified curvature image in partial volume estimation"},    
    {"-weight",ARGV_FLOAT, (char *) 4, (char *) curve_params,
     "Parameters to use to weight curvature information"},
    {"-iterations",ARGV_INT,(char *) 0, (char *) &num_iterations,
     "Number of iterations to use unless convergence is reached"},
    {NULL, ARGV_HELP, (char *) NULL, (char *) NULL, "General options:"},
    {"-clobber",ARGV_CONSTANT,(char *) TRUE, (char *) &clobber,  
      "Overwrite any existing output files"},  
    {"-noclobber",ARGV_CONSTANT, (char *) FALSE, (char *) &clobber,
       "Do not overwrite any existing output files (default)"},
    {NULL, ARGV_END, NULL, NULL, NULL}
  };
  
  double mean[PURE_CLASSES + 1]; /*   nuisance parameters */
  double var[PURE_CLASSES + 1]; 
  double var_measurement;                    /* measurent noise */
  double pr_prior;

  double  slice_width[MAX_DIMENSIONS];
  double  width_stencil[MAX_DIMENSIONS*MAX_DIMENSIONS*MAX_DIMENSIONS];

  double val[CLASSES],mrf_probability[CLASSES],value; /*  Some temporary variables */
  int specified;

  /* Volumes to be read from files: in is the input image, and mask is the
     image used for masking non brain tissues out. */
  Volume volume_in, volume_mask, volume_restrict = NULL;

  Volume volume_likelihood[CLASSES]; /* Volumes for storing tissue type likelihoods. */

  Volume volume_classified;  /* Classifications */
  Volume volume_subcort = NULL;

  Volume volume_pve[PURE_CLASSES]; /* Volumes for partial volume estimates */
  Volume volume_pve_ml[PURE_CLASSES];  

  /* Intialize nuisance parameters */
  pr_prior = 0;

  /* For producing necessary info to the outputfiles. */
  history = time_stamp(argc, argv);

  if (ParseArgv(&argc, argv, argTable, 0))  {
    Usage_info(argv[0]);
    exit(EXIT_FAILURE);
  }

  if (em_t == TRUE) {
    em = TRUE;
  }

  //check whether output files exist
  if (!clobber) {
    filename = strcpy(filename,argv[2]); 
    if (file_exists(strcat(filename,"_disc.mnc"))) {
      fprintf(stderr,"Output file(s) exist!\n\n");
      return (2);  
    }
    if (mlestimates_t == ML || mlestimates_t == NOML) {      
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_csf.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_gm.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_wm.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      if( subcort_image ) {
        filename = strcpy(filename,argv[2]);  
        if (file_exists(strcat(filename,"_sc.mnc"))) {
          fprintf(stderr,"Output file(s) exist!\n\n");
          return (2);
        }
      }
    }
    if (mlestimates_t == MLONLY) {      
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_exactcsf.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_exactgm.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_exactwm.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
      if( subcort_image ) {
        filename = strcpy(filename,argv[2]);  
        if (file_exists(strcat(filename,"_exactsc.mnc"))) {
          fprintf(stderr,"Output file(s) exist!\n\n");
          return (2);
        }
      }
    }
    if (classify) {
      filename = strcpy(filename,argv[2]);  
      if (file_exists(strcat(filename,"_classify.mnc"))) {
        fprintf(stderr,"Output file(s) exist!\n\n");
        return (2);
      }
    }
  }

  if (argc != 3) {
    (void) fprintf(stderr,"Exactly one input file and one output prefix must be specified\n");
    Usage_info(pname);
    return(2);
  }

  if (mask_image == NULL) {
    (void) fprintf(stderr,"A mask volume must be specified\n");
    Usage_info(argv[0]);
    return(2);
  }

  specified = (param_file != NULL) + (seg_image != NULL) + (tag_file != NULL);
  if (specified != 1) {
    (void) fprintf(stderr,"Exactly one of -file, -image or -tags must be specified\n");
    Usage_info(pname);
    return(2);
  }

  if( input_volume( argv[1], 3, NULL, NC_UNSPECIFIED, FALSE, 0.0, 0.0,
                    TRUE, &volume_in, (minc_input_options *) NULL) != OK ) {
    fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
    fprintf(stderr, "Cannot open file %s.\n", argv[1] );
    return(2);
  } else if( get_volume_n_dimensions(volume_in) != 3 ) {
    return(2);
  }

  if( input_volume( mask_image, 3, NULL, NC_BYTE, FALSE, 0.0, 2.0, TRUE,
                    &volume_mask, (minc_input_options *) NULL) != OK ) {
    fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
    fprintf(stderr, "Cannot open file %s.\n", mask_image );
    return(2);
  } else if( get_volume_n_dimensions(volume_mask) != 3 ) {
    return(2);
  }

  if( mask_restrict ) {
    if( input_volume( mask_restrict, 3, NULL, NC_BYTE, FALSE, 0.0, 2.0, TRUE,
                      &volume_restrict, (minc_input_options *) NULL) != OK ) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
      fprintf(stderr, "Cannot open file %s.\n", mask_restrict );
      return(2);
    } else if( get_volume_n_dimensions(volume_restrict) != 3 ) {
      return(2);
    }
  }

  /* subcortical can be available in all cases: segmented image is given,
     with tags or with a parameter file. */

  if (subcort_image != NULL) {
    use_subcort = TRUE;
    if (input_volume(subcort_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0,
                     TRUE, &volume_subcort, (minc_input_options *) NULL) != OK)
      return (1); 
  }

  if (param_file != NULL) {
    error_code = Get_params_from_file(param_file,mean,var,&var_measurement);
    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
      fprintf(stderr,"Subfunction Get_params_from_file returned errorcode %d .\n",error_code);
      return(2);
    }
  } else if (seg_image != NULL) {

    Volume volume_seg;
    if(input_volume(seg_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0,
                    TRUE, &volume_seg, (minc_input_options *) NULL) != OK) {
       return(1);
    }

    error_code = Estimate_params_from_image(volume_in,volume_restrict,volume_subcort,
                                            volume_seg, mean, var, &var_measurement);   

    delete_volume( volume_seg );

    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_IMAGE);
      fprintf(stderr,"Subfunction Estimate_params_from_image returned errorcode %d .\n",error_code);
      return(2);   
    }
  } else if (tag_file != NULL) {
    error_code = Estimate_params_from_tags(tag_file,volume_in,mean,var);
    var_measurement = 0;
    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_TAG);
      fprintf(stderr,"Subfunction Estimate_params_from_tag returned errorcode %d .\n",error_code);
      return(2);  
    }
  }

  /* Set boundary points for extreme high voxel values. */
  char max_class = 0; 
  for(c = 1; c < PURE_CLASSES; c++) {
    if( mean[c] > mean[max_class] ) max_class = c;
  }

  /* Curvature volume used for biasing CSF estimation */
  Volume volume_curve;
  if (curve_image != NULL) {
    use_curve = TRUE;
    if(input_volume(curve_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                    TRUE, &volume_curve, (minc_input_options *) NULL) != OK) {
      return(1);
    }
  }

  /* Initialize required volumes and set their ranges for 
     getting rid of unnecessary surprises. */
  for( c = 0;c < CLASSES;c++) {
    // Use NC_BYTE instead of NC_UNSPECIFIED to save memory (a lot!) and
    // get nearly identical convergence.
    // volume_likelihood[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
    //volume_likelihood[c] = copy_volume_definition(volume_in, NC_BYTE, 
    volume_likelihood[c] = copy_volume_definition(volume_in, NC_SHORT, 
                                                  FALSE, 0.0 , 0.0);

    if( !volume_likelihood[c] || !volume_is_alloced( volume_likelihood[c] ) ) {
      fprintf(stderr,"Error: not enough memory to store volume likelihoods[%d].\n\n", c );
      exit(EXIT_FAILURE);
    }
    set_volume_real_range( volume_likelihood[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX );
  }

  /* Calculate the likelihoods */

  get_volume_sizes( volume_in, sizes );

  if (num_iterations <= 0) {
    num_iterations = MAX_ITERATIONS;
  }

  int pve_symmetric = 1;

  get_volume_separations( volume_in, slice_width );  /* Get slice separations in each direction */ 
  value = MIN3(slice_width[0],slice_width[1],slice_width[2]);/* And normalize them so that the   */ 
  slice_width[0] = slice_width[0]/value;              /* minimum is set to one.            */
  slice_width[1] = slice_width[1]/value;              /* All this for simplification of    */ 
  slice_width[2] = slice_width[2]/value;              /* usage of the MRFs.                */
  ii = 0;
  for(i = -1; i < 2; i++) {
    for(j = -1; j < 2; j++) {
      for(k = -1; k < 2; k++) {
        if( i == 0 && j == 0 && k == 0 ) {
          width_stencil[ii] = 0.0;
        } else {
          width_stencil[ii] = 1.0 / sqrt( pow(slice_width[0] * abs(i),2) +
                                          pow(slice_width[1] * abs(j),2) +
                                          pow(slice_width[2] * abs(k),2) );
        }
        ii++;
      }
    }
  }

  volume_classified = copy_volume_definition(volume_in, NC_BYTE, 
                                             FALSE, 0 , CLASSES); 
  set_volume_real_range( volume_classified, 0, CLASSES);

  printf("Beta %lf \n",mrf_params[BETA]);
  printf("Same %lf \n",mrf_params[SAME]);
  printf("Similar %lf \n",mrf_params[SMLR]);
  printf("Different %lf \n",mrf_params[DIFF]);
  printf("Exact ml estimates: %d \n",(int)(mlestimates_t==MLONLY||mlestimates_t==ML));
  printf("Simplified ml estimates: %d \n",(int)(mlestimates_t==NOML||mlestimates_t==ML));
  printf("Parameter updates: %d \n", em);

  initialize_Potts_table();

  iteration = 1;   /* Start the iterative algorithm */
  changed = TRUE;
  changed_num = 0;
  changed_num_last = -1;

  while(changed && (iteration <= num_iterations)) {
    printf("Iteration %d \n",iteration);

    /* Update the likelihoods with current parameters only if that is necessary */
    if(em || est_params || iteration == 1) {

      for( i = 0; i < sizes[0]; ++i) {
        for( j = 0; j < sizes[1]; ++j) {
          for( k = 0; k < sizes[2]; ++k ) {
            short mask_val = floor( get_volume_real_value(volume_mask,i,j,k,0,0) + 0.50 );
            if( mask_val > MASK_TR ) {
              if( mask_val < MASK_CHANGED_TR ) {
                if (use_subcort) {
                  if (get_volume_real_value(volume_subcort,i,j,k,0,0) > 0.5)
                    sc_region = TRUE;  // == 1 or 2 (could contain some CSF)
                  else
                    sc_region = FALSE;
                }
                value  = get_volume_real_value(volume_in,i,j,k,0,0);

                for(c = 1;c < PURE_CLASSES + 1; c++) {
                  if ((c == SCLABEL)&&(!sc_region)) {
                    val[c - 1] = 0;
                  } else {
                    val[c - 1] = Compute_Gaussian_likelihood(value,mean[c],
                                                             var[c] + var_measurement );
                  }
                }
                if (!sc_region) {
                  val[WMGMLABEL - 1] = Compute_marginalized_likelihood(value,mean[WMLABEL], mean[GMLABEL],
                                                                       var[WMLABEL], var[GMLABEL], 
                                                                       var_measurement);
                  val[WMSCLABEL - 1] = 0;
                  val[SCGMLABEL - 1] = 0;
                } else {
                  val[WMSCLABEL - 1] = Compute_marginalized_likelihood(value,mean[WMLABEL], mean[SCLABEL],
                                                                       var[WMLABEL], var[SCLABEL], 
                                                                       var_measurement);
                  val[SCGMLABEL - 1] = Compute_marginalized_likelihood(value,mean[SCLABEL], mean[GMLABEL],
                                                                       var[SCLABEL], var[GMLABEL], 
                                                                       var_measurement);
                  val[WMGMLABEL-1] = 0;
                }
                // maybe should have CSFSC class inside sc_region.
                val[GMCSFLABEL - 1] = Compute_marginalized_likelihood(value,mean[GMLABEL], mean[CSFLABEL],
                                                                      var[GMLABEL], var[CSFLABEL], 
                                                                      var_measurement );
                val[CSFBGLABEL - 1] = 0.0;  // don't allow BG inside mask
                if( value < mean[CSFLABEL] ) {
                  val[CSFLABEL-1] = 1.0;
                }

                if( Normalize(val,CLASSES) ) {
                  // All values are VERY_SMALL so pick something.
                  // ignore min since it will be BG.
                  if( value > mean[max_class] ) {
                    val[max_class-1] = 1.0;
                    Normalize(val,CLASSES);
                  } else {
                    val[CSFLABEL-1] = 1.0;
                    // printf( "Warning: Voxel (%d,%d,%d) out of range with value %g\n",
                    //         i, j, k, value );
                  }
                }

                for(c = 0;c < CLASSES;c++) {
                  set_volume_real_value(volume_likelihood[c], i , j , k, 0, 0, val[c]);
                }
                /* Note: no need to change volume_classified here for est_params.
                         It doesn't make any difference after convergence. */
                if( iteration == 1 ) {
                  c = Maxarg(val,CLASSES);
                  set_volume_real_value(volume_classified, i , j , k, 0, 0, c);
                }
              }
            } else {    /* if voxel is in the background, computations are not needed */
              set_volume_real_value(volume_classified, i , j ,k , 0 , 0, 0);
              for(c = 0;c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c], i , j , k, 0, 0, 0.0);
              }
            }
          }
        }
      }
    }

    /* ICM step */
    printf("ICM step \n");
    changed = FALSE;
    changed_num = 0;

    // Don't allow boundary voxels to change (should be outside the brain mask
    // anyways). This is to avoid checking for boundaries in Compute_mrf_probability
    // with the 3x3x3 stencil.
    for(i = 1; i < sizes[0]-1; ++i) {
      int iii = (pve_symmetric && iteration%2==0) ? sizes[0]-i-1 : i;
      for(j = 1; j < sizes[1]-1; ++j) {
        int jjj = (pve_symmetric && iteration%2==0) ? sizes[1]-j-1 : j;
        for(k = 1; k < sizes[2]-1; ++k) {
          int kkk = (pve_symmetric && iteration%2==0) ? sizes[2]-k-1 : k;
          short mask_val = floor( get_volume_real_value(volume_mask,iii,jjj,kkk,0,0) + 0.50 );
          if( mask_val > MASK_TR && mask_val < MASK_CHANGED_TR ) {

            if (use_subcort) {
              subcort_val = get_volume_real_value(volume_subcort,iii,jjj,kkk,0,0);
              if (subcort_val > 0.5 )
                sc_region = TRUE;
              else
                sc_region = FALSE;
            }

            // It's actually faster to recompute these values on the fly
            // at every ICM iteration as opposed to computing these values
            // only and saving them in a vector. The extra memory usage
            // makes the program run much slower at 0.5mm voxel resolution.
            double mrf_similar = mrf_params[SMLR];
            if( use_curve ) {
              double curve_val = get_volume_real_value(volume_curve,iii,jjj,kkk,0,0);
              if (curve_val < 0) {
                mrf_similar = curve_params[0]/
                              (1+exp(curve_params[1]*(fabs(curve_val)-curve_params[2])))
                              - curve_params[3];
              }
            }

            Compute_mrf_probability(mrf_probability, volume_classified,iii,jjj,kkk,
                                    width_stencil, mrf_params[BETA], mrf_params[SAME],
                                    mrf_similar, mrf_params[DIFF], pr_prior, sc_region );

            for(c = 0; c < CLASSES;c++) {
              if ((sc_region)&&(c+1 == WMGMLABEL))
                continue;
              else if ((!sc_region)&&((c+1 == WMSCLABEL)||(c+1 == SCGMLABEL)||(c+1 == SCLABEL)))
                continue;
              mrf_probability[c] *= get_volume_real_value(volume_likelihood[c],iii,jjj,kkk,0,0);
            }
            new_tissue_class = Maxarg(mrf_probability,CLASSES);

            current_tissue_class = get_volume_real_value(volume_classified,iii,jjj,kkk,0,0);
            if(new_tissue_class != current_tissue_class) {

              // Flag all neighbouring voxels as candidates to change on the
              // next iteration.
              int ii, jj, kk;
              for(ii = -1; ii <= 1; ii++) {
                for(jj = -1; jj <= 1; jj++) {
                  for(kk = -1; kk <= 1; kk++) {
                    if( get_volume_real_value(volume_mask,iii+ii,jjj+jj,kkk+kk,0,0) > MASK_TR ) {
                      set_volume_real_value(volume_mask,iii+ii,jjj+jj,kkk+kk,0,0,MASK_CHANGED);
                    }
                  }
                }
              }

              set_volume_real_value(volume_classified, iii, jjj, kkk, 0 ,0, new_tissue_class);
              changed = TRUE;
              changed_num++;
            } else {
              set_volume_real_value(volume_mask,iii,jjj,kkk,0,0,MASK_NOTCHANGED);
            }
            if(em) {   /* Store necessary probabilities for the parameter estimation step */
              set_volume_real_value(volume_mask,iii,jjj,kkk,0,0,MASK_CHANGED);
              Normalize(mrf_probability,CLASSES);
              for(c = 0; c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c],iii,jjj,kkk,0,0,mrf_probability[c]);
              }
            }
          }
        }
      }
    }

    /* And finally re-estimate the nuisance parameters if necessary for each pure-tissue class */
    if(em) {
      Parameter_estimation(volume_in,volume_restrict,volume_likelihood,mean,var,&var_measurement);
    } else if(est_params) {
      est_params = Parameter_estimation_classified(volume_in,volume_restrict,volume_subcort,
                                                   volume_classified,mean,var,&var_measurement);
    }

    if ((changed_num_last > -1)&&(changed_num >= changed_num_last)&&(!est_params)) {
      changed = FALSE;
    } else {
      changed_num_last = changed_num;
    }

    printf("changed_%d: %d\n",iteration,changed_num );
    fflush(stdout);
    iteration++;
  }
  
  /* Free memory that is no longer needed beyond this point. */
  
  for(c = 0;c < CLASSES;c++) {
    delete_volume(volume_likelihood[c]);
  }
  if (use_subcort) {
    delete_volume(volume_subcort);
  }
  if (use_curve) {
    delete_volume(volume_curve);
  }
  delete_volume(volume_mask);

  if( volume_restrict ) {
    delete_volume(volume_restrict);
  }

  /* write necessary files */

  if( classify ) {
    filename = strcpy(filename,argv[2]); 

    int max_range = (use_subcort) ? PURE_CLASSES : PURE_CLASSES-1;

    Volume final_cls = copy_volume_definition(volume_in, NC_BYTE, TRUE, 0, max_range); 
    set_volume_real_range( final_cls, 0, max_range);
    /* Final classification is based on MLONLY. */
    Compute_final_classification(volume_in,volume_classified,final_cls,mean,
                                 var, var_measurement );

    output_modified_volume(strcat(filename,"_classify.mnc"),
                           NC_BYTE, FALSE, 0, max_range, final_cls, argv[1],
                           history, (minc_output_options *) NULL);
    delete_volume(final_cls);
  }

  if(mlestimates_t==NOML || mlestimates_t==ML) {
    /* Allocate the memory for partial volume fractions */
    for(c = 0;c < PURE_CLASSES;c++) {
      volume_pve[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
                                             FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve[c],
                             LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
    printf("Computing fast partial volume fractions \n");
    Compute_partial_volume_vectors(volume_in,volume_classified,volume_pve,mean);

    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_wm.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve[WMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_gm.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve[GMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_csf.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve[CSFLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[2]); 
      output_modified_volume(strcat(filename,"_sc.mnc"),
                             NC_SHORT, FALSE, 0,0,volume_pve[SCLABEL - 1],argv[1],history,
                             (minc_output_options *) NULL);
    }

    for(c = 0;c < PURE_CLASSES;c++) { 
      delete_volume(volume_pve[c]);
    }
  }

  if(mlestimates_t==MLONLY || mlestimates_t==ML) {
    /* Allocate the memory for partial volume fractions */
    for(c = 0;c < PURE_CLASSES;c++) {
      volume_pve_ml[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
                                                FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve_ml[c],
                             LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
    printf("Computing ml partial volume fractions \n");
    Compute_partial_volume_vectors_ml(volume_in,volume_classified,volume_pve_ml,mean,
                                      var,var_measurement);

    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_exactwm.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve_ml[WMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_exactgm.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve_ml[GMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_exactcsf.mnc"),
                           NC_SHORT, FALSE, 0,0,volume_pve_ml[CSFLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[2]); 
      output_modified_volume(strcat(filename,"_exactsc.mnc"),
                             NC_SHORT, FALSE, 0,0,volume_pve_ml[SCLABEL - 1],argv[1],history,
                             (minc_output_options *) NULL);
    }
    for(c = 0;c < PURE_CLASSES;c++) { 
      delete_volume(volume_pve_ml[c]);
    }
  }

  filename = strcpy(filename,argv[2]); 
  output_modified_volume(strcat(filename,"_disc.mnc"),
                         NC_BYTE, FALSE, 0,CLASSES,volume_classified,argv[1],history,
                         (minc_output_options *) NULL);

  /* Free memory */
  delete_volume(volume_in);
  delete_volume(volume_classified);

  return(0);
}

