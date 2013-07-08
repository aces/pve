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
const char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1] = {{0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                                                           {0, 0, 0, 0, 0, 1, 0, 1, 0, 0},
                                                           {0, 0, 0, 0, 0, 1, 1, 0, 0, 1},
                                                           {0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
                                                           {0, 0, 0, 0, 0, 0, 0, 0, 1, 1},
                                                           {0, 1, 1, 0, 0, 0, 0, 0, 0, 0},
                                                           {0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
                                                           {1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                                                           {0, 0, 0, 1, 1, 0, 0, 0, 0, 0},
                                                           {0, 0, 1, 0, 1, 0, 0, 0, 0, 0}};

int main(int argc, char** argv)

{
 
  char *history;
  char *pname = argv[0];
  char outfilename[255]; 
  char* filename = outfilename;
  char tag_filename[255];
  char *ptag_filename = tag_filename;


  int i,j,k,ii;                                    /* Loop variables */
  int error_code;  
  int sizes[MAX_DIMENSIONS];
  int changed_num;
  int changed_num_last;
  double cnt_num;
  char c, new_tissue_class,current_tissue_class;
  int iteration;
  double curve_val;
  double smlr;
  double subcort_val;

  BOOLEAN changed;
  BOOLEAN em = FALSE;
  BOOLEAN use_curve = FALSE;
  BOOLEAN use_steady_state = TRUE;
  BOOLEAN use_subcort = FALSE;
  BOOLEAN sc_region = FALSE;

  static double mrf_params[4] = { 0.1, -2, -1, 1};
  static double curve_params[4] = { -2.0, -25, 1, 1.0 };
  static char* param_file = NULL;
  static char* seg_image = NULL;
  static char* tag_file = NULL;
  static char* mask_image = NULL;
  static char* curve_image = NULL;
  static char* subcort_image = NULL;
  static int est_params = FALSE;
  static int em_t = FALSE;
  static int mlestimates_t = MLONLY;
  static int classify = FALSE;
  static int clobber = FALSE;
  static int use_counter = FALSE;
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
    {"-count",ARGV_CONSTANT,(char *) TRUE, (char *) &use_counter,
     "Output volume containing change count"},
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
#if defined(NOT_IMPL)
  double pr_wm,pr_gm,pr_csf,pr_wmgm,pr_gmcsf,pr_csfbg; /* For later use */
#endif /* NOT_IMPL defined */

  double  slice_width[MAX_DIMENSIONS];
  double  width_stencil[MAX_DIMENSIONS*MAX_DIMENSIONS*MAX_DIMENSIONS];

  double val[CLASSES],mrf_probability[CLASSES],value; /*  Some temporary variables */
  int specified;

  Volume volume_in,volume_mask; /* Volumes to be read from files: 
                                   in is the input image, and mask is the
                                   image used for masking non brain tissues out. */
  
  
  Volume volume_likelihood[CLASSES]; /* Volumes for storing tissue type likelihoods. */

  Volume volume_classified;  /* Classifications */
  Volume volume_curve; /* Curvature volume used for biasing CSF estimation */
  Volume volume_count; /*Volume containing a count of # of times voxel class has changed */
  Volume volume_subcort = NULL;

  Volume volume_pve[PURE_CLASSES]; /* Volumes for partial volume estimates */
  Volume volume_pve_ml[PURE_CLASSES];  

  /* Intialize nuisance parameters */
#if defined(NOT_IMPL)
  pr_wm = 0.16667;      /* All tissue types are equally likely */
  pr_gm = 0.16667;
  pr_csf = 0.16667;
  pr_wmgm = 0.16667;
  pr_gmcsf = 0.16667;
  pr_csfbg = 0.16667;
#endif /* NOT_IMPL defined */
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

  error_code = Open_images(argv[1],mask_image,&volume_in, &volume_mask);
  if (error_code != 0) {
    fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
    fprintf(stderr,"Subfunction Open_images returned errorcode %d .\n",error_code);
    return(2);
  }
  get_volume_sizes( volume_in, sizes ); 

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

    /* subcortical only available when segmented image is given. */
    if (subcort_image != NULL) {
      use_subcort = TRUE;
      if (input_volume(subcort_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0,
                       TRUE, &volume_subcort, (minc_input_options *) NULL) != OK)
        return (1); 
    }

    error_code = Estimate_params_from_image(volume_in,volume_mask,volume_subcort,
                                            volume_seg,mean,var,&var_measurement);   

    delete_volume( volume_seg );

    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_IMAGE);
      fprintf(stderr,"Subfunction Estimate_params_from_image returned errorcode %d .\n",error_code);
      return(2);   
    }
  } else if (tag_file != NULL) {
    if(strcmp(tag_file,"default")) {
      ptag_filename = strcpy(ptag_filename,tag_file);
    }
    else {
      ptag_filename = strcpy(ptag_filename,DEFAULT_TAG_FILENAME);
    }
    error_code = Estimate_params_from_tags(tag_filename,volume_in,mean,var);
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

  if (curve_image != NULL) {
    use_curve = TRUE;
    if(input_volume(curve_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                    TRUE, &volume_curve, (minc_input_options *) NULL) != OK)
      return(1);

  }

  /* Initialize required volumes and set their ranges for 
     getting rid of unnecessary surprises. */
  for( c = 0;c < CLASSES;c++) {
    volume_likelihood[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
                                                  FALSE, 0.0 , 0.0);

    if( !volume_likelihood[c] || !volume_is_alloced( volume_likelihood[c] ) ) {
      fprintf(stderr,"Error: not enough memory to store volume likelihoods[%d].\n\n", c );
      exit(EXIT_FAILURE);
    }
    set_volume_real_range( volume_likelihood[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX );
  }
 
  if (use_counter) {    
    volume_count = copy_volume_definition(volume_in,NC_BYTE,FALSE,0,100);
    set_volume_real_range(volume_count,0,100);

    for( i = 0; i < sizes[0]; ++i) {
      for( j = 0; j < sizes[1]; ++j) {
        for( k = 0; k < sizes[2]; ++k ) {
          set_volume_real_value(volume_count,i,j,k,0,0,0.0);
        }
      }
    }
  }

  if (num_iterations > 0)
    use_steady_state = FALSE;
  else
    num_iterations = MAX_ITERATIONS;

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
          TRUE, 0 , CLASSES); 
  set_volume_real_range( volume_classified, 0, CLASSES);

  printf("Same %lf \n",mrf_params[SAME]);
  printf("Similar %lf \n",mrf_params[SMLR]);
  printf("Different %lf \n",mrf_params[DIFF]);
  printf("Exact ml estimates: %d \n",(int)(mlestimates_t==MLONLY||mlestimates_t==ML));
  printf("Simplified ml estimates: %d \n",(int)(mlestimates_t==NOML||mlestimates_t==ML));
  printf("Parameter updates: %d \n", em);

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
            if( i == 0 || j == 0 || k == 0 || i == (sizes[0]-1) ||
                j == (sizes[1]-1) || k == (sizes[2]-1) ) {
              set_volume_real_value(volume_classified, i , j , k, 0, 0, 0 );
              for(c = 0;c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c], i , j , k, 0, 0, 0.0);
              }
            } else {
              if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) { 

	        if (use_subcort) {
		  if (get_volume_real_value(volume_subcort,i,j,k,0,0) > 0.5)
		    sc_region = TRUE;
		  else
		    sc_region = FALSE;
	        }
                value  = get_volume_real_value(volume_in,i,j,k,0,0);
                for(c = 1;c < PURE_CLASSES + 1; c++) {
		  if ((c == SCLABEL)&&(!sc_region))
		    val[c - 1] = 0;
		  else
		    val[c - 1] = Compute_Gaussian_likelihood(value,mean[c],
							     var[c] + var_measurement );
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

	        val[GMCSFLABEL - 1] = Compute_marginalized_likelihood(value,mean[GMLABEL], mean[CSFLABEL],
								      var[GMLABEL], var[CSFLABEL], 
								      var_measurement );
//              val[CSFBGLABEL - 1] = Compute_marginalized_likelihood(value, mean[CSFLABEL], mean[BGLABEL], 
//                                                                    var[CSFLABEL] ,var[BGLABEL], 
//                                                                    var_measurement );   
                val[CSFBGLABEL - 1] = 0.0;  // don't allow BG inside mask

                if( Normalize(val,CLASSES) ) {
                  // All values are VERY_SMALL so pick something.
                  // ignore min since it will be BG.
                  if( value > mean[max_class] ) {
                    val[max_class-1] = 1.0;
                    Normalize(val,CLASSES);
                  } else {
                    val[CSFLABEL-1] = 1.0;
                    printf( "Warning: Voxel (%d,%d,%d) out of range with value %g\n",
                            i, j, k, value );
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
    }

    /* ICM step */
    printf("ICM step \n");
    changed = FALSE;
    changed_num = 0;
    for(i = 0; i < sizes[0]; ++i) {
      for(j = 0; j < sizes[1]; ++j) {
        for(k = 0; k < sizes[2]; ++k) {
          current_tissue_class = get_volume_real_value(volume_classified,i,j,k,0,0);

// NOTE: A BG voxel can never change to a new type this way!!!! (use mask instead)
          if(current_tissue_class != 0) {

            if (use_subcort) {
              subcort_val = get_volume_real_value(volume_subcort,i,j,k,0,0);
              if (subcort_val > 0.5 )
                sc_region = TRUE;
              else
                sc_region = FALSE;
            }

            //use curvature info to weight the MRF spatially if provided
            double mrf_similar = mrf_params[SMLR];
            if (use_curve) {
              curve_val = get_volume_real_value(volume_curve,i,j,k,0,0);
              if (curve_val < 0) {
                mrf_similar = curve_params[0]/(1+exp(curve_params[1]*(fabs(curve_val)-curve_params[2]))) 
                              - curve_params[3];
              }
            }

            Compute_mrf_probability(mrf_probability, &volume_classified,i,j,k,
                                    width_stencil, mrf_params[BETA], mrf_params[SAME],
                                    mrf_similar, mrf_params[DIFF],pr_prior );

            if (sc_region) {
              mrf_probability[WMGMLABEL-1] = 0;
            } else {
              mrf_probability[WMSCLABEL-1] = 0;
              mrf_probability[SCGMLABEL-1] = 0;
              mrf_probability[SCLABEL-1] = 0;
            }
 
            Normalize(mrf_probability,CLASSES);

            for(c = 0; c < CLASSES;c++) {

              if ((sc_region)&&(c+1 == WMGMLABEL))
                continue;
              else if ((!sc_region)&&((c+1 == WMSCLABEL)||(c+1 == SCGMLABEL)||(c+1 == SCLABEL)))
                continue;
              mrf_probability[c] = get_volume_real_value(volume_likelihood[c],i,j,k,0,0) * 
                                 mrf_probability[c];
            }
            Normalize(mrf_probability,CLASSES);
            new_tissue_class = Maxarg(mrf_probability,CLASSES);

            if(new_tissue_class != current_tissue_class) {
              set_volume_real_value(volume_classified, i, j, k, 0 ,0, new_tissue_class);
              if (use_counter) {
                cnt_num = get_volume_real_value(volume_count,i,j,k,0,0);
                cnt_num += 1.0;
                set_volume_real_value(volume_count,i,j,k,0,0,cnt_num);
              }
              changed = TRUE;
              changed_num++;
            }
            if(em) {   /* Store necessary probabilities for the parameter estimation step */
              for(c = 0; c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c],i,j,k,0,0,mrf_probability[c]);
              }
            }
          }
        }
      }
    }

    /* And finally re-estimate the nuisance parameters if necessary for each pure-tissue class */
    if(em) {
      Parameter_estimation(volume_in,volume_mask,volume_likelihood,mean,var,&var_measurement);
    } else if(est_params) {
      est_params = Parameter_estimation_classified(volume_in,volume_mask,volume_subcort,
                                                   volume_classified,mean,var,&var_measurement);
    }

    if (use_steady_state) {
      if ((changed_num_last > -1)&&(changed_num >= changed_num_last)&&(!est_params))
        changed = FALSE;
      else 
        changed_num_last = changed_num;
    }

    printf("changed_%d: %d\n",iteration,changed_num);
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

  /* write necessary files */

  //output count of of the number of times every voxel has changed class
  if (use_counter) {
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_count.mnc"),
                           NC_UNSPECIFIED, FALSE, 0,0,volume_count,argv[1],history,
                           (minc_output_options *) NULL); 
    delete_volume(volume_count);
  }

  if( classify ) {
    filename = strcpy(filename,argv[2]); 

    Volume final_cls = copy_volume_definition(volume_in, NC_BYTE, 
                                              TRUE, 0, 255); 
                                              //TRUE, 0, PURE_CLASSES); 
    //set_volume_real_range( final_cls, 0, PURE_CLASSES);
    set_volume_real_range( final_cls, 0, 255);
    Compute_final_classification(volume_in,volume_classified,final_cls,mean,var);

    output_modified_volume(strcat(filename,"_classify.mnc"),
                           //NC_BYTE, FALSE, 0, PURE_CLASSES, final_cls, argv[1],
                           NC_BYTE, FALSE, 0, 255, final_cls, argv[1],
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
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[WMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_gm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[GMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_csf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[CSFLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[2]); 
      output_modified_volume(strcat(filename,"_sc.mnc"),
                             NC_UNSPECIFIED, FALSE, 0,0,volume_pve[SCLABEL - 1],argv[1],history,
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
                           NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[WMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_exactgm.mnc"),
                           NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[GMLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    filename = strcpy(filename,argv[2]); 
    output_modified_volume(strcat(filename,"_exactcsf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[CSFLABEL - 1],argv[1],history,
                           (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[2]); 
      output_modified_volume(strcat(filename,"_exactsc.mnc"),
                             NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[SCLABEL - 1],argv[1],history,
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

