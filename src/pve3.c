/* ----------------------------- MNI Header -----------------------------------
@NAME       :  pve3
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Partial volume classification of the multi-channel MRI.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Feb 2002 (Jussi Tohka, jupeto@cs.tut.fi) 
@MODIFIED   : Last modification in Aug 2004 by Vivek Singh
---------------------------------------------------------------------------- */

/* Partial volume estimation for multi channel data 
   Look for pve3_aux.c for auxiliary functions. 
   This file contains only main.
   Given all the necessary inputs , this program generates 4 - 7 minc volumes, 
   1 containing discrete classification (including mixed classes) and 3 containing 
   true ML PV estimates for pure classes and/or 3 containing approximative PV 
   estimates for pure classes. At some point there should be some documentation/user's
   guide also.
   Use pve for singlespectral data.
*/


#include "pve3_aux.h"
#include <time_stamp.h>
#include "ParseArgv.h"

int main(int argc, char** argv)

{
 
  char *history;
  char *pname = argv[0];
  char outfilename[255]; 
  char* filename = outfilename;
  char tag_filename[255];
  char *ptag_filename = tag_filename;

  int i,j,k;                                    /* Loop variables */
  int error_code;  
  int sizes[MAX_DIMENSIONS];
  int changed_num;
  int changed_num_last;
  int cnt_num;
  char c, new_tissue_class,current_tissue_class;
  int iteration;
  double curve_val;
  double smlr;
  double subcort_val;

  BOOLEAN changed;
  BOOLEAN em = FALSE;
  BOOLEAN mlestimates = FALSE;
  BOOLEAN mlestimates_only = FALSE;
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
  static int em_t = FALSE;
  static int mlestimates_t = NOML;
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
    {"-em",ARGV_CONSTANT, (char *) TRUE, (char *) &em_t,
     "Use expectation maximization in parameter estimation"},
    {"-noem",ARGV_CONSTANT, (char *) FALSE, (char *) &em_t,
     "Do not use expectation maximization in parameter estimation (default)"},
    {"-noml",ARGV_CONSTANT, (char *) NOML, (char *) &mlestimates_t,
     "Provide partial volume estimates based on simplified model"},
    {"-mlonly",ARGV_CONSTANT, (char *) MLONLY, (char *) &mlestimates_t,
     "Provide partial volume estimates based on complete model"},
    {"-ml", ARGV_CONSTANT, (char *) ML, (char *) &mlestimates_t, 
     "Provide both simplified and complete model estimates"},
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

  pVector mean[PURE_CLASSES + 1]; /*   nuisance parameters */
  pMatrix var[PURE_CLASSES + 1];  /* The variable mean contains four pointers to mean vectors */
                   /* and the variable var contains four pointers to var matrices. */  
                   /* This trick is made to simplify the use of labels.    */

  

  Vector3D mean_wm,mean_gm,mean_csf,mean_bg,mean_sc;
  Matrix3D var_wm,var_gm,var_csf,var_bg,var_sc;
  Matrix3D var_measurement;	/* measurent noise */
  double pr_prior;
#if defined(NOT_IMPL)
  double pr_wm,pr_gm,pr_csf,pr_wmgm,pr_gmcsf,pr_csfbg; /* For later use */
#endif /* NOT_IMPL defined */
  double beta = 0.1;

  double  slice_width[MAX_DIMENSIONS];

  double val[CLASSES],mrf_probability[CLASSES],value; /*  Some temporary variables */
  Matrix3D varsum;
  Vector3D intensity;
  int specified;

  Volume volume_inT1, volume_inT2, volume_inPD;
  Volume volume_mask;                       /* Volumes to be read from files: 
                                               ins are the input images, and mask is the
                                               image used for masking non brain tissues out. */


  Volume volume_likelihood[CLASSES]; /* Volumes for storing tissue type likelihoods. */
 

  Volume volume_classified;  /* Classifications */
  Volume volume_pve[PURE_CLASSES]; /* Volumes for partial volume estimates */
  Volume volume_pve_ml[PURE_CLASSES];  
  Volume volume_curve;  /* Curvature volume used for biasing CSF estimation */
  Volume volume_count; /*Volume containing a count of # of times voxel class has changed */
  Volume volume_subcort = NULL;

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

  if (mlestimates_t == MLONLY) {
    mlestimates_only = TRUE;
  }
  else if (mlestimates_t == ML) {
    mlestimates = TRUE;
  }
	
  if (em_t == TRUE) {
    em = TRUE;
  }
  

  //check wheter output files exist
  if (!clobber) {
    filename = strcpy(filename,argv[4]); 
    if (file_exists(strcat(filename,"_disc.mnc"))) {
      fprintf(stderr,"Output file(s) exist!\n\n");
      return (2);  
    }
    if (!mlestimates_only) {      
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_csf.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_gm.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_wm.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_sc.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
    }
    if (mlestimates) {
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_exactcsf.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_exactgm.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_exactwm.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
      filename = strcpy(filename,argv[4]);  
      if (file_exists(strcat(filename,"_exactsc.mnc"))) {
	fprintf(stderr,"Output file(s) exist!\n\n");
	return (2);
      }
    }
  }

  /* Initialize pointers to vectors and matrices holding nuisance parameters */ 
  var[WMLABEL] = &var_wm[0];
  var[GMLABEL] = &var_gm[0];
  var[CSFLABEL] = &var_csf[0];
  var[BGLABEL] = &var_bg[0];
  var[SCLABEL] = &var_sc[0];
  mean[WMLABEL] = &mean_wm[0];
  mean[GMLABEL] = &mean_gm[0];
  mean[CSFLABEL] = &mean_csf[0];
  mean[BGLABEL] = &mean_bg[0];
  mean[SCLABEL] = &mean_sc[0];

  if (argc != 5) {
    (void) fprintf(stderr,"Exactly three input files and one output prefix must be specified\n");
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

  error_code = Open_images3(argv[1],argv[2],argv[3],mask_image,
			    &volume_inT1,&volume_inT2,&volume_inPD, &volume_mask);
  if( error_code != 0) {
    printf("\n %s \n",ERROR_INPUT_FILE);
    printf("Subfunction Open_images3 returned errorcode %d .\n",error_code);
    return(3);
  }
  
  if (subcort_image != NULL) {
    use_subcort = TRUE;
    if (input_volume(subcort_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0,
		     TRUE, &volume_subcort, (minc_input_options *) NULL) != OK)
      return (1); 
  }
  
  if (param_file != NULL) {
    error_code = Get_params_from_file3(param_file,mean,var,var_measurement);
    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_FILE);
      fprintf(stderr,"Subfunction Get_params_from_file returned errorcode %d .\n",error_code);
      return(2);
    }
  }
  else if (seg_image != NULL) {

     error_code = Estimate_params_from_image3(volume_inT1,volume_inT2,volume_inPD,volume_mask,volume_subcort,
					    seg_image,mean,var, var_measurement);

     if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_IMAGE);
      fprintf(stderr,"Subfunction Estimate_params_from_image returned errorcode %d .\n",error_code);
        return(2);   
    }
  }
  else if (tag_file != NULL) {
    if(strcmp(tag_file,"default")) {
      ptag_filename = strcpy(ptag_filename,tag_file);
    }
    else {
      ptag_filename = strcpy(ptag_filename,DEFAULT_TAG_FILENAME);
    }

    error_code = Estimate_params_from_tags3(tag_filename,volume_inT1,
					    volume_inT2,volume_inPD,mean,var);

    SetZero(var_measurement);	
    if(error_code != 0) {
      fprintf(stderr,"\n %s \n",ERROR_PARAMS_TAGS);
      fprintf(stderr,"Subfunction Estimate_params_from_tag returned errorcode %d .\n",error_code);
      return(2);  
    }
  }

  if (curve_image != NULL) {
    use_curve = TRUE;
    if(input_volume(curve_image,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
		    TRUE, &volume_curve, (minc_input_options *) NULL) != OK)
      return(1);

  }


  /* Initialize required volumes and set their ranges for 
     getting rid of unncessary surprises. */
  for( c = 0;c < CLASSES;c++) {
    volume_likelihood[c] = copy_volume_definition(volume_inT1, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);  
    set_volume_real_range( volume_likelihood[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX );
  }
 
  volume_classified = copy_volume_definition(volume_inT1, NC_BYTE, 
	  TRUE, 0 , CLASSES); 
  set_volume_real_range( volume_classified, 0, CLASSES);

  /* Calculate the likelihoods */
  
  get_volume_sizes( volume_inT1, sizes ); 

  if (use_counter) {    
    volume_count = copy_volume_definition(volume_inT1,NC_BYTE,FALSE,0,100);
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

  get_volume_separations( volume_inT1 , slice_width );  /* Get slice separations in each direction */ 
  value = MIN3(slice_width[0],slice_width[1],slice_width[2]);/* And normalize them so that the   */ 
  slice_width[0] = slice_width[0]/value;                    /* minimum is set to one.            */
  slice_width[1] = slice_width[1]/value;                    /* All this for simplification of    */ 
  slice_width[2] = slice_width[2]/value;                    /* usage of the MRFs.                */
  printf("Same %lf \n",mrf_params[SAME]);
  printf("Similar %lf \n",mrf_params[SMLR]);
  printf("Different %lf \n",mrf_params[DIFF]);
  printf("Only ml estimates: %d \n",mlestimates_only);
  printf("ml estimates: %d \n",mlestimates);
  printf("Parameter updates: %d \n", em);
  iteration = 1;  /*Start the iterative algorithm */  
  changed = TRUE;
  changed_num  = 0;
  changed_num_last = -1;
  
  while(changed && (iteration < num_iterations)) {
    printf("Iteration %d \n",iteration);

    /* Update the likelihoods with current parameters only if that is necessary */
    if(em || iteration == 1) {
      for( i = 0; i < sizes[0]; ++i) {
        for( j = 0; j < sizes[1]; ++j) {
          for( k = 0; k < sizes[2]; ++k ) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) { 

	      if (use_subcort) {
		if (get_volume_real_value(volume_subcort,i,j,k,0,0) > 0.5)
		  sc_region = TRUE;
		else
		  sc_region = FALSE;
	      }
              intensity[0] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              intensity[1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              intensity[2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              for(c = 1;c < PURE_CLASSES + 1;c++) { 

		if ((c == SCLABEL)&&(!sc_region))
		  val[c - 1] = 0;
		else {
		  AddMatrices(var[c],var_measurement,varsum);     
		  val[c - 1] = Compute_Gaussian_likelihood3(intensity,mean[c],varsum );
		}
	      }

	      if (!sc_region) {
		
		val[WMGMLABEL - 1] = Compute_marginalized_likelihood3(intensity,mean[WMLABEL], mean[GMLABEL],
								      var[WMLABEL], var[GMLABEL], 
								      var_measurement, INTERVALS );
		val[WMSCLABEL - 1] = 0;
		val[SCGMLABEL - 1] = 0;

	      }
	      else {
		val[WMSCLABEL - 1] = Compute_marginalized_likelihood3(intensity,mean[WMLABEL], mean[SCLABEL],
								     var[WMLABEL], var[SCLABEL], 
								     var_measurement, INTERVALS );
		val[SCGMLABEL - 1] = Compute_marginalized_likelihood3(intensity,mean[SCLABEL], mean[GMLABEL],
								     var[SCLABEL], var[GMLABEL], 
								     var_measurement, INTERVALS );
		val[WMGMLABEL-1] = 0;
		

	      }
              val[GMCSFLABEL - 1] = Compute_marginalized_likelihood3(intensity,mean[GMLABEL], mean[CSFLABEL],
                                                                var[GMLABEL], var[CSFLABEL], 
                                                                var_measurement, INTERVALS );
              val[CSFBGLABEL - 1] = Compute_marginalized_likelihood3(intensity, mean[CSFLABEL], mean[BGLABEL], 
                                                            var[CSFLABEL] ,var[BGLABEL], 
                                                            var_measurement, INTERVALS );         
              Normalize(val,CLASSES);
             
              for(c = 0;c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c], i , j , k, 0, 0, val[c]);
              }
              if(iteration == 1) {
                c = Maxarg(val,CLASSES);
                set_volume_real_value(volume_classified, i , j , k, 0, 0, c);
              }
             
            }
            else {    /* if voxel is in the background, computations are not needed */
              set_volume_real_value(volume_classified, i , j ,k , 0 , 0, 0);
              for(c = 0;c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c], i , j , k, 0, 0, 0.0);
              }
            }
          }
        }
      }

    }
    /* ICM-step */     
    printf("ICM step \n");
    changed = FALSE;
    changed_num = 0;
    for(i = 0; i < sizes[0]; ++i) {
      for(j = 0; j < sizes[1]; ++j) {
        for(k = 0; k < sizes[2]; ++k) {
	  mrf_params[SMLR] = -1;
          current_tissue_class = get_volume_real_value(volume_classified,i,j,k,0,0);
	  if(current_tissue_class != 0) {

            if (use_subcort) {
              subcort_val = get_volume_real_value(volume_subcort,i,j,k,0,0);
              if (subcort_val > 0.5 )
                sc_region = TRUE;
              else
                sc_region = FALSE;
            }

	    //use curvature info to weight the MRF spatially if provided
	    if (use_curve) {
	      curve_val = get_volume_real_value(volume_curve,i,j,k,0,0);
	      if (curve_val < 0) {
		mrf_params[SMLR] = curve_params[0]/(1+exp(curve_params[1]*(fabs(curve_val)-curve_params[2]))) 
		  - curve_params[3];
	      }
	    }
            for(c = 0;c < CLASSES ;c++) {
	      
	      if (sc_region) {
		if (c+1 == WMGMLABEL) {
		  mrf_probability[c] = 0;
		  continue;
		}
	      }
	     
	      if (!sc_region) {
		if ((c+1 == WMSCLABEL)||(c+1 == SCGMLABEL)||(c+1 == SCLABEL)) {
		  mrf_probability[c] = 0;
		  continue;
		}
	      }

	      if (((c+1) == GMCSFLABEL)||((c+1) == CSFLABEL))
		smlr = mrf_params[SMLR];
	      else
		smlr = -1;
	      
	      mrf_probability[c] = Compute_mrf_probability(c + 1,&volume_classified,i,j,k,
							   slice_width, mrf_params[BETA], mrf_params[SAME], 
							   smlr, mrf_params[DIFF],pr_prior,sizes);  
	      
	      
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
            if(em) {    /* Store necessary probabilities for the parameter estimation step */
              for(c = 0; c < CLASSES;c++) {
                set_volume_real_value(volume_likelihood[c],i,j,k,0,0,mrf_probability[c]);
              }
            }
          }
        }
      }
    }
    if(em) {
      Parameter_estimation3(volume_inT1,volume_inT2,volume_inPD,volume_mask,
                           volume_likelihood,mean,var,var_measurement);
    } 

    if (use_steady_state) {
      if ((changed_num_last > -1)&&(changed_num >= changed_num_last))
	changed = FALSE;
      else 
	changed_num_last = changed_num;  
    }
    printf("changed_%d: %d\n",iteration,changed_num);
    iteration++;
   
  }
  
  /* Finally estimate the partial volume vectors 
     but first some memory is freed.*/
  
  for(c = 0;c < CLASSES;c++) {
    delete_volume(volume_likelihood[c]);
  }
  
  /* Allocate the memory for partial volume fractions */
  printf("Computing partial volume fractions \n");
  if(!mlestimates_only) {
    for(c = 0;c < PURE_CLASSES;c++) {
      volume_pve[c] = copy_volume_definition(volume_inT1, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
   
  
    Compute_partial_volume_vectors3(volume_inT1,volume_inT2,volume_inPD,
                                 volume_classified,volume_pve,
                                 mean,var,var_measurement);
    
  }
  if(mlestimates) {
    for(c = 0;c < PURE_CLASSES;c++) {
      volume_pve_ml[c] = copy_volume_definition(volume_inT1, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve_ml[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
  
    Compute_partial_volume_vectors3_ml(volume_inT1,volume_inT2,volume_inPD,
                                 volume_classified,volume_pve_ml,
                                 mean,var,var_measurement);
  }
 
 

  /* write necessary files */

  filename = strcpy(filename,argv[4]); 
  output_modified_volume(strcat(filename,"_disc.mnc"),
                          NC_BYTE, FALSE, 0,CLASSES,volume_classified,argv[1],history,
			  (minc_output_options *) NULL); 

  
  if(!mlestimates_only) {
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_wm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[WMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_gm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[GMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_csf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[CSFLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[4]); 
      output_modified_volume(strcat(filename,"_sc.mnc"),
			     NC_UNSPECIFIED, FALSE, 0,0,volume_pve[SCLABEL - 1],argv[1],history,
			     (minc_output_options *) NULL);
    }
  }
  if(mlestimates) {
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactwm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[WMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactgm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[GMLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactcsf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[CSFLABEL - 1],argv[1],history,
                         (minc_output_options *) NULL);
    if (use_subcort) {
      filename = strcpy(filename,argv[4]); 
      output_modified_volume(strcat(filename,"_exactsc.mnc"),
			     NC_UNSPECIFIED, FALSE, 0,0,volume_pve[SCLABEL - 1],argv[1],history,
			     (minc_output_options *) NULL);
    }
  } 

  //output count of of the number of times every voxel has changed class
  if (use_counter) {
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_count.mnc"),
			   NC_UNSPECIFIED, FALSE, 0,0,volume_count,argv[1],history,
			   (minc_output_options *) NULL); 
  }
  
  /* Free memory */
  
  delete_volume(volume_inT1);
  delete_volume(volume_inT2);
  delete_volume(volume_inPD);
  delete_volume(volume_mask);
  delete_volume(volume_classified);
  if(!mlestimates_only) {
    for(c = 0;c < PURE_CLASSES;c++) { 
      delete_volume(volume_pve[c]);
    }
  }
  if(mlestimates) {
    for(c = 0;c < PURE_CLASSES;c++) { 
      delete_volume(volume_pve_ml[c]);
    }
  }  
  return(0);
}

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
