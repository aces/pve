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
@MODIFIED   : The last modifacation in May 2002 by Jussi Tohka
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


int main(int argc, char** argv)

{
 
  char *history;
  char outfilename[255]; 
  char* filename = outfilename;
  char tag_filename[255];
  char *ptag_filename = tag_filename;

  int i,j,k;                                    /* Loop variables */
  int error_code;  
  int sizes[MAX_DIMENSIONS];
  char c, new_tissue_class,current_tissue_class;
  int iteration;

  BOOLEAN changed;
  BOOLEAN em = FALSE;
  BOOLEAN mlestimates = FALSE;
  BOOLEAN mlestimates_only = FALSE;

  pVector mean[PURE_CLASSES + 1]; /*   nuisance parameters */
  pMatrix var[PURE_CLASSES + 1];  /* The variable mean contains four pointers to mean vectors */
                   /* and the variable var contains four pointers to var matrices. */  
                   /* This trick is made to simplify the use of labels.    */

  

  Vector3D mean_wm,mean_gm,mean_csf,mean_bg;
  Matrix3D var_wm,var_gm,var_csf,var_bg;
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

  int same = -2;                          /* For defining the Potts model */
  int similar = -1;
  int different = 1;

  Volume volume_inT1, volume_inT2, volume_inPD;
  Volume volume_mask;                       /* Volumes to be read from files: 
                                               ins are the input images, and mask is the
                                               image used for masking non brain tissues out. */


  Volume volume_likelihood[CLASSES]; /* Volumes for storing tissue type likelihoods. */
 

  Volume volume_classified;  /* Classifications */
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

  /* Initialize pointers to vectors and matrices holding nuisance parameters */ 
  var[WMLABEL] = &var_wm[0];
  var[GMLABEL] = &var_gm[0];
  var[CSFLABEL] = &var_csf[0];
  var[BGLABEL] = &var_bg[0];
  mean[WMLABEL] = &mean_wm[0];
  mean[GMLABEL] = &mean_gm[0];
  mean[CSFLABEL] = &mean_csf[0];
  mean[BGLABEL] = &mean_bg[0];

  if(argc < 2) {
    printf("%s", ERROR_TOO_FEW);
    return(1);
  }

  if(!(strcmp(argv[1],"-help"))) {
    Display_help();
    return(0);
  }

  else if(!(strcmp(argv[1], "-file"))) {
    if(argc < 8) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Get_params_from_file3(argv[7],mean,var,var_measurement); 
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_FILE);
        printf("Subfunction Get_params_from_file3 returned errorcode %d .\n",error_code);
        return(2);
      }
     /* Read the image data for input volume and also for brainmask volume. */

      error_code = Open_images3(argv[2],argv[3],argv[4],argv[5],
                           &volume_inT1,&volume_inT2,&volume_inPD, &volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images3 returned errorcode %d .\n",error_code);
        return(3);
      }
    }
  }
  else if(!(strcmp(argv[1], "-image"))) {
    if(argc < 8) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Open_images3(argv[2],argv[3],argv[4],argv[5],
                           &volume_inT1,&volume_inT2,&volume_inPD, &volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images3 returned errorcode %d .\n",error_code);
        return(3);
      }

      error_code = Estimate_params_from_image3(volume_inT1,volume_inT2,volume_inPD,volume_mask,argv[7],mean,
                                               var, var_measurement);
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_IMAGE);
        printf("Subfunction Estimate_params_from_image3 returned errorcode %d .\n",error_code);
        return(2);  
      }
    }
  }
   else if(!(strcmp(argv[1], "-tags"))) {
    if(argc < 8) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Open_images3(argv[2],argv[3],argv[4],argv[5],
                                &volume_inT1,&volume_inT2,&volume_inPD,&volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images3 returned errorcode %d .\n",error_code);
        return(3);
      }
      if(strcmp(argv[7],"default")) {
        ptag_filename = strcpy(ptag_filename,argv[7]);
      }
      else {
        ptag_filename = strcpy(ptag_filename,DEFAULT_TAG_FILENAME);
      }
      error_code = Estimate_params_from_tags3(tag_filename,volume_inT1,
                                              volume_inT2,volume_inPD,mean,var);
      SetZero(var_measurement);
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_TAGS);
        printf("Subfunction Estimate_params_from_tags3 returned errorcode %d .\n",error_code);
        return(2);  
      }
    }
  }
  else {
    printf("%s",ERROR_SECOND_ARG);
    return(4);
  }


  if(argc > 8) {
    if(!(strcmp(argv[8], "-em"))) {
      em = TRUE;                  /* Turn on the parameter updates */
    }
  }
  if(argc > 9) {
    if(!(strcmp(argv[9], "-ml"))) {
      mlestimates = TRUE;
    }
    if(!(strcmp(argv[9], "-mlonly"))) {
      mlestimates_only = TRUE;
      mlestimates = TRUE;
    }
  } 

  if(argc > 14) {
    if(!(strcmp(argv[10], "-mrf"))) {
      sscanf(argv[11],"%lf",&beta);
      sscanf(argv[12],"%d",&same);
      sscanf(argv[13],"%d",&similar);
      sscanf(argv[14],"%d",&different);
    }
  } 

  
 

  /* For producing necessary info to the outputfiles. */
  history = time_stamp(argc, argv);
  
  /* Initialize required volumes and set their ranges for 
     getting rid of unncessary surprises. */

  for( c = 0;c < CLASSES;c++) {
    volume_likelihood[c] = copy_volume_definition(volume_inT1, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);  
    set_volume_real_range( volume_likelihood[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX );
  }
 
  volume_classified = copy_volume_definition(volume_inT1, NC_BYTE, 
	  TRUE, 0 , 6); 
  set_volume_real_range( volume_classified, 0,6 );

  /* Calculate the likelihoods */
  
  get_volume_sizes( volume_inT1, sizes ); 

  get_volume_separations( volume_inT1 , slice_width );  /* Get slice separations in each direction */ 
  value = MIN3(slice_width[0],slice_width[1],slice_width[2]);/* And normalize them so that the   */ 
  slice_width[0] = slice_width[0]/value;                    /* minimum is set to one.            */
  slice_width[1] = slice_width[1]/value;                    /* All this for simplification of    */ 
  slice_width[2] = slice_width[2]/value;                    /* usage of the MRFs.                */
  printf("Same %d \n",same);
  printf("Similar %d \n",similar);
  printf("Different %d \n",different);
  printf("Only ml estimates: %d \n",mlestimates_only);
  printf("Ml estimates: %d \n",mlestimates);
  printf("Parameter updates: %d \n", em);

  iteration = 1;  /*Start the iterative algorithm */
  changed = TRUE;

  while(changed && (iteration < MAX_ITERATIONS)) {
    printf("Iteration %d \n",iteration);

    /* Update the likelihoods with current parameters only if that is necessary */
    if(em || iteration == 1) {
      for( i = 0; i < sizes[0]; ++i) {
        for( j = 0; j < sizes[1]; ++j) {
          for( k = 0; k < sizes[2]; ++k ) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) { 
              intensity[0] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              intensity[1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              intensity[2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              for(c = 1;c < PURE_CLASSES + 1;c++) { 
                AddMatrices(var[c],var_measurement,varsum);     
                val[c - 1] = Compute_Gaussian_likelihood3(intensity,mean[c],varsum );
              }
              val[WMGMLABEL - 1] = Compute_marginalized_likelihood3(intensity,mean[WMLABEL], mean[GMLABEL],
                                                               var[WMLABEL], var[GMLABEL], 
                                                               var_measurement, INTERVALS );
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
    for(i = 0; i < sizes[0]; ++i) {
      for(j = 0; j < sizes[1]; ++j) {
        for(k = 0; k < sizes[2]; ++k) {
          current_tissue_class = get_volume_real_value(volume_classified,i,j,k,0,0);
	  if(current_tissue_class != 0) {
            for(c = 0;c < CLASSES ;c++) {
              mrf_probability[c] = Compute_mrf_probability(c + 1,&volume_classified,i,j,k,
                                                             slice_width, beta, same , 
                                                             similar, different,pr_prior,sizes);
            }
            Normalize(mrf_probability,CLASSES);
            for(c = 0; c < CLASSES;c++) {
              mrf_probability[c] = get_volume_real_value(volume_likelihood[c],i,j,k,0,0) * 
                                 mrf_probability[c];
            }
            Normalize(mrf_probability,CLASSES);
            new_tissue_class = Maxarg(mrf_probability,CLASSES);
            if(new_tissue_class != current_tissue_class) {
              set_volume_real_value(volume_classified, i, j, k, 0 ,0, new_tissue_class);
              changed = TRUE; 
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

  filename = strcpy(filename,argv[6]); 
  output_modified_volume(strcat(filename,"_disc.mnc"),
                          NC_BYTE, FALSE, 0,6,volume_classified,argv[2],history,
			  (minc_output_options *) NULL); 

  
  if(!mlestimates_only) {
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_wm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[WMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_gm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[GMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_csf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[CSFLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
  }
  if(mlestimates) {
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_exactwm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[WMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_exactgm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[GMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[6]); 
    output_modified_volume(strcat(filename,"_exactcsf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[CSFLABEL - 1],argv[2],history,
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
const char POTTS_LOOKUP_TABLE[CLASSES + 1][CLASSES + 1] = {{0, 0, 0, 0, 0, 0, 1},
                                                           {0, 0, 0, 0, 1, 0, 1},
                                                           {0, 0, 0, 0, 1, 1, 0},
                                                           {0, 0, 0, 0, 0, 1, 0},
                                                           {0, 1, 1, 0, 0, 0, 0},
                                                           {0 ,0, 1, 1, 0, 0, 0},
                                                           {1, 1, 0, 0, 0, 0, 0}};

