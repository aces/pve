/* ----------------------------- MNI Header -----------------------------------
@NAME       :  pve
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: Partial volume calssification of the single channel MRI.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Jan 2002 (Jussi Tohka, jupeto@cs.tut.fi) 
@MODIFIED   : Last modification in May 2002 by Jussi Tohka
---------------------------------------------------------------------------- */

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

  double mean[PURE_CLASSES + 1]; /*   nuisance parameters */
  double var[PURE_CLASSES + 1]; 
  double var_measurement;                    /* measurent noise */
  double pr_prior;
#if defined(NOT_IMPL)
  double pr_wm,pr_gm,pr_csf,pr_wmgm,pr_gmcsf,pr_csfbg; /* For later use */
#endif /* NOT_IMPL defined */
  double beta = 0.1;

  double  slice_width[MAX_DIMENSIONS];

  double val[CLASSES],mrf_probability[CLASSES],value; /*  Some temporary variables */

  int same = -2;                          /* For defining the Potts model */
  int similar = -1;
  int different = 1;
 

  Volume volume_in,volume_mask; /* Volumes to be read from files: 
                                   in is the input image, and mask is the
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

  if(argc < 2) {
    printf("%s", ERROR_TOO_FEW);
    return(1);
  }

  if(!(strcmp(argv[1],"-help"))) {
    Display_help();
    return(0);
  }

  else if(!(strcmp(argv[1], "-file"))) {
    if(argc < 6) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Get_params_from_file(argv[5],mean,var,&var_measurement); 
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_FILE);
        printf("Subfunction Get_params_from_file returned errorcode %d .\n",error_code);
        return(2);
      }
      error_code = Open_images(argv[2],argv[3],&volume_in, &volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images returned errorcode %d .\n",error_code);
        return(3);
      }
    }
  }
  else if(!(strcmp(argv[1], "-image"))) {
    if(argc < 6) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Open_images(argv[2],argv[3],&volume_in, &volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images returned errorcode %d .\n",error_code);
        return(3);
      }
      error_code = Estimate_params_from_image(volume_in,volume_mask,argv[5],mean,var, 
                                             &var_measurement);
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_IMAGE);
        printf("Subfunction Estimate_params_from_image returned errorcode %d .\n",error_code);
        return(2);  
      }
    }
  }

  else if(!(strcmp(argv[1], "-tags"))) {
    if(argc < 6) {
      printf("%s", ERROR_TOO_FEW);
      return(1);
    }
    else {
      error_code = Open_images(argv[2],argv[3],&volume_in, &volume_mask);
      if( error_code != 0) {
        printf("\n %s \n",ERROR_INPUT_FILE);
        printf("Subfunction Open_images returned errorcode %d .\n",error_code);
        return(3);
      }
      if(strcmp(argv[5],"default")) {
        ptag_filename = strcpy(ptag_filename,argv[5]);
      }
      else {
        ptag_filename = strcpy(ptag_filename,DEFAULT_TAG_FILENAME);
      }
      error_code = Estimate_params_from_tags(tag_filename,volume_in,mean,var);
      var_measurement = 0;
      if(error_code != 0) {
        printf("\n %s \n",ERROR_PARAMS_TAG);
        printf("Subfunction Estimate_params_from_tag returned errorcode %d .\n",error_code);
        return(2);  
      }
    }
  }
  else {
    printf("%s",ERROR_SECOND_ARG);
    return(4);
  }

  /* Now see if there are more options given */

  if(argc > 6) {
    if(!(strcmp(argv[6], "-em"))) {
      em = TRUE;                  /* Turn on the parameter updates */
    }
  }
  if(argc > 7) {
    if(!(strcmp(argv[7], "-ml"))) {     /* Choices between simplified final PV estimation */
      mlestimates = TRUE;               /* and full ML PV estimation                      */
    }
    if(!(strcmp(argv[7], "-mlonly"))) {
      mlestimates_only = TRUE;
      mlestimates = TRUE;
    }
  } 
  if(argc > 12) {                       /* Tuning for Potts model  */
    if(!(strcmp(argv[8], "-mrf"))) {
      sscanf(argv[9],"%lf",&beta);
      sscanf(argv[10],"%d",&same);
      sscanf(argv[11],"%d",&similar);
      sscanf(argv[12],"%d",&different);
    }
  } 
  

  /* For producing necessary info to the outputfiles. */
  history = time_stamp(argc, argv);
  
  /* Initialize required volumes and set their ranges for 
     getting rid of unncessary surprises */
  for( c = 0;c < CLASSES;c++) {
    volume_likelihood[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);  
    set_volume_real_range( volume_likelihood[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX );
  }
 
  volume_classified = copy_volume_definition(volume_in, NC_BYTE, 
	  TRUE, 0 , 6); 
  set_volume_real_range( volume_classified, 0,CLASSES );

  get_volume_sizes( volume_in, sizes ); 

  get_volume_separations( volume_in , slice_width );  /* Get slice separations in each direction */ 
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

  iteration = 1;   /* Start the iterative algorithm */
  changed = TRUE;
 
  while(changed && (iteration < MAX_ITERATIONS)) {
    printf("Iteration %d \n",iteration);


  /* Update the likelihoods with current parameters only if that is necessary */
    if(em || iteration == 1) {
      printf("Computing likelihoods. \n");
      for( i = 0; i < sizes[0]; ++i) {
        for( j = 0; j < sizes[1]; ++j) {
          for( k = 0; k < sizes[2]; ++k ) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) { 
              value  = get_volume_real_value(volume_in,i,j,k,0,0);
              for(c = 1;c < PURE_CLASSES + 1;c++) { 
                val[c - 1] = Compute_Gaussian_likelihood(value,mean[c],
                                                       var[c] + var_measurement );
              }
              val[WMGMLABEL - 1] = Compute_marginalized_likelihood(value,mean[WMLABEL], mean[GMLABEL],
                                                               var[WMLABEL], var[GMLABEL], 
                                                               var_measurement, INTERVALS );
              val[GMCSFLABEL - 1] = Compute_marginalized_likelihood(value,mean[GMLABEL], mean[CSFLABEL],
                                                                var[GMLABEL], var[CSFLABEL], 
                                                                var_measurement, INTERVALS );
              val[CSFBGLABEL - 1] = Compute_marginalized_likelihood(value, mean[CSFLABEL], mean[BGLABEL], 
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
  

     /* ICM step */  
    printf("ICM step \n");
    changed = FALSE;
    for(i = 0; i < sizes[0]; ++i) {
      for(j = 0; j < sizes[1]; ++j) {
        for(k = 0; k < sizes[2]; ++k) {
          current_tissue_class = get_volume_real_value(volume_classified,i,j,k,0,0);
	  if(current_tissue_class != 0) {
            for( c = 0;c < CLASSES;c++) {
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
            if(em) {    /* Store necessary probalities for the parameter estimation step */
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
      volume_pve[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
  
    Compute_partial_volume_vectors(volume_in,volume_classified,volume_pve,mean);
  }
  if(mlestimates) {
    for(c = 0;c < PURE_CLASSES;c++) {
      volume_pve_ml[c] = copy_volume_definition(volume_in, NC_UNSPECIFIED, 
	  FALSE, 0.0 , 0.0);
      set_volume_real_range( volume_pve_ml[c],
                         LIKELIHOOD_RANGE_MIN , LIKELIHOOD_RANGE_MAX  );
    }
    Compute_partial_volume_vectors_ml(volume_in,volume_classified,volume_pve_ml,mean,var,var_measurement);
  }
  


  /* write necessary files */
  filename = strcpy(filename,argv[4]); 
  output_modified_volume(strcat(filename,"_disc.mnc"),
                          NC_BYTE, FALSE, 0,6,volume_classified,argv[2],history,
                         (minc_output_options *) NULL);
  if(!mlestimates_only) {
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_wm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[WMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_gm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[GMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_csf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve[CSFLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
  }
  if(mlestimates) {
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactwm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[WMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactgm.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[GMLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);
    filename = strcpy(filename,argv[4]); 
    output_modified_volume(strcat(filename,"_exactcsf.mnc"),
                          NC_UNSPECIFIED, FALSE, 0,0,volume_pve_ml[CSFLABEL - 1],argv[2],history,
                         (minc_output_options *) NULL);

  } 
 
  /* Free memory */

  delete_volume(volume_in);
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
