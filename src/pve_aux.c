/* pve_aux.c: auxilary functions for pve.c: */

#include "pve_aux.h"
#include "minvarellipsoid.h"
#include "least_trimmed_squares.h"

/* -----------------------------------------------------
   Functions related to output. Help, errormessages , etc. 
   ----------------------------------------------------- */

/* Displys some help about the usage of the program */

void Display_help() 
{  
printf("\n Usage 1: \n");
printf(" pve  -file input_file.mnc brainmask.mnc output_file_name parametersfile.pve \n");
printf(" Give output_filename without mnc extension. \n");
printf(" Look at file pve_aux.c to get an idea how parameters file should look out. \n");
printf(" \n");
printf(" Usage 2: \n");
printf(" pve -image input_file.mnc brainmask.mnc output_filename segmented_image.mnc \\ \n");
printf(" Again give output_filename without mnc extension. \n");
printf(" segmented_image is hard segmentation from input_file.mnc created by e.g. \n");
printf(" classify_clean \n");
printf(" \n");
printf(" Usage 3:\n");
printf(" pve -tags input_file.mnc brainmask.mnc output_filename tagfilename.tag \n");
printf(" Again give output_filename without mnc extension. \n");
printf(" Instead of giving tagfile name you can say default in which case the default \n");
printf(" tagfile is used. \n");
printf(" MORE OPTIONS :\n");
printf(" [-em/-noem] [-ml/-mlonly/-noml] [-mrf beta same similar different ] \n");
printf(" Order of these is important!! \n");

}

void Usage_info(char* pname) {
  
  (void) fprintf(stderr,
		 "\nUsage: %s [<options>] -<np estimator> <estimator_file> -mask <maskfile> <infile> <outfile_prefix>\n", pname);
  (void) fprintf(stderr,"        (where -<np estimator> must be either -file, -image or -tags)\n\n");
  (void) fprintf(stderr,"       %s [-help]\n\n", pname);

}

/* -----------------------------------------------------------
   Functions related to parameter input
-------------------------------------------------------------- */


/* --------------------------------------------------------------------- */
/* Gets paramaters from a ASCII file of the type:
pve
double repsesenting white matter mean
 "                  gray     "
 "                  csf      "
 "                  background "
 "                  white matter variance
"                   gray    "
"                   csf     "
"                   background "
                    variance of measurent noise

Note that function does NOT check whether the input is reasonable. 
For example you can give a negative variance and the function does not notice it.
So be careful about the parameters you input... 

returns 0 if everything is ok.
*/

int Get_params_from_file(char* fn, double* mean , double* var, double* pvmeasurement)
{
  FILE *params_file;
  char tmp[255];
  
  params_file = fopen(fn,"r");
  if(params_file == NULL) return(1);

  /* Check that first line is "pve" */    
  if(fscanf(params_file,"%s",tmp) != 1) {
    fclose(params_file);
    return(2);
  }
 
  if(strcmp(tmp,"pve") != 0) { 
    fclose(params_file);
    return(3);
  }  

  /* get values */
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &mean[WMLABEL], &mean[GMLABEL], &mean[CSFLABEL], &mean[BGLABEL],
            &var[WMLABEL], &var[GMLABEL], &var[CSFLABEL], &var[BGLABEL],
            pvmeasurement) != 9) {
    fclose(params_file);
    return(2);
  }

  return(0);
}

/* Gets the parameters from the command line. This option is obsolete. */


/* Estimate parameters based on segmented image:
   Estimation is maximum_likelihood or minimum variance ellipsoid
   estimation based on trimmed image. i.e depending on the value of the constant
   NEIGHBOURHOOD , some voxels on the borders of tissue types may be exluded from
   estimation.  */

int Estimate_params_from_image(Volume volume_in, Volume volume_mask, Volume volume_subcort,
                               char* segmentation_filename, double* mean,
                               double* var, double* pvmeasurement) 
{
  Volume volume_seg;
  int sizes[MAX_DIMENSIONS];
  int sizes_seg[MAX_DIMENSIONS];
  char c;

  if(input_volume(segmentation_filename,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, &volume_seg, (minc_input_options *) NULL) != OK)
    return(1);

  if(get_volume_n_dimensions(volume_seg) != 3)
    return(2);
  get_volume_sizes( volume_in, sizes ); 
  get_volume_sizes( volume_seg, sizes_seg );  
  if(!((sizes[0] == sizes_seg[0]) && (sizes[1] == sizes_seg[1]) &&
      (sizes[2] == sizes_seg[2]))) return(3);
  
  /* Then start parameter estimation */ 

  for(c = 1; c <= PURE_CLASSES; c++) {
    if(!strcmp(ESTIMATOR,"ML")) {
      if(Estimate_ml(volume_in,volume_mask,volume_seg,volume_subcort,c,&mean[c],&var[c]) != 0) 
        return(4);
    } else if(!strcmp(ESTIMATOR,"MVE")) {
      if(Estimate_mve(volume_in,volume_mask,volume_seg,volume_subcort,c,&mean[c],&var[c]) != 0)
        return(4);
    } else {
      if(Estimate_mcd(volume_in,volume_mask,volume_seg, volume_subcort,c,&mean[c],&var[c]) != 0)
        return(4);
    }
  }
  mean[BGLABEL] = 0;
  var[BGLABEL] = 0.1* MIN3(var[WMLABEL],var[GMLABEL],var[CSFLABEL]);
  *pvmeasurement = 0;
  printf("White matter mean: %f \n", mean[WMLABEL]);
  printf("Gray  matter mean: %f \n", mean[GMLABEL]);
  printf("CSF mean: %f \n", mean[CSFLABEL]);
  if (volume_subcort != NULL)
    printf("SC mean: %f \n", mean[SCLABEL]);
  printf("White matter variance: %f \n", var[WMLABEL]);
  printf("Gray matter variance: %f \n", var[GMLABEL]);
  printf("CSF variance: %f \n", var[CSFLABEL]);
  if (volume_subcort != NULL)
    printf("SC variance: %f \n", var[SCLABEL]);
  delete_volume(volume_seg);
  return(0);
}

/* ML estimation of means and variances. Returns 0 if successful. */

int Estimate_ml(Volume volume_in,Volume volume_mask,Volume volume_seg,
                Volume volume_subcort, char ref_label,
                double* mean,double* var)  
{ 
  double* samples = NULL;
  long int nofsamples;
  long int i;  

  *mean = 0.0;
  *var = 0.0;

  if ( ref_label == SCLABEL ) {
    if( volume_subcort ) {
      samples = Collect_values_subcortical(volume_in,
                                           volume_subcort,ref_label,&nofsamples);
    } else {
      return(0);
    }
  } else {
    samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples);
  }
  if(samples == NULL) return(1);
  for(i = 0;i < nofsamples;i++) {
    *mean = *mean + samples[i];
  }
  *mean = *mean/nofsamples; 
   for(i = 0;i < nofsamples;i++) {
    *var = *var + pow(samples[i] - *mean,2);
  }
   
  *var = *var/nofsamples;
  free(samples);

  return(0);
}

/* Mve estimation of the means and variances. Returns 0 if successful. */

int Estimate_mve(Volume volume_in,Volume volume_mask,Volume volume_seg,
                 Volume volume_subcort, char ref_label, double* mean,double* var)
{
  double* samples = NULL;
  long int nofsamples;
  int noftrials;
 
  *mean = 0.0;
  *var = 0.0;

  if ( ref_label == SCLABEL ) {
    if( volume_subcort ) {
      samples = Collect_values_subcortical(volume_in, volume_subcort, ref_label, &nofsamples);
    } else {
      return(0);
    }
  } else {
    samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples); 
  }
  if(samples == NULL) return(1);
  noftrials = MIN(2*nofsamples,MAXTRIALS);
  if(minvarellipsoid(samples, nofsamples, noftrials,mean,var) != 0) {
    free(samples);
    return(2);
  }

  free(samples);
  return(0);
}  

/* MCD estimation of the means and variances. Returns 0 if succesful. */

int Estimate_mcd(Volume volume_in,Volume volume_mask,Volume volume_seg,
		 Volume volume_subcort, char ref_label, double* mean,double* var)
{
  double* samples = NULL;
  long int nofsamples;
 
  *mean = 0.0;
  *var = 0.0;

  if ( ref_label == SCLABEL ) {
    if( volume_subcort ) {
      samples = Collect_values_subcortical(volume_in, volume_subcort, ref_label, &nofsamples);
    } else {
      return(0);
    }
  } else {
    samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples); 
  }
  if(samples == NULL) return(1);
  if(least_trimmed_squares(samples, nofsamples, mean,var) != 0) {
    free(samples);
    return(2);
  }

  free(samples);
  return(0);
}

/* Function that collects intensity values for parameter estimation into one vector.
   Returns the pointer to the vector where the samples are collected.
   This function supports also MVE estimator */

double* Collect_values(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                long int* pcount) 
{
  int sizes[MAX_DIMENSIONS];
  long int i;
  int j,k,i1,j1,k1;
  long int count = 0;
  long int sampleno;
  BOOLEAN on_the_border;
  double *sample_vector_tmp;
  double *sample_vector;
  double voxel_value_interval;
  Volume volume_tmp; /* volume for storing information about voxels status i.e. 
                        is it to be included into the vector of intensity values from 
                        which the parameters are to be estimated */

  volume_tmp = copy_volume_definition(volume_in, NC_BYTE, TRUE, 0 , 1); 
  set_volume_real_range( volume_tmp, 0,1 );
  get_volume_sizes( volume_in, sizes ); 
  voxel_value_interval = (get_volume_real_max(volume_in) - get_volume_real_min(volume_in))
                         /DATATYPE_SIZE;
  srand(time(0)); /* Set seed */
  
  for(i = 1; i < sizes[0] - 1;i++) {
    for(j = 1;j < sizes[1] - 1;j++) {
      for( k = 1; k < sizes[2] - 1;k++) {
        set_volume_real_value(volume_tmp,i,j,k,0,0,0);
        if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {     /* Voxel in the brain */
          if(rint(get_volume_real_value(volume_seg,i,j,k,0,0)) == ref_label) {  /* And of right type */
	    /* Check neigbourhood: Choices for neighbourhood are 6 and 26 neignbourhood or 
               no neighbourhood checking. 
               The choice is defined by the constant NEIGHBOURHOOD  */
            on_the_border = FALSE;
            if(NEIGHBOURHOOD == 26) { /* 26 - neighbourhood */
              for(i1 = -1;i1 < 2;i1++) {
                for(j1 = -1;j1 < 2;j1++) {
                  for(k1 = -1;k1 < 2;k1++) {
                    if(rint(get_volume_real_value(volume_seg,i + i1,j + j1,k + k1,0,0)) 
                        != ref_label)  
                      on_the_border = TRUE;
		  }
                }
	      }
            } else if(NEIGHBOURHOOD == 6) { /* 6-neighbourhood */
              for(i1 = -1;i1 < 2;i1++) {
                if(rint(get_volume_real_value(volume_seg,i + i1,j,k,0,0)) 
                        != ref_label)  
                  on_the_border = TRUE;
              }
              for(i1 = -1;i1 < 2;i1++) {
                if(rint(get_volume_real_value(volume_seg,i,j + i1,k,0,0)) 
                        != ref_label)  
                  on_the_border = TRUE;
              }
              for(i1 = -1;i1 < 2;i1++) {
                if(rint(get_volume_real_value(volume_seg,i,j,k + i1,0,0)) 
                        != ref_label)  
                  on_the_border = TRUE;
              }
            }
            if(!on_the_border) {
              count = count + 1;
              set_volume_real_value(volume_tmp,i,j,k,0,0,1);
            }
          }
	}
      }
    }
  }
  
  if(count > 0) {
    if(!(strcmp(ESTIMATOR,"MVE") || strcmp(ESTIMATOR,"MCD") ) && (count > MAXSAMPLES)) {
    /* if we are using MVE or MCD estimator so 
     1) we must get rid of some samples;
     2) we must add little bit of noise to the other samples; */  
      sample_vector_tmp = malloc(count * sizeof (double));
      if(sample_vector_tmp == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_tmp,i,j,k,0,0) > 0.5) {
              sample_vector_tmp[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
              sampleno = sampleno + 1;
	    }
          }
        }
      }
      count = sampleno;
      sample_vector = malloc(MAXSAMPLES * sizeof(double));
      if( sample_vector == NULL ) return(NULL);
      for(i =0;i < MAXSAMPLES;i++) {
        sampleno = rint(rand() *( (double) (count - 1) / RAND_MAX));
        sample_vector[i] = sample_vector_tmp[sampleno];
      }
      free(sample_vector_tmp);
      count = MAXSAMPLES;
    } else {
      sample_vector = malloc(count * sizeof (double));
      if(sample_vector == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_tmp,i,j,k,0,0) > 0.5) {
              sample_vector[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
              sampleno = sampleno + 1;
	    }
          }
        }
      } 
    }
    if(!(strcmp(ESTIMATOR,"MVE") || strcmp(ESTIMATOR,"MCD") )) {
      for(i = 0; i < count; i++) {
        /* Add a bit of uniform random noise to the samples */
        sample_vector[i] = sample_vector[i] + (rand() - RAND_MAX / 2)*  
                                               voxel_value_interval/RAND_MAX;
    
      } 
    }           
  } else {
    sample_vector = NULL;
  }
  delete_volume(volume_tmp);
  *pcount = count;    
  return(sample_vector);
}


double* Collect_values_subcortical(Volume volume_in,Volume volume_subcort, char ref_label,
		       long int* pcount) {

  int sizes[MAX_DIMENSIONS];
  long int i;
  int j,k;
  long int count = 0;
  long int sampleno;
  double *sample_vector_tmp;
  double *sample_vector;
  double voxel_value_interval;

  voxel_value_interval = (get_volume_real_max(volume_in) - get_volume_real_min(volume_in))
    /DATATYPE_SIZE;
  srand(time(0)); /* Set seed */
  get_volume_sizes(volume_subcort,sizes);

  for(i = 1; i < sizes[0] - 1;i++) {
    for(j = 1;j < sizes[1] - 1;j++) {
      for( k = 1; k < sizes[2] - 1;k++) {
	if (get_volume_real_value(volume_subcort,i,j,k,0,0) > SC_TR )
	  count++;
      }
    }
  }

  if(count > 0) {
    if(!(strcmp(ESTIMATOR,"MVE") || strcmp(ESTIMATOR,"MCD") ) && (count > MAXSAMPLES)) {
    /* if we are using MVE or MCD estimator so 
      1) we must get rid of some samples;
      2) we must add little bit of noise to the other samples; */  
      sample_vector_tmp = malloc(count * sizeof (double));
      if(sample_vector_tmp == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
	for(j = 1;j < sizes[1] - 1;j++) {
	  for( k = 1; k < sizes[2] - 1;k++) {
	    if(get_volume_real_value(volume_subcort,i,j,k,0,0) > SC_TR ) {
	      sample_vector_tmp[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
	      sampleno = sampleno + 1;
	    }
	  }
	}
      }
      sample_vector = malloc(MAXSAMPLES * sizeof(double));
      if(sample_vector == NULL) return(NULL);
      for(i =0;i < MAXSAMPLES;i++) {
	sampleno = rint(rand() *( (double) (count - 1) / RAND_MAX));
	sample_vector[i] = sample_vector_tmp[sampleno];
      }
      free(sample_vector_tmp);
      count = MAXSAMPLES;
    } else {
      sample_vector = malloc(count * sizeof (double));
      if(sample_vector == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_subcort,i,j,k,0,0) > SC_TR ) {
              sample_vector[sampleno] = get_volume_real_value(volume_in,i,j,k,0,0);
              sampleno = sampleno + 1;
	    }
          }
        }
      }
    }

    if(!(strcmp(ESTIMATOR,"MVE") || strcmp(ESTIMATOR,"MCD") )) {
      for(i = 0; i < count; i++) {
        /* Add a bit of uniform random noise to the samples */
        sample_vector[i] = sample_vector[i] + (rand() - RAND_MAX / 2)*  
                                               voxel_value_interval/RAND_MAX;
      } 
    }
  } else {
    sample_vector = NULL;
  }

  *pcount = count;    
  return(sample_vector);
}

/* Get parameters based on a tag_file. */
/* Sub-cortical not available for this option. */

int Estimate_params_from_tags(char* tag_filename,Volume volume_in,
                              double* means, double* vars) 
{
   int n_tag_volumes, num_samples,i,c;
   int adj_num_samples[PURE_CLASSES];
   STRING *labels;
   double **tags;
   double wx,wy,wz;
   double v1,v2,v3;
   int sizes[MAX_DIMENSIONS];   

   get_volume_sizes( volume_in, sizes ); 

   if ( input_tag_file(tag_filename, &n_tag_volumes, &num_samples,
		      &tags, NULL, NULL, NULL, NULL, &labels ) != OK ) 
      return(1);

   for(i = 0;i < PURE_CLASSES; i++) {
     adj_num_samples[i] = 0;
     means[i + 1] = 0.0;
     vars[i + 1] = 0.0;
   }

   /* Compute means */
   for( i = 0;i < num_samples;i++) {
     
     /* get the world coordinate from the tag file */
     wx = tags[i][0];
     wy = tags[i][1];
     wz = tags[i][2];
    
    /* convert world into voxel coordinates */
     convert_3D_world_to_voxel(volume_in, wx, wy, wz, &v1, &v2, &v3);

     if( v1 > - 0.5 && v2 > -0.5 && v3 > -0.5 && v1 < sizes[0] - 0.5 && 
       v2 < sizes[1] - 0.5  && v3 < sizes[2] - 0.5) {
       c = atoi( labels[i]);
       if( c > 0 && c < PURE_CLASSES + 1) {
         means[c] = means[c] + 
                    get_volume_real_value(volume_in,rint(v1),rint(v2),rint(v3),0,0);
         adj_num_samples[c - 1] = adj_num_samples[c - 1] + 1;
       }
       else {
	 printf("\n Warning there is something curious with the labels \n");
       }
     }
   }
   for(i = 1;i < PURE_CLASSES + 1;i++) {
     means[i] = means[i]/adj_num_samples[i - 1];
   }
   /* Then same for variances */
    for( i = 0;i < num_samples;i++) {
     
     /* get the world coordinate from the tag file */
     wx = tags[i][0];
     wy = tags[i][1];
     wz = tags[i][2];
    
    /* convert world into voxel coordinates */
     convert_3D_world_to_voxel(volume_in, wx, wy, wz, &v1, &v2, &v3);

     if( v1 > - 0.5 && v2 > -0.5 && v3 > -0.5 && v1 < sizes[0] - 0.5 && 
       v2 < sizes[1] - 0.5  && v3 < sizes[2] - 0.5) {
       c = atoi( labels[i]);
       if( c > 0 && c < PURE_CLASSES + 1) {
         vars[c] = vars[c] + 
                    pow(get_volume_real_value(volume_in,rint(v1),rint(v2),rint(v3),0,0) 
                        - means[c],2);
       }
       else {
         printf("\n Warning there is something curious with the labels \n");
       }
     }
   }
   for(i = 1;i < PURE_CLASSES + 1;i++) {
     vars[i] = vars[i]/adj_num_samples[i - 1];
   }
   means[BGLABEL] = 0.0;  /* Set background mean  */
    /* Set background variance */
   vars[BGLABEL] = 0.1* MIN3(vars[WMLABEL],vars[GMLABEL],vars[CSFLABEL]);
      
   free_tag_points(n_tag_volumes, num_samples,
		  tags, NULL, NULL, NULL, NULL, labels );
   return(0);
}

/* -------------------------------------------------------------------- */
/* Opens the input image and mask image. Returns 0 if ok.              */

int Open_images(char* in_fn, char* mask_fn, Volume* pvolume_in, 
                Volume* pvolume_mask)
{
  if(input_volume(in_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_in, (minc_input_options *) NULL) != OK)
    return(1);
  if(input_volume(mask_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_mask, (minc_input_options *) NULL) != OK)
    return(2);

  /* Check that each image really has three dimensions */

  if(get_volume_n_dimensions(*pvolume_in) != 3)
    return(3);
  if(get_volume_n_dimensions(*pvolume_mask) != 3)
    return(4);

  /* If we have got this far everything is ok. */

  return(0);
}

/* -------------------------------------------------------------------
  Functions that do tiny tasks for other functions from this file.
  They are not called from main
  -------------------------------------------------------------------- */
/* No explanation needed ... */

/* Limits x to the range from 0 to 1. */ 

void Limit_0_1(double* x) 
{
  if( *x < 0) *x = 0.0;
  else if( *x > 1) *x = 1.0 ;

}

/* ---------------------------------------------------------------------------
   And finally functions that are called from main to do the actual computations.
  ---------------------------------------------------------------------------- */

/* Normalizes n probalities to sum to one */

int Normalize(double* pval, char n)

{ double sca = 0.0;
  char i;

  for(i = 0;i < n;i++) {
    sca = sca + pval[i];
  }
  if(fabs(sca) >  VERY_SMALL) {       /* To avoid divisions by zero */
    for(i = 0;i < n;i++) {
      pval[i] = pval[i]/sca;
    }
    return 0;
  } else {
    return 1;
  }
}


/* Finds maximum argument out of the n possibilities */

char Maxarg(double* pval,char n)
{
  double maximum,index;
  char i;
  
  maximum = pval[0];
  index = 1;
  for(i = 1;i < n;i++) {
    if(pval[i] > maximum) {
      index = i + 1;
      maximum = pval[i];
    }
  }
  return(index);
}


/* -------------------------------------------------------------------
Computes likelihood of value given parameters mean and variance. 
Returns the likelihood.                                                */

double Compute_Gaussian_likelihood(double value, double mean , double var)

{ 
  return(exp(-(pow((value - mean),2))/(2 * var))/(sqrt(2 * PI * var)));

}

/* -------------------------------------------------------------------
Computes the likelihoods for the mixed classes. Returns the likelihood.
var1,var2 are the variances of pdfs representing pure classes. measurement_var 
is the measurement noise. So the model for the variable y (representing the 
intensity value) that is composed of t * tissue1 and (1 - t)* tissue2 becomes :

y = t*x1 + (1 - t)*x2 + xm,
x1 ~ N(mean1,var1) , x2 ~ N(mean2,var2) , xm ~ N(0,measurement_var).

Note: Use Gauss-Legendre quadrature of order 7, with four intervals.

*/

double Compute_marginalized_likelihood(double value, double mean1 , double mean2, 
                                       double var1, double var2, 
                                       double measurement_var ) {

  double gl4_pt[] = { -0.861136311594, -0.339981043585,
                       0.339981043585,  0.861136311594 };
  double gl4_wgt[] = { 0.347854845137, 0.652145154863,
                       0.652145154863, 0.347854845137 };

  double lh, tmean , tvar, t, t1, t2, interval_len;
  int i, k;

  int nof_intervals = 4; 
  interval_len = (double) 1.0 / nof_intervals;
  lh = 0.0;
  t1 = 0.0;
  for(i = 0; i < nof_intervals; i++) {
    t2 = t1 + interval_len;
    for( k = 0; k < 4; k++ ) {
      t = 0.5 * ( interval_len * gl4_pt[k] + t1 + t2 );
      tmean = mean2 + t * ( mean1 - mean2 );
      tvar = measurement_var + var2 + t * ( -var2 - var2 + t * ( var1 + var2 ) );
      lh += gl4_wgt[k] * Compute_Gaussian_likelihood( value, tmean, tvar );
    }
    t1 += interval_len;
  }
  lh *= 0.5 * interval_len;

  return(lh);
}

/* Computes prior probalibility of a voxel have a certain label given the neighbourhood labels.
   Takes classified volume and three coordinates.
    The first argument is of course the label under investigation.
   prior is the prior probability of the label to occur.
   Integers same , similar and different and double beta specify the MRF.
   No need to check for on_the_border as this is checked in the likelyhood loop.
 */

void Compute_mrf_probability(double* mrf_probability, Volume* pvolume, int x, int y , int z, 
                             double* width_stencil, double beta, double same, double similar, 
                             double different, double prior ) {

  int i,j,k,ii;
  char label, label2;  
  double similarity_value;

  for( label = 0; label < CLASSES; label++ ) {
    mrf_probability[label] = 0;
  }

  ii = 0; 
  for(i = -1; i < 2; i++) {
    for(j = -1; j < 2; j++) {
      for(k = -1; k < 2; k++) {
        if( i == 0 && j == 0 && k == 0 ) {
          for( label = 0; label < CLASSES; label++ ) {
            mrf_probability[label] += prior;
          }
        } else {
          label2 = get_volume_real_value(*pvolume, x + i, y + j , z + k,0,0);
          for( label = 0; label < CLASSES; label++ ) {
            if(Are_same(label+1,label2)) {
              similarity_value = same;
            } else if(Are_similar(label+1,label2)) {
              if (((label+1) == GMCSFLABEL)||((label+1) == CSFLABEL)) {
                similarity_value = similar;
              } else {
                similarity_value = -1;
              }
            } else {
              similarity_value = different;
            }
            mrf_probability[label] += similarity_value * width_stencil[ii];
          }
        }
        ii++;
      }
    }
  } 
  for( label = 0; label < CLASSES; label++ ) {
    mrf_probability[label] = exp( -beta * mrf_probability[label] );
  }

} 



/* Function that handles the re-estimation of means and variances for pure 
   tissue classes. Currently the variance of measurement noise is estimated 
   to be a fraction - defined by MEASUREMENT_FRACTION (see the header file) -
   of the lowest tissue class variance. HOWEVER IF THE VARIANCE OF THE MEASUREMENT
   NOISE WAS INITIALLY ZERO, THEN IT IS ASSUMED THAT THE MODEL WHICH DOES NOT 
   CONTAIN MEASUREMENT NOISE COMPONENT IS USED AND HENCE VARIANCE OF MEASUREMENT
   IS ZERO. This is implemented in a way that if parameter var_measurement is zero, it 
   is treated as a constant. Also if all tissue class variances are zero, then it is 
   assumed that the model with only measurement noise component is used and tissue class 
   variances remain zero. Also parameters related to the backgraound class are not updated.
*/

void Parameter_estimation(Volume volume_in, Volume volume_mask, 
                         Volume probabilities[CLASSES],double* mean, double* var,
                         double* var_measurement)
{
  int i,j,k,sizes[MAX_DIMENSIONS];
  char c;
  double total_probability;
  int adapt_var_measurement = 0;

  get_volume_sizes(volume_in,sizes);
 
  /* Print old parameters for the reference */
  printf("Mean, WM: %f GM: %f CSF: %f \n",mean[WMLABEL],mean[GMLABEL],mean[CSFLABEL]);
  printf("Variance, WM: %f GM: %f CSF: %f \n",var[WMLABEL],var[GMLABEL],var[CSFLABEL]);
  /* First re-estimate means */
 
  for(c = 1;c < PURE_CLASSES + 1;c++) {
    total_probability = 0;
    mean[c] = 0;
    for( i = 0; i < sizes[0];i++) {
      for( j = 0;j < sizes[1];j++) {
        for( k = 0; k < sizes[2]; k++) {
          if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
            total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
            mean[c] = mean[c] + get_volume_real_value(probabilities[c - 1],i,j,k,0,0)
                              * get_volume_real_value(volume_in,i,j,k,0,0);
          }
	}
      }
    }
    mean[c] = mean[c] / total_probability;
  }

  printf("Mean, WM: %f GM: %f CSF: %f \n",mean[WMLABEL],mean[GMLABEL],mean[CSFLABEL]);

  /* Then re-estimate the variances.
     Case 1: If the model contains physiological noise components. */
  if((fabs(var[1]) > VERY_SMALL) || (fabs(var[2]) > VERY_SMALL) || /* Noise model contains tissue-wise */
     (fabs(var[3]) > VERY_SMALL)) {                                /* components.                      */
    for(c = 1;c < 4;c++) {
      total_probability = 0;
      var[c] = 0;
      for( i = 0; i < sizes[0];i++) {
        for( j = 0;j < sizes[1];j++) {
          for( k = 0; k < sizes[2]; k++) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
              total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
              var[c] = var[c] + get_volume_real_value(probabilities[c - 1],i,j,k,0,0)
                              * pow(get_volume_real_value(volume_in,i,j,k,0,0) - mean[c],2);
            }
          }
        }
      }
      var[c] = var[c] / total_probability;
      
    }
    if(fabs(*var_measurement) > VERY_SMALL) {
      if(adapt_var_measurement) 
        *var_measurement = MIN3(var[1],var[2],var[3])*MEASUREMENT_FRACTION;
      for(c = 1;c < 4;c++) {
        var[c] = var[c] - *var_measurement;
        if(var[c] < 0) var[c] = 0; /* To ensure that variances are not negative */ 
      }
    }
    printf("Variance, WM: %f GM: %f CSF: %f \n",var[WMLABEL],var[GMLABEL],var[CSFLABEL]);
  } else {    /* The model contains only measurement noise. We update only
               measurement variance.          */
    *var_measurement = 0;
    total_probability = 0;
    for(c = 1;c < 4;c++) {
      for( i = 0; i < sizes[0];i++) {
        for( j = 0;j < sizes[1];j++) {
          for( k = 0; k < sizes[2]; k++) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
              total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
              *var_measurement = *var_measurement 
                                 + get_volume_real_value(probabilities[c - 1],i,j,k,0,0)
                                 * pow(get_volume_real_value(volume_in,i,j,k,0,0) - mean[c],2);  
            }
          }
        }
      }
    }
    *var_measurement = *var_measurement / total_probability;
  }
}

/* Function that re-estimates the means and variances for pure tissue classes
   based on the current classified image and partial estimates. Self-correcting
   the means and variances improves tissue classification.
*/

int Parameter_estimation_classified(Volume volume_in, Volume volume_mask, 
                                    Volume volume_subcort, Volume classified,
                                    double * mean, double * var,
                                    double * var_measurement ) {

  int i,j,k,count,sizes[MAX_DIMENSIONS];
  char c, label;
  double total_probability;
  double t;  /* mixing proportion to be estimated */
  int adapt_var_measurement = 0;
  Volume pve_cls;

  get_volume_sizes(volume_in,sizes);
  pve_cls = copy_volume_definition(volume_in, NC_BYTE, TRUE, 0 , PURE_CLASSES); 
  set_volume_real_range( pve_cls, 0, PURE_CLASSES );
 
  /* Compute new classification based on partial estimates (simplified model). */

  for( i = 0; i < sizes[0];i++) {
    for( j = 0;j < sizes[1];j++) {
      for( k = 0; k < sizes[2]; k++) {
        c = get_volume_real_value(classified,i,j,k,0,0);
	label = BGLABEL;
	switch(c) {
	case BGLABEL:
        case WMLABEL:
        case GMLABEL:
        case CSFLABEL:
	case SCLABEL: label = c; break;
        case WMGMLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[GMLABEL]) / 
                                       (mean[WMLABEL] - mean[GMLABEL]);
             Limit_0_1(&t);
             if( t > 0.5 ) {
               label = WMLABEL;
             } else {
               label = GMLABEL;
             }
             break; 
        case GMCSFLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[CSFLABEL]) 
                                     / (mean[GMLABEL] - mean[CSFLABEL]);
             Limit_0_1(&t);
             if( t > 0.5 ) {
               label = GMLABEL;
             } else {
               label = CSFLABEL;
             }
             break; 
        case CSFBGLABEL:
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[BGLABEL]) / 
                                       (mean[CSFLABEL] - mean[BGLABEL]);
             Limit_0_1(&t);
             if( t > 0.5 ) {
               label = CSFLABEL;
             }
             break;
	case WMSCLABEL:
	     t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[SCLABEL]) / 
                                       (mean[WMLABEL] - mean[SCLABEL]);
	     Limit_0_1(&t);
             if( t > 0.5 ) {
               label = WMLABEL;
             } else {
               label = SCLABEL;
             }
             break;
	case SCGMLABEL:
	     t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[GMLABEL]) / 
                                       (mean[SCLABEL] - mean[GMLABEL]);
	     Limit_0_1(&t);
             if( t > 0.5 ) {
               label = SCLABEL;
             } else {
               label = GMLABEL;
             }
             break;	    
        default: break;
        }
        set_volume_real_value(pve_cls,i,j,k,0,0,label);
      }
    }
  }

  /* Print old parameters for the reference */
  printf("Old Mean, WM: %f GM: %f CSF: %f\n",
         mean[WMLABEL],mean[GMLABEL],mean[CSFLABEL]);
  printf("Old Variance, WM: %f GM: %f CSF: %f\n",
         var[WMLABEL],var[GMLABEL],var[CSFLABEL]);

  double old_mean[PURE_CLASSES + 1]; /*   old nuisance parameters */
  double old_var[PURE_CLASSES + 1];
  for( c = 1; c <= PURE_CLASSES; c++ ) {
    old_mean[c] = mean[c];
    old_var[c] = var[c];
  }

  /* Then start parameter estimation */ 
  for( c = 1; c <= PURE_CLASSES; c++ ) {
    if(!strcmp(ESTIMATOR,"ML")) {
      if( Estimate_ml( volume_in, volume_mask, pve_cls, volume_subcort, c, &mean[c], &var[c] ) != 0)
        return(4);
    } else if(!strcmp(ESTIMATOR,"MVE")) {
      if( Estimate_mve( volume_in, volume_mask, pve_cls, volume_subcort, c, &mean[c], &var[c] ) != 0)
        return(4);
    } else {
      if( Estimate_mcd( volume_in, volume_mask, pve_cls, volume_subcort, c,
                        &mean[c], &var[c] ) != 0)
        return(4);
    }
  }
  mean[BGLABEL] = 0;
  var[BGLABEL] = 0.1* MIN3(var[WMLABEL],var[GMLABEL],var[CSFLABEL]);

  printf("New Mean, WM: %f GM: %f CSF: %f\n",
         mean[WMLABEL],mean[GMLABEL],mean[CSFLABEL]);
  printf("New Variance, WM: %f GM: %f CSF: %f\n",
         var[WMLABEL],var[GMLABEL]);

  double err_mean = 0.0;
  double err_var = 0.0;
  for( c = 1; c < PURE_CLASSES; c++ ) {
    err_mean += ( old_mean[c] - mean[c] ) * ( old_mean[c] - mean[c] ) / 
                ( 1.0e-10 + mean[c] * mean[c] );
    err_var += ( old_var[c] - var[c] ) * ( old_var[c] - var[c] ) / 
               ( 1.0e-10 + var[c] * var[c] );
  }
  err_mean = sqrt( err_mean );
  err_var = sqrt( err_var );
  printf("Change in mean = %f  Change in variance = %f\n", err_mean, err_var );

  *var_measurement = 0;

  delete_volume(pve_cls);

  return( ( err_mean > 0.001 ) || ( err_var > 0.005 ) );
}

/* And finally, a function that calculates the partial volume vectors:
   Made for whole images (as oppesed to everything else) with the hope that 
   someone someday invents a good way to use spatial information in this task 
   also. Note also that estimator used here 
   is just a heuristic (actually it's MLE for measurent noise only 
   model). However, with good enough data, it seems
   to approximate MLE well. MLE can be also computed by analytical means, but that results
   quite complicated expressions, so I don't think it's worthwhile at the moment. 

 */

int Compute_partial_volume_vectors(Volume volume_in,Volume volume_classified,
                                   Volume volume_pve[PURE_CLASSES], double* mean)
{
  int sizes[MAX_DIMENSIONS];
  int i,j,k;
  char c;
  double t;  /* mixing proportion to be estimated */

  get_volume_sizes(volume_in,sizes);
   
  for( i = 0; i < sizes[0];i++) {
    for( j = 0;j < sizes[1];j++) {
      for( k = 0; k < sizes[2]; k++) {
        c = get_volume_real_value(volume_classified,i,j,k,0,0);
	switch(c) {
	case BGLABEL: 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case WMLABEL: 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case GMLABEL: 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;                
        case CSFLABEL: 
	     set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
	case SCLABEL:
	     set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,1.0);
	     break;
        case WMGMLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[GMLABEL]) / 
                                       (mean[WMLABEL] - mean[GMLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case GMCSFLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[CSFLABEL]) 
                                     / (mean[GMLABEL] - mean[CSFLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0, t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case CSFBGLABEL:
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[BGLABEL]) / 
                                       (mean[CSFLABEL] - mean[BGLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,t); 
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
	case WMSCLABEL:
	     t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[SCLABEL]) / 
                                       (mean[WMLABEL] - mean[SCLABEL]);
	     Limit_0_1(&t);
	     set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0); 
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,1 - t);	    
             break;
	case SCGMLABEL:
	     t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[GMLABEL]) / 
                                       (mean[SCLABEL] - mean[GMLABEL]);
	     Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0); 
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,t);	    
             break;	    
        default: return(1); break;
        }
      }
    }
  }
  return(0);
}
/* And finally, a function that calculates the partial volume vectors based on the exact 
   Maximum likelihood criterion. Uses grid search, but should not be too time taking. 
   However, is more computationally intensive than function above.

 */

int Compute_partial_volume_vectors_ml(Volume volume_in, Volume volume_classified,
                                      Volume volume_pve_ml[PURE_CLASSES],
                                      double* mean, double* var, double var_measurement)
{
  int sizes[MAX_DIMENSIONS];
  int i,j,k,n;
  char c;
  double t;  /* mixing proportion to be estimated */
  double prob;
  double value;
  double mean_tmp,var_tmp; 

  double maxvalue = - 1.0;
  int maxindex = 0;

  get_volume_sizes(volume_in,sizes);
   
 
  for( i = 0; i < sizes[0];i++) {
    for( j = 0;j < sizes[1];j++) {
      for( k = 0; k < sizes[2]; k++) {       
        c = get_volume_real_value(volume_classified,i,j,k,0,0);
	switch(c) {
	case BGLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case WMLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case GMLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
	     set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case CSFLABEL:
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
	case SCLABEL:
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,1.0);
             break;
        case WMGMLABEL:
             value = get_volume_real_value(volume_in,i,j,k,0,0);
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               mean_tmp = (((double) n)/100)*mean[WMLABEL] + (1 - ((double) n)/100)*mean[GMLABEL];
               var_tmp = pow(((double) n) / 100,2)*var[WMLABEL] + 
                         pow(1- ((double) n) / 100,2)*var[GMLABEL] + 
                         var_measurement;
               prob = Compute_Gaussian_likelihood(value,mean_tmp,var_tmp);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;    
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case GMCSFLABEL: 
             value = get_volume_real_value(volume_in,i,j,k,0,0);
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               mean_tmp = (((double) n)/100)*mean[GMLABEL] + (1 - ((double) n)/100)*mean[CSFLABEL];
               var_tmp = pow(((double) n) / 100,2)*var[GMLABEL] + 
                         pow(1- ((double) n) / 100,2)*var[CSFLABEL] + 
                         var_measurement;
               prob = Compute_Gaussian_likelihood(value,mean_tmp,var_tmp);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;    
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0, t);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case CSFBGLABEL:
             value = get_volume_real_value(volume_in,i,j,k,0,0);
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               mean_tmp = (((double) n)/100)*mean[CSFLABEL] + (1 - ((double) n)/100)*mean[BGLABEL];
               var_tmp = pow(((double) n) / 100,2)*var[CSFLABEL] + 
                         pow(1- ((double) n) / 100,2)*var[BGLABEL] + 
                         var_measurement;
               prob = Compute_Gaussian_likelihood(value,mean_tmp,var_tmp);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;    
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,t); 
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
	case WMSCLABEL:
             value = get_volume_real_value(volume_in,i,j,k,0,0);
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               mean_tmp = (((double) n)/100)*mean[WMLABEL] + (1 - ((double) n)/100)*mean[SCLABEL];
               var_tmp = pow(((double) n) / 100,2)*var[WMLABEL] + 
                         pow(1- ((double) n) / 100,2)*var[SCLABEL] + 
                         var_measurement;
               prob = Compute_Gaussian_likelihood(value,mean_tmp,var_tmp);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100; 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0); 
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,1 - t);	    
             break;
	case SCGMLABEL:
             value = get_volume_real_value(volume_in,i,j,k,0,0);
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               mean_tmp = (((double) n)/100)*mean[SCLABEL] + (1 - ((double) n)/100)*mean[GMLABEL];
               var_tmp = pow(((double) n) / 100,2)*var[SCLABEL] + 
                         pow(1- ((double) n) / 100,2)*var[GMLABEL] + 
                         var_measurement;
               prob = Compute_Gaussian_likelihood(value,mean_tmp,var_tmp);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;    
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0); 
             set_volume_real_value(volume_pve_ml[SCLABEL - 1],i,j,k,0,0,t);	    
             break;
        default: return(1); break;
        }
      }
    }
  }
  return(0);
}
