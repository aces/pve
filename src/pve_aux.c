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

int Estimate_params_from_image(Volume volume_in, Volume volume_mask, 
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
  for(c = 1; c < PURE_CLASSES + 1;c++) {
    if(!strcmp(ESTIMATOR,"ML")) {
      if(Estimate_ml(volume_in,volume_mask,volume_seg,c,&mean[c],&var[c]) != 0) 
        return(4);
    }
    else if(!strcmp(ESTIMATOR,"MVE")) {
      if(Estimate_mve(volume_in,volume_mask,volume_seg,c,&mean[c],&var[c]) != 0)
        return(4);
    }
    else {
      if(Estimate_mcd(volume_in,volume_mask,volume_seg,c,&mean[c],&var[c]) != 0)
        return(4);
    }
  }
  mean[BGLABEL] = 0;
  var[BGLABEL] = 0.1* MIN3(var[WMLABEL],var[GMLABEL],var[CSFLABEL]);
  *pvmeasurement = 0;
  printf("White matter mean: %f \n", mean[WMLABEL]);
  printf("Gray  matter mean: %f \n", mean[GMLABEL]);
  printf("CSF mean: %f \n", mean[CSFLABEL]);
  printf("White matter variance: %f \n", var[WMLABEL]);
  printf("Gray matter variance: %f \n", var[GMLABEL]);
  printf("CSF variance: %f \n", var[CSFLABEL]);
  delete_volume(volume_seg);
  return(0);
}

/* ML estimation of means and variances. Returns 0 if succesful. */

int Estimate_ml(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var)  
{ 
  double* samples = NULL;
  long int nofsamples;
  long int i;  

  samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples);
  if(samples == NULL) return(1);
  *mean = 0.0;
  for(i = 0;i < nofsamples;i++) {
    *mean = *mean + samples[i];
  }
  *mean = *mean/nofsamples; 
  *var = 0.0;
   for(i = 0;i < nofsamples;i++) {
    *var = *var + pow(samples[i] - *mean,2);
  }
   
  *var = *var/nofsamples;
  free(samples);
  return(0);
}

/* Mve estimation of the means and variances. Returns 0 if succesful. */

int Estimate_mve(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var)
{
  double* samples = NULL;
  long int nofsamples;
  int noftrials;
 
  samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples); 
  if(samples == NULL) return(1);
  noftrials = MIN(2*nofsamples,MAXTRIALS);
  if(minvarellipsoid(samples, nofsamples, noftrials,mean,var) != 0) return(2);
  return(0);
  
}  

/* MCD estimation of the means and variances. Returns 0 if succesful. */

int Estimate_mcd(Volume volume_in,Volume volume_mask,Volume volume_seg,char ref_label,
                 double* mean,double* var)
{
  double* samples = NULL;
  long int nofsamples;
 
 
  samples = Collect_values(volume_in,volume_mask,volume_seg,ref_label,&nofsamples); 
  if(samples == NULL) return(1);
  if(least_trimmed_squares(samples, nofsamples, mean,var) != 0) return(2);
  return(0);
}  


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

  volume_tmp = copy_volume_definition(volume_in, NC_BYTE, 
	  TRUE, 0 , 1); 
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
            if(NEIGHBOURHOOD == 26) { /* 26 - neighbourhood */
              on_the_border = FALSE;
              for(i1 = -1;i1 < 2;i1++) {
                for(j1 = -1;j1 < 2;j1++) {
                  for(k1 = -1;k1 < 2;k1++) {
                    if(rint(get_volume_real_value(volume_seg,i + i1,j + j1,k + k1,0,0)) 
                        != ref_label)  
                      on_the_border = TRUE;
		  }
                }
	      }
            }
            else if(NEIGHBOURHOOD == 6) { /* 6-neighbourhood */
              on_the_border = FALSE;
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
            else { on_the_border = FALSE; } /* no neighbourhood checking */
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
      sample_vector = malloc(MAXSAMPLES * sizeof(double));
      for(i =0;i < MAXSAMPLES;i++) {
        sampleno = rint(rand() *( (double) (MAXSAMPLES - 1) / RAND_MAX));
        sample_vector[i] = sample_vector_tmp[sampleno];
      }
      free(sample_vector_tmp);
      count = MAXSAMPLES;
    }
    else {
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
  }
  delete_volume(volume_tmp);
  *pcount = count;    
  return(sample_vector);
}

/* Get parameters based on a tag_file */

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

void Normalize(double* pval, char n)

{ double sca = 0.0;
  char i;

  for(i = 0;i < n;i++) {
    sca = sca + pval[i];
  }
  if(fabs(sca) >  VERY_SMALL) {       /* To avoid divisions by zero */
    for(i = 0;i < n;i++) {
      pval[i] = pval[i]/sca;
    }
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

Note: The numerical integration routine used by the 
function is primitive , but so is the mankind...

*/

double Compute_marginalized_likelihood(double value, double mean1 , double mean2, 
                                       double var1, double var2, 
                                       double measurement_var, 
                                       unsigned int nof_intervals)

{ 
  double lh, tmean , tvar, t, interval_len;
  int i;  
  
  interval_len = (double) 1 / nof_intervals;
  lh = 0;
  for(i = 0; i < nof_intervals; i++) {
    t = (i + 0.5) * interval_len;
    tmean = t * mean1 + ( 1 - t ) * mean2;
    tvar = pow(t,2) * var1 + pow((1 - t),2) * var2 + measurement_var;
    lh = lh + Compute_Gaussian_likelihood(value, tmean,tvar)*interval_len;
  }
  return(lh);
 }

/* Computes prior probalibility of a voxel have a certain label given the neighbourhood labels.
   Takes classified volume and three coordinates.
    The first argument is of course the label under investigation.
   prior is the prior probability of the label to occur.
   Integers same , similar and different and double beta specify the MRF.
   on_the_border and sizes are used detect when we are on image borders and prevent nasty
   segmentation faults.

 */

double Compute_mrf_probability(char label, Volume* pvolume, int x, int y , int z, 
                               double* slice_width, double beta, double same, double similar, double different, 
                               double prior, int* sizes)
{
  int i,j,k;
  char label2;  
  double exponent, distance;
  double similarity_value;
  char on_the_border;

  on_the_border = (x == 0) || (y == 0) || (z == 0) || (x == sizes[0] - 1) 
                  || (y == sizes[1] - 1) || (z == sizes[2] - 1); 
  /* To determine if it's possible to get out of image limits. 
     If not true (as it usually is) this saves the trouble calculating this 27 times */
 
  exponent = 0;
  for(i = -1; i < 2; i++) {
    for(j = -1; j < 2; j++) {
      for(k = -1; k < 2; k++) {
         if( i == 0 && j == 0 && k == 0 ) 
           exponent = exponent + prior;
         else {
           if(!on_the_border) { /* There's no danger of getting out of the image limits */
             label2 = get_volume_real_value(*pvolume, x + i, y + j , z + k,0,0);
           }
           else {
             if(x + i < 0 || y + j < 0 || z + k < 0 || x + i > sizes[0] - 1
                || y + j > sizes[1] - 1 || z + k > sizes[2] - 1) {
               label2 = BGLABEL;
             }
             else {
               label2 = get_volume_real_value(*pvolume, x + i, y + j , z + k,0,0);
             }
           }    
           if(Are_same(label,label2)) similarity_value = same;
           else if(Are_similar(label,label2)) similarity_value = similar;
           else similarity_value = different;

           distance = sqrt(pow(slice_width[0] * abs(i),2) +
                           pow(slice_width[1] * abs(j),2) +
                           pow(slice_width[2] * abs(k),2));

           exponent = exponent + similarity_value/distance;
	 }
      }
    }
  } 
  
  return(exp(- ( beta* exponent)));  
 

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
	      /*  if( k == 80 && j == 100 && i == 100) {
                 printf("Test: %f %f %f %f \n",total_probability,
                                             get_volume_real_value(probabilities[c - 1],i,j,k,0,0),
					     get_volume_real_value(volume_in,i,j,k,0,0),var[c]); } */
              total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
              var[c] = var[c] + get_volume_real_value(probabilities[c - 1],i,j,k,0,0)
                              * pow(get_volume_real_value(volume_in,i,j,k,0,0) - mean[c],2);
	      /*  if( k == 80 && j == 100 && i == 100) {
                 printf("Test: %f %f %f %f \n",total_probability,
                                             get_volume_real_value(probabilities[c - 1],i,j,k,0,0),
					     get_volume_real_value(volume_in,i,j,k,0,0),var[c]); } */
               
            }
          }
        }
      }
      var[c] = var[c] / total_probability;
      
    }
    printf("Variance, WM: %f GM: %f CSF: %f \n",var[WMLABEL],var[GMLABEL],var[CSFLABEL]);
    if(fabs(*var_measurement) > VERY_SMALL) {
      if(adapt_var_measurement) 
        *var_measurement = MIN3(var[1],var[2],var[3])*MEASUREMENT_FRACTION;
      for(c = 1;c < 4;c++) {
        var[c] = var[c] - *var_measurement;
        if(var[c] < 0) var[c] = 0; /* To ensure that variances are not negative */ 
      }
    }
  }
  else {    /* The model contains only measurement noise. We update only
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
    printf("Variance: %f \n",*var_measurement);
  }
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
             break;
        case WMLABEL: 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             break;
        case GMLABEL: 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             break;                
        case CSFLABEL: 
	     set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1.0);
             break;
        case WMGMLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[GMLABEL]) / 
                                       (mean[WMLABEL] - mean[GMLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case GMCSFLABEL: 
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[CSFLABEL]) 
                                     / (mean[GMLABEL] - mean[CSFLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0, t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1 - t);
             break; 
        case CSFBGLABEL:
             t = (get_volume_real_value(volume_in,i,j,k,0,0) - mean[BGLABEL]) / 
                                       (mean[CSFLABEL] - mean[BGLABEL]);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,t); 
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
             break;
        case WMLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             break;
        case GMLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,1.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,0.0);
             break;                
        case CSFLABEL: 
             set_volume_real_value(volume_pve_ml[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve_ml[CSFLABEL - 1],i,j,k,0,0,1.0);
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
             break;
        default: return(1); break;
        }
      }
    }
  }
  return(0);
}
