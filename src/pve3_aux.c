/* pve3_aux.c: auxilary functions for pve.c: */

#include "pve3_aux.h"

/* -----------------------------------------------------
   Functions related to output. Help, errormessages , etc. 
   ----------------------------------------------------- */

/* Displys some help about the usage of the program */

void Display_help() 
{  
printf("\n Usage1: \n");
printf(" pve  -file infile_T1.mnc infile_T2.mnc infile_PD.mnc brainmask.mnc output_file_name parametersfile.pve \n");
printf(" Give output_filename without mnc extension. \n");
printf(" Look at file pve_aux.c to get an idea how parameters file should look out. \n");
printf(" \n");
printf(" Usage2: \n");
printf(" pve -image infile_T1.mnc infile_T2.mnc infile_PD.mnc brainmask.mnc output_filename segmented_image.mnc \\ \n");
printf(" Again give output_filename without mnc extension. \n");
printf(" segmented_image is hard segmentation from input_file.mnc created by e.g. \n");
printf(" classify_clean \n");
printf(" \n");
printf(" Usage 3:\n");
printf(" pve -tags infile_T1.mnc  infile_T2.mnc infile_PD.mnc brainmask.mnc output_filename tagfilename.tag \n");
printf(" Again give output_filename without mnc extension. \n");
printf(" Instead of giving tagfile name you can say default in which case the default \n");
printf(" tagfile is used. \n");
printf(" MORE OPTIONS :\n");
printf(" [-em/-noem] [-ml/-mlonly/-noml] [-mrf beta same similar different ] \n");
printf(" Order of these is important \n");

}

void Usage_info(char* pname) {
  
  (void) fprintf(stderr,
		 "\nUsage: %s [<options>] -<np estimator> <estimator_file> -mask <maskfile> <t1> <t2> <pd> <outfile_prefix>\n", pname);
  (void) fprintf(stderr,"        (where -<np estimator> must be either -file, -image or -tags)\n\n");
  (void) fprintf(stderr,"       %s [-help]\n\n", pname);

}



/* -----------------------------------------------------------
   Functions related to parameter input
-------------------------------------------------------------- */

/* Checks that input matrices (five of them) are proper. i.e. that they are 
   symmetric and positive definite or zero. All matrices cannot be zero so 
   least one matrix must be non-zero. However, you can give zero matrix for 
   white matter variance and non-zero for gray matter variance because in 
   this way it is possible to remove a class from consideration. */

int Are_proper(pMatrix var_wm,pMatrix var_gm,pMatrix var_csf, 
               pMatrix var_bg,pMatrix var_measurement)
{
  int  return_value = 1;

  if((!(IsPositiveDefinite(var_wm) || IsZero(var_wm))) || 
      (!IsSymmetric(var_wm)))
    return_value = 0;
  if((!(IsPositiveDefinite(var_gm) || IsZero(var_gm))) || 
      (!IsSymmetric(var_gm)))
    return_value = 0;
  if((!(IsPositiveDefinite(var_csf) || IsZero(var_csf))) || 
      (!IsSymmetric(var_csf)))
    return_value = 0;
  if((!(IsPositiveDefinite(var_bg) || IsZero(var_bg))) || 
      (!IsSymmetric(var_bg)))
    return_value = 0;
  if((!(IsPositiveDefinite(var_measurement) || IsZero(var_measurement))) || 
      (!IsSymmetric(var_measurement)))
    return_value = 0;  
  if(IsZero(var_wm) && IsZero(var_gm) && IsZero(var_csf) &&
     IsZero(var_bg) && IsZero(var_measurement))
    return_value = 0;

  return(return_value);
} 


/* --------------------------------------------------------------------- */
/* Gets paramaters from an ASCII file of the type:
pve3
3 doubles repsesenting white matter mean
 "                  gray     "
 "                  csf      "
 "                  background "
9 doubles representing white matter variance
"                   gray    "
"                   csf     "
"                   background "
"                   variance of measurent noise



returns 0 if everything is ok.
*/

int Get_params_from_file3(char* fn, pVector means[PURE_CLASSES + 1],
                          pMatrix vars[PURE_CLASSES + 1], pMatrix pvmeasurement)
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
 
  if(strcmp(tmp,"pve3") != 0) { 
    fclose(params_file);
    return(3);
  }  

  /* get values */
  if(fscanf(params_file,"%lf%lf%lf", &means[WMLABEL][0],&means[WMLABEL][1],&means[WMLABEL][2]) != 3) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf", &means[GMLABEL][0],&means[GMLABEL][1],&means[GMLABEL][2]) != 3) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf", &means[CSFLABEL][0],&means[CSFLABEL][1],&means[CSFLABEL][2]) != 3) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf", &means[BGLABEL][0],&means[BGLABEL][1],&means[BGLABEL][2]) != 3) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &vars[WMLABEL][0],&vars[WMLABEL][1],&vars[WMLABEL][2],&vars[WMLABEL][3],
            &vars[WMLABEL][4],&vars[WMLABEL][5],&vars[WMLABEL][6],&vars[WMLABEL][7],
            &vars[WMLABEL][8]) != 9) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &vars[GMLABEL][0],&vars[GMLABEL][1],&vars[GMLABEL][2],&vars[GMLABEL][3],
            &vars[GMLABEL][4],&vars[GMLABEL][5],&vars[GMLABEL][6],&vars[GMLABEL][7],
            &vars[GMLABEL][8]) != 9) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &vars[CSFLABEL][0],&vars[CSFLABEL][1],&vars[CSFLABEL][2],&vars[CSFLABEL][3],&vars[CSFLABEL][4],
            &vars[CSFLABEL][5],&vars[CSFLABEL][6],&vars[CSFLABEL][7],&vars[CSFLABEL][8]) != 9) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &vars[BGLABEL][0],&vars[BGLABEL][1],&vars[BGLABEL][2],&vars[BGLABEL][3],&vars[BGLABEL][4],
            &vars[BGLABEL][5],&vars[BGLABEL][6],&vars[BGLABEL][7],&vars[BGLABEL][8]) != 9) {
    fclose(params_file);
    return(2);
  }
  if(fscanf(params_file,"%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
            &pvmeasurement[0],&pvmeasurement[1],&pvmeasurement[2],
            &pvmeasurement[3],&pvmeasurement[4],
            &pvmeasurement[5],&pvmeasurement[6],
            &pvmeasurement[7],&pvmeasurement[8]) != 9) {
    fclose(params_file);
    return(2);
  } 
  

  /* Check that all the covariance matrices are proper. */
  if(!Are_proper(vars[WMLABEL],vars[GMLABEL],vars[CSFLABEL],vars[BGLABEL],pvmeasurement))
    return(4);

  printf("WM mean");
  PrintVector(means[WMLABEL]);
  printf("GM mean");
  PrintVector(means[GMLABEL]);
  printf("CSF mean");
  PrintVector(means[CSFLABEL]);
  printf("WM var");
  PrintMatrix(vars[WMLABEL]);
  printf("GM var");
  PrintMatrix(vars[GMLABEL]);
  printf("CSF var");
  PrintMatrix(vars[CSFLABEL]);
 

  /* All seems to be ok so return 0 */
  return(0);
}

int Estimate_params_from_image3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
                                Volume volume_mask, Volume volume_subcort, 
				char* segmentation_filename, pVector mean[PURE_CLASSES + 1] ,
				pMatrix var[PURE_CLASSES + 1], pMatrix pvmeasurement) 
{
  Volume volume_seg;
  int sizes[MAX_DIMENSIONS];
  int sizes_seg[MAX_DIMENSIONS];
  char c,i;


  if(input_volume(segmentation_filename,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, &volume_seg, (minc_input_options *) NULL) != OK)
    return(1);


  if(get_volume_n_dimensions(volume_seg) != 3)
    return(2);
  get_volume_sizes( volume_inT1, sizes ); 
  get_volume_sizes( volume_seg, sizes_seg );  
  if(!((sizes[0] == sizes_seg[0]) && (sizes[1] == sizes_seg[1]) &&
      (sizes[2] == sizes_seg[2]))) return(3);
  
  /* Then start parameter estimation */ 


  for(c = 1; c < PURE_CLASSES + 1;c++) {
    if(!strcmp(ESTIMATOR,"ML")) {

      if ((c == SCLABEL)&&(volume_subcort == NULL))
	continue;

      if (Estimate_ml3(volume_inT1,volume_inT2,volume_inPD,volume_mask,
		       volume_seg,volume_subcort,c,mean[c],var[c]) != 0) return(3);

    }   
    else {
      if(Estimate_mve3(volume_inT1,volume_inT2,volume_inPD,volume_mask,
                    volume_seg,c,mean[c],var[c]) != 0) return(3);
    }
  }

  
  SetZeroVector(mean[BGLABEL]);
  SetZero(var[BGLABEL]);
  for(i = 1;i < 4;i++) {
    SetElement(var[BGLABEL],i,i, 0.1* 
               MIN3(GetElement(var[WMLABEL],i,i),
                    GetElement(var[GMLABEL],i,i),
                    GetElement(var[CSFLABEL],i,i)));
  } 

  SetZero(pvmeasurement);
  printf("White matter mean:");
  PrintVector(mean[WMLABEL]);
  printf("Gray  matter mean:");
  PrintVector(mean[GMLABEL]);
  printf("CSF mean:");
  PrintVector(mean[CSFLABEL]);
  if (volume_subcort != NULL) {
    printf("Subcortical mean:");
    PrintVector(mean[SCLABEL]);
  }
  printf("White matter variance:");
  PrintMatrix(var[WMLABEL]);
  printf("Gray matter variance:" );
  PrintMatrix(var[GMLABEL]);
  printf("CSF variance:");
  PrintMatrix(var[CSFLABEL]);
  if (volume_subcort != NULL) {
    printf("Subcortical variance:");
    PrintMatrix(var[SCLABEL]);
  }
  delete_volume(volume_seg);
  return(0);
}

/* ML estimation of means and variances */

int Estimate_ml3(Volume volume_inT1, Volume volume_inT2,Volume volume_inPD,
		 Volume volume_mask,Volume volume_seg,Volume volume_subcort,
		 char ref_label,pVector mean,pMatrix var)  
{ 
  double* samples = NULL;
  long int nofsamples;
  long int i,j;  
  Vector3D standardized_sample;
  Matrix3D tmpcov;

  if (ref_label == SCLABEL) {
    samples = Collect_values3_subcortical(volume_inT1,volume_inT2,volume_inPD,
			      volume_subcort,ref_label,&nofsamples);   
  }
  else {
    samples = Collect_values3(volume_inT1,volume_inT2,volume_inPD,
			      volume_mask,volume_seg,ref_label,&nofsamples);
  }
  if(samples == NULL) return(1);
  SetZeroVector(mean);
  for(i = 0;i < nofsamples;i++) {
    for(j = 0;j < 3;j++) {
      mean[j] = mean[j] + samples[3*i + j];
    }
  }

  ScalarMultiplyVector(mean,(double) 1/nofsamples, mean); 
  SetZero(var);
  for(i = 0;i < nofsamples;i++) {
    for(j = 0;j < 3;j++) { 
      standardized_sample[j] = samples[3*i + j] - mean[j];
    }
    VectorProduct(standardized_sample,standardized_sample,tmpcov);
    AddMatrices(var,tmpcov,var); 
  }
   
  ScalarMultiply(var,(double) 1 / nofsamples, var);
  free(samples);

  return(0);


}

/* Mve estimation of the means and variances. */

int Estimate_mve3(Volume volume_inT1, Volume volume_inT2,Volume volume_inPD,
                 Volume volume_mask,Volume volume_seg,char ref_label,
                  pVector mean,pMatrix var)
{
  double* samples = NULL;
  long int nofsamples;
  int noftrials;
 
  samples = Collect_values3(volume_inT1,volume_inT2,volume_inPD,
                            volume_mask,volume_seg,ref_label,&nofsamples);
  if(samples == NULL) return(1); 
  noftrials = MIN(2*nofsamples,MAXTRIALS);
  if(minvarellipsoid3(samples, nofsamples, noftrials,mean,var) != 0) return(1);
  return(0);
}  


/* Function that collects intensity values for parameter estimation into one vector.
   Returns the pointer to the vector where the samples are collected.
   This function supports also MVE estimator */

double* Collect_values3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD,
                        Volume volume_mask,Volume volume_seg,char ref_label,
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
  double voxel_value_interval[3];
  Volume volume_tmp; /* volume for storing information about voxels status i.e. 
                        is it to be included into the vector of intensity values from 
                        which the parameters are to be estimated */

  volume_tmp = copy_volume_definition(volume_inT1, NC_BYTE, 
	  TRUE, 0 , 1); 
  set_volume_real_range( volume_tmp, 0,1 );
  get_volume_sizes( volume_inT1, sizes ); 
  voxel_value_interval[0] = (get_volume_real_max(volume_inT1) - 
                             get_volume_real_min(volume_inT1))/DATATYPE_SIZE;
  voxel_value_interval[1] = (get_volume_real_max(volume_inT2) - 
                             get_volume_real_min(volume_inT2))/DATATYPE_SIZE;
  voxel_value_interval[2] = (get_volume_real_max(volume_inPD) - 
                             get_volume_real_min(volume_inPD))/DATATYPE_SIZE;
  srand(time(0)); /* Set seed */
  
  for(i = 1; i < sizes[0] - 1;i++) {
    for(j = 1;j < sizes[1] - 1;j++) {
      for( k = 1; k < sizes[2] - 1;k++) {
        set_volume_real_value(volume_tmp,i,j,k,0,0,0);
        if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {     /* Voxel in the brain */
          if(rint(get_volume_real_value(volume_seg,i,j,k,0,0)) == ref_label) {  /* And of right type */
	    /* Check neigbourhood: Choices for neighbourhood are 6 and 26 neighbourhood or 
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
    /* if we are using MVE or MCD (not implemented) estimator so    
     1) we must get rid of some samples;
     2) we must add little bit of noise to the other samples; */  
      sample_vector_tmp = malloc(count * 3 * sizeof (double)); /* reserve room for 3 values (T1,T2,PD) for each sample */
      if(sample_vector_tmp == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_tmp,i,j,k,0,0) > 0.5) {
              sample_vector_tmp[sampleno] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              sample_vector_tmp[sampleno + 1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              sample_vector_tmp[sampleno + 2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              sampleno = sampleno + 3;
	    }
          }
        }
      }
      sample_vector = malloc(MAXSAMPLES * 3 * sizeof(double));
      if(sample_vector == NULL) return(NULL);
      for(i =0;i < MAXSAMPLES;i++) {
        sampleno = rint(rand() *( (double) (MAXSAMPLES - 1) / RAND_MAX));
        sample_vector[3 * i] = sample_vector_tmp[3 * sampleno];
        sample_vector[3 * i + 1] = sample_vector_tmp[3 * sampleno + 1];
        sample_vector[3 * i + 2] = sample_vector_tmp[3 * sampleno + 2];
      }
      free(sample_vector_tmp);
      count = MAXSAMPLES;
    }
    else {
      sample_vector = malloc(3 * count * sizeof (double));
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_tmp,i,j,k,0,0) > 0.5) {
              sample_vector[sampleno] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              sample_vector[sampleno + 1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              sample_vector[sampleno + 2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              sampleno = sampleno + 3;              
	    }
          }
        }
      } 
    }
    if(!(strcmp(ESTIMATOR,"MVE")  || strcmp(ESTIMATOR,"MCD") )) {
      for(i = 0; i < count; i++) {
        for(j = 0;j < 3;j++) {
        /* Add a bit of uniform random noise to the samples */
          sample_vector[3*i + j] = sample_vector[3*i + j] + (rand() - RAND_MAX / 2)*  
                                               voxel_value_interval[j]/RAND_MAX;
        }
      } 
    }           
  }
  delete_volume(volume_tmp);
  *pcount = count;    
  return(sample_vector);
}

double* Collect_values3_subcortical(Volume volume_inT1,Volume volume_inT2,
				    Volume volume_inPD, Volume volume_subcort,
				    char ref_label,long int* pcount) 
{
  int sizes[MAX_DIMENSIONS];
  long int i;
  int j,k,i1,j1,k1;
  long int count = 0;
  long int sampleno;
  BOOLEAN on_the_border;
  double *sample_vector_tmp;
  double *sample_vector;
  double voxel_value_interval[3];

  voxel_value_interval[0] = (get_volume_real_max(volume_inT1) - 
                             get_volume_real_min(volume_inT1))/DATATYPE_SIZE;
  voxel_value_interval[1] = (get_volume_real_max(volume_inT2) - 
                             get_volume_real_min(volume_inT2))/DATATYPE_SIZE;
  voxel_value_interval[2] = (get_volume_real_max(volume_inPD) - 
                             get_volume_real_min(volume_inPD))/DATATYPE_SIZE;

  srand(time(0)); /* Set seed */
  get_volume_sizes( volume_subcort, sizes ); 

  for(i = 1; i < sizes[0] - 1;i++) {
    for(j = 1;j < sizes[1] - 1;j++) {
      for( k = 1; k < sizes[2] - 1;k++) {
	if (get_volume_real_value(volume_subcort,i,j,k,0,0) > 1.5)
	  count++;
      }
    }
  }

  if(count > 0) {
    if(!(strcmp(ESTIMATOR,"MVE") || strcmp(ESTIMATOR,"MCD") ) && (count > MAXSAMPLES)) {
    /* if we are using MVE or MCD (not implemented) estimator so    
     1) we must get rid of some samples;
     2) we must add little bit of noise to the other samples; */  
      sample_vector_tmp = malloc(count * 3 * sizeof (double)); /* reserve room for 3 values (T1,T2,PD) for each sample */
      if(sample_vector_tmp == NULL) return(NULL);
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_subcort,i,j,k,0,0) > 1.5) {
              sample_vector_tmp[sampleno] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              sample_vector_tmp[sampleno + 1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              sample_vector_tmp[sampleno + 2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              sampleno = sampleno + 3;
	    }
          }
        }

      }
      sample_vector = malloc(MAXSAMPLES * 3 * sizeof(double));
      if(sample_vector == NULL) return(NULL);
      for(i =0;i < MAXSAMPLES;i++) {
        sampleno = rint(rand() *( (double) (MAXSAMPLES - 1) / RAND_MAX));
        sample_vector[3 * i] = sample_vector_tmp[3 * sampleno];
        sample_vector[3 * i + 1] = sample_vector_tmp[3 * sampleno + 1];
        sample_vector[3 * i + 2] = sample_vector_tmp[3 * sampleno + 2];
      }
      free(sample_vector_tmp);
      count = MAXSAMPLES;
    }
    else {
      sample_vector = malloc(3 * count * sizeof (double));
      sampleno = 0;
      for(i = 1; i < sizes[0] - 1;i++) {
        for(j = 1;j < sizes[1] - 1;j++) {
          for( k = 1; k < sizes[2] - 1;k++) {
            if(get_volume_real_value(volume_subcort,i,j,k,0,0) > 1.5) {
              sample_vector[sampleno] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              sample_vector[sampleno + 1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              sample_vector[sampleno + 2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              sampleno = sampleno + 3;              
	    }
          }
        }

      }
    }

    if(!(strcmp(ESTIMATOR,"MVE")  || strcmp(ESTIMATOR,"MCD") )) {
      for(i = 0; i < count; i++) {
        for(j = 0;j < 3;j++) {
        /* Add a bit of uniform random noise to the samples */
          sample_vector[3*i + j] = sample_vector[3*i + j] + (rand() - RAND_MAX / 2)*  
                                               voxel_value_interval[j]/RAND_MAX;
        }
      } 
    }
  }

  *pcount = count;    
  return(sample_vector);
}


/* Get parameters based on a tag_file */

int Estimate_params_from_tags3(char* tag_filename,Volume volume_inT1, Volume volume_inT2, Volume volume_inPD,
                                 pMatrix means[PURE_CLASSES + 1], pMatrix vars[PURE_CLASSES + 1]) 
{
   int n_tag_volumes, num_samples,i,c;
   int adj_num_samples[PURE_CLASSES];
   STRING *labels;
   double **tags;
   double wx,wy,wz;
   double v1,v2,v3;
   int sizes[MAX_DIMENSIONS];   
   Vector3D value;
   Matrix3D tmp_matrix;

   get_volume_sizes( volume_inT1, sizes ); 

   if ( input_tag_file(tag_filename, &n_tag_volumes, &num_samples,
		      &tags, NULL, NULL, NULL, NULL, &labels ) != OK ) 
      return(1);

   for(i = 0;i < PURE_CLASSES; i++) {
     adj_num_samples[i] = 0;
     SetZeroVector(means[i + 1]);
     SetZero(vars[i + 1]);
   }

   /* Compute means */
   for( i = 0;i < num_samples;i++) {
     
     /* get the world coordinate from the tag file */
     wx = tags[i][0];
     wy = tags[i][1];
     wz = tags[i][2];
    
    /* convert world into voxel coordinates. At this point T1 - volume is used
       for this. It may become necessary to test all volumes separately. */
     convert_3D_world_to_voxel(volume_inT1, wx, wy, wz, &v1, &v2, &v3);

     if( v1 > - 0.5 && v2 > -0.5 && v3 > -0.5 && v1 < sizes[0] - 0.5 && 
       v2 < sizes[1] - 0.5  && v3 < sizes[2] - 0.5) {
       c = atoi( labels[i]);
       if( c > 0 && c < PURE_CLASSES + 1) {
         value[0] =  get_volume_real_value(volume_inT1,rint(v1),rint(v2),rint(v3),0,0);
         value[1] =  get_volume_real_value(volume_inT2,rint(v1),rint(v2),rint(v3),0,0);
         value[2] =  get_volume_real_value(volume_inPD,rint(v1),rint(v2),rint(v3),0,0);
         adj_num_samples[c - 1] = adj_num_samples[c - 1] + 1;
         AddVectors(means[c],value,means[c]);
       }
       else {
         printf("\n Warning: There is something curious with the labels. \n");
       }
     }
   }
   for(i = 1;i < PURE_CLASSES + 1;i++) {
     ScalarMultiplyVector(means[i],((double) 1)/ adj_num_samples[i - 1],means[i]);
   }
   /* Then same for variances */
    for( i = 0;i < num_samples;i++) {
     
     /* get the world coordinate from the tag file */
     wx = tags[i][0];
     wy = tags[i][1];
     wz = tags[i][2];
    
    /* convert world into voxel coordinates */
     convert_3D_world_to_voxel(volume_inT1, wx, wy, wz, &v1, &v2, &v3);

     if( v1 > - 0.5 && v2 > -0.5 && v3 > -0.5 && v1 < sizes[0] - 0.5 && 
       v2 < sizes[1] - 0.5  && v3 < sizes[2] - 0.5) {
       c = atoi( labels[i]);
       if( c > 0 && c < PURE_CLASSES + 1) {
         value[0] =  get_volume_real_value(volume_inT1,rint(v1),rint(v2),rint(v3),0,0);
         value[1] =  get_volume_real_value(volume_inT2,rint(v1),rint(v2),rint(v3),0,0);
         value[2] =  get_volume_real_value(volume_inPD,rint(v1),rint(v2),rint(v3),0,0);
         SubtractVectors(value,means[c],value);
         VectorProduct(value,value,tmp_matrix);
         AddMatrices(vars[c],tmp_matrix,vars[c]);
       }
       else {
         printf("\n Warning there is something curious with the labels \n");
       }
     }
   }
   for(i = 1;i < PURE_CLASSES + 1;i++) {
     ScalarMultiply(vars[i],((double) 1)/adj_num_samples[i - 1],  vars[i]);
   }
   SetZeroVector(means[BGLABEL]);  /* Set background mean  */
   SetZero(vars[BGLABEL]);         /* Set background variance */
   SetElement(vars[BGLABEL],1,1,10.0);
   SetElement(vars[BGLABEL],2,2,10.0);
   SetElement(vars[BGLABEL],3,3,10.0);
   free_tag_points(n_tag_volumes, num_samples,
		  tags, NULL, NULL, NULL, NULL, labels );
   return(0);
}




/* -------------------------------------------------------------------- */
/* Opens the input images and mask image. Returns 0 if ok.              */


int Open_images3(char* inT1_fn, char* inT2_fn, char* inPD_fn, 
                 char* mask_fn, Volume* pvolume_inT1, 
                 Volume* pvolume_inT2,Volume* pvolume_inPD,
                 Volume* pvolume_mask)
{
  if(input_volume(inT1_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_inT1, (minc_input_options *) NULL) != OK)
    return(1);
  if(input_volume(inT2_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_inT2, (minc_input_options *) NULL) != OK)
    return(2);
  if(input_volume(inPD_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_inPD, (minc_input_options *) NULL) != OK)
    return(3);
  if(input_volume(mask_fn,3,NULL,NC_UNSPECIFIED,FALSE,0.0, 0.0, 
                 TRUE, pvolume_mask, (minc_input_options *) NULL) != OK)
    return(4);



  /* Check that each image really has three dimensions */

  if(get_volume_n_dimensions(*pvolume_inT1) != 3)
    return(5);
  if(get_volume_n_dimensions(*pvolume_inT2) != 3)
    return(6);
  if(get_volume_n_dimensions(*pvolume_inPD) != 3)
    return(7);
  if(get_volume_n_dimensions(*pvolume_mask) != 3)
    return(8);

  /* If we have got this far everything is ok. */

  return(0);
}



/* -------------------------------------------------------------------
  Functions that do tiny tasks for other functions from this file.
  They are not called from main 
  -------------------------------------------------------------------- */ 
/* No explanation needed ... */

/*
char Are_same(char label1, char label2) {return(label1 == label2);} */

/* Figures out whether given tissue types share a component.
   Note that if tissue types are same this function returns zero. 
   Required for Markov Random field computations 
 */

/* char Are_similar(char label1,char label2)
{ 
  switch(label1) {
    case BGLABEL: if (label2 == CSFBGLABEL) return(1);
       else return(0);
       break;
    case WMLABEL: if( label2 == WMGMLABEL) return(1);
       else return(0);
       break;
    case GMLABEL: if( (label2 == WMGMLABEL) || (label2 == GMCSFLABEL)) return(1);
       else return(0);
       break;
    case CSFLABEL: if( (label2 == GMCSFLABEL) || (label2 == CSFBGLABEL)) return(1);
       else return(0);
       break;
    case WMGMLABEL: if( (label2 == WMLABEL) || (label2 == GMLABEL)) return(1);
       else return(0);
       break;
    case GMCSFLABEL: if( (label2 == GMLABEL) || (label2 == CSFLABEL)) return(1);
       else return(0);
       break;
    case CSFBGLABEL: if( (label2 == CSFLABEL) || (label2 == BGLABEL)) return(1);
       else return(0);
       break;
    default : return(0);  break;
  } 
} 

*/


/* Limits x to the range from 0 to 1. */ 

void Limit_0_1(double* x) 
{
  if( *x < 0) *x = 0.0;
  else if( *x > 1) *x = 1.0 ;

}

/* ---------------------------------------------------------------------------
   Functions that are called from main to do the actual computations.
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

double Compute_Gaussian_likelihood3(pVector value, pVector mean , pMatrix var)

{ 
  double d,exponent;
  Matrix3D invvar;
  Vector3D v;

  d = Determinant(var);
  Invert(var,invvar);
  SubtractVectors(value,mean,v);
  exponent = -0.5 * QuadraticForm(invvar,v);
   
  
  return(exp(exponent)/(sqrt(pow(2 * PI,3) * d)));

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

double Compute_marginalized_likelihood3(pVector value, pVector mean1 , pVector mean2, 
                                       pMatrix var1, pMatrix var2, 
                                       pMatrix measurement_var, 
                                       unsigned int nof_intervals)

{ 
  double lh, t, interval_len;
  int i;  
  Vector3D tmean;
  Matrix3D tvar;  
  Vector3D tmpv1,tmpv2;
  Matrix3D tmpm1,tmpm2,tmpm3;

  interval_len = (double) 1 / nof_intervals;
  lh = 0;
  for(i = 0; i < nof_intervals; i++) {
    t = (i + 0.5) * interval_len;

    ScalarMultiplyVector(mean1,t,tmpv1);       /* Find mean vector for */
    ScalarMultiplyVector(mean2,1 - t,tmpv2);   /* current t.           */
    AddVectors(tmpv1,tmpv2,tmean);

    ScalarMultiply(var1,pow(t,2),tmpm1);       /* Find the covariance   */
    ScalarMultiply(var2,pow(1 - t,2),tmpm2);   /* matrix for current t. */
    AddMatrices(tmpm1,tmpm2,tmpm3);            /* tmpm3 = tmpm1 + tmpm2 */
    AddMatrices(tmpm3,measurement_var,tvar);   /* tvar = tmpm3 + measurement_var */
 
    lh = lh + Compute_Gaussian_likelihood3(value, tmean,tvar)*interval_len;
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
     If false (as it usually is) this saves the trouble calculating this 27 times */
 
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

void Parameter_estimation3(Volume volume_inT1,Volume volume_inT2,Volume volume_inPD, 
                           Volume volume_mask, 
                           Volume probabilities[6],pVector mean[4], pMatrix var[4],
                           pMatrix var_measurement)
{
  int i,j,k,sizes[MAX_DIMENSIONS];
  char c;
  double total_probability;
  /* int adapt_var_measurement = 0; */
  Vector3D tmpvec, intensity;
  Matrix3D tmpmatrix;
  /* BOOLEAN posdef; */
  Matrix3D oldvar[3];

  printf("WM mean");
  PrintVector(mean[WMLABEL]);
  printf("GM mean");
  PrintVector(mean[GMLABEL]);
  printf("CSF mean");
  PrintVector(mean[CSFLABEL]);
  printf("WM var");
  PrintMatrix(var[WMLABEL]);
  printf("GM var");
  PrintMatrix(var[GMLABEL]);
  printf("CSF var");
  PrintMatrix(var[CSFLABEL]);

  get_volume_sizes(volume_inT1,sizes);
  
  CopyMatrix(var[1],oldvar[0]);
  CopyMatrix(var[2],oldvar[1]);
  CopyMatrix(var[3],oldvar[2]);

  /* First re-estimate means */
 
  for(c = 1;c < 4;c++) {
    total_probability = 0;
    SetZeroVector(mean[c]);
    for( i = 0; i < sizes[0];i++) {
      for( j = 0;j < sizes[1];j++) {
        for( k = 0; k < sizes[2]; k++) {
          if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
            total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
            intensity[0] = get_volume_real_value(volume_inT1,i,j,k,0,0);
            intensity[1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
            intensity[2] = get_volume_real_value(volume_inPD,i,j,k,0,0);

            ScalarMultiplyVector(intensity,get_volume_real_value(probabilities[c - 1],i,j,k,0,0),tmpvec);
            AddVectors(mean[c],tmpvec,mean[c]);
            
          }
	}
      }
    }
    ScalarMultiplyVector(mean[c], 1 / total_probability, mean[c]);
  }
  
  /* Then re-estimate the variances.
     Case 1: If the model contains physiological noise components. */
  if( (!IsZero(var[1])) || (!IsZero(var[2])) || /* Noise model contains tissue-wise */
     (!IsZero(var[3]))) {                        /* components.                      */
    for(c = 1;c < 4;c++) {
      total_probability = 0;
      SetZero(var[c]);
      for( i = 0; i < sizes[0];i++) {
        for( j = 0;j < sizes[1];j++) {
          for( k = 0; k < sizes[2]; k++) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
              total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
              intensity[0] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              intensity[1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              intensity[2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              SubtractVectors(intensity,mean[c],intensity);
              VectorProduct(intensity,intensity,tmpmatrix);
              ScalarMultiply(tmpmatrix,get_volume_real_value(probabilities[c - 1],i,j,k,0,0),
                             tmpmatrix);
              AddMatrices(var[c],tmpmatrix,var[c]);
            }
          }
        }
      }
      ScalarMultiply(var[c],1/total_probability,var[c]);
    }
    if(!IsZero(var_measurement)) {
      /* I don't know how to this at the moment, so adapting measurement 
         variance is not an option */
      /* if(adapt_var_measurement) 
       *var_measurement = MIN3(var[1],var[2],var[3])*MEASUREMENT_FRACTION; */
      for(c = 1;c < 4;c++) {
        SubtractMatrices(var[c],var_measurement,var[c]);
        if(!IsPositiveDefinite(var[c])) {
          CopyMatrix(oldvar[c - 1],var[c]);
        }          
      }
    }
  }
  else {    /* The model contains only measurement noise. We update only
               measurement variance.          */
    SetZero(var_measurement);
    total_probability = 0;
    for(c = 1;c < 4;c++) {
      for( i = 0; i < sizes[0];i++) {
        for( j = 0;j < sizes[1];j++) {
          for( k = 0; k < sizes[2]; k++) {
            if(get_volume_real_value(volume_mask,i,j,k,0,0) > MASK_TR) {
              total_probability = total_probability + 
                                get_volume_real_value(probabilities[c - 1],i,j,k,0,0);
              intensity[0] = get_volume_real_value(volume_inT1,i,j,k,0,0);
              intensity[1] = get_volume_real_value(volume_inT2,i,j,k,0,0);
              intensity[2] = get_volume_real_value(volume_inPD,i,j,k,0,0);
              SubtractVectors(intensity,mean[c],intensity);
              VectorProduct(intensity,intensity,tmpmatrix);
              ScalarMultiply(tmpmatrix,get_volume_real_value(probabilities[c - 1],i,j,k,0,0),
                             tmpmatrix);
              AddMatrices(var_measurement,tmpmatrix,var_measurement); 
            }
          }
        }
      }
    }
    ScalarMultiply(var_measurement, 1 / total_probability,var_measurement);
 
  }
  printf("WM mean");
  PrintVector(mean[WMLABEL]);
  printf("GM mean");
  PrintVector(mean[GMLABEL]);
  printf("CSF mean");
  PrintVector(mean[CSFLABEL]);
  printf("WM var");
  PrintMatrix(var[WMLABEL]);
  printf("GM var");
  PrintMatrix(var[GMLABEL]);
  printf("CSF var");
  PrintMatrix(var[CSFLABEL]);
  
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

int Compute_partial_volume_vectors3(Volume volume_inT1, Volume volume_inT2, 
                                    Volume volume_inPD, Volume volume_classified,
                                    Volume volume_pve[PURE_CLASSES],
                                    pVector means[PURE_CLASSES + 1], 
                                    pMatrix vars[PURE_CLASSES + 1], pMatrix var_measurement )
{
  int sizes[MAX_DIMENSIONS];
  int i,j,k;
  char c;
  double t;  /* mixing proportion to be estimated */
  Vector3D value;
  Vector3D diff1,diff2;
  Matrix3D var_tmp;
  Matrix3D inv_var_wmgm,inv_var_gmcsf,inv_var_csfbg,inv_var_scgm,inv_var_wmsc;

  get_volume_sizes(volume_inT1,sizes);
   
  /* Approximate partial volume covariances by choosing t = 0.5 */
  /* Covariances are treated as constants although they are not */

 
  AddMatrices(vars[WMLABEL],vars[GMLABEL],var_tmp);
  ScalarMultiply(var_tmp,0.25,var_tmp);
  AddMatrices(var_tmp,var_measurement,var_tmp);
  Invert(var_tmp,inv_var_wmgm);

  AddMatrices(vars[GMLABEL],vars[CSFLABEL],var_tmp);
  ScalarMultiply(var_tmp,0.25,var_tmp);
  AddMatrices(var_tmp,var_measurement,var_tmp);
  Invert(var_tmp,inv_var_gmcsf);

  AddMatrices(vars[BGLABEL],vars[CSFLABEL],var_tmp);
  ScalarMultiply(var_tmp,0.25,var_tmp);
  AddMatrices(var_tmp,var_measurement,var_tmp);
  Invert(var_tmp,inv_var_csfbg);

  AddMatrices(vars[SCLABEL],vars[GMLABEL],var_tmp);
  ScalarMultiply(var_tmp,0.25,var_tmp);
  AddMatrices(var_tmp,var_measurement,var_tmp);
  Invert(var_tmp,inv_var_scgm);

  AddMatrices(vars[SCLABEL],vars[GMLABEL],var_tmp);
  ScalarMultiply(var_tmp,0.25,var_tmp);
  AddMatrices(var_tmp,var_measurement,var_tmp);
  Invert(var_tmp,inv_var_scgm);

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
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             SubtractVectors(value,means[GMLABEL],diff1);
             SubtractVectors(means[WMLABEL],means[GMLABEL],diff2);
             t = QuadraticForm2(inv_var_gmcsf,diff1,diff2)/QuadraticForm(inv_var_wmgm,diff2);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case GMCSFLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             SubtractVectors(value,means[CSFLABEL],diff1);
             SubtractVectors(means[GMLABEL],means[CSFLABEL],diff2);
             t = QuadraticForm2(inv_var_gmcsf,diff1,diff2)/QuadraticForm(inv_var_gmcsf,diff2);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0, t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case CSFBGLABEL:
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             SubtractVectors(value,means[BGLABEL],diff1);
             SubtractVectors(means[CSFLABEL],means[BGLABEL],diff2);
             t = QuadraticForm2(inv_var_csfbg,diff1,diff2)/QuadraticForm(inv_var_csfbg,diff2);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,t); 
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break;
        case WMSCLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             SubtractVectors(value,means[SCLABEL],diff1);
             SubtractVectors(means[WMLABEL],means[SCLABEL],diff2);
             t = QuadraticForm2(inv_var_gmcsf,diff1,diff2)/QuadraticForm(inv_var_gmcsf,diff2);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case SCGMLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             SubtractVectors(value,means[GMLABEL],diff1);
             SubtractVectors(means[SCLABEL],means[GMLABEL],diff2);
             t = QuadraticForm2(inv_var_gmcsf,diff1,diff2)/QuadraticForm(inv_var_wmgm,diff2);
             Limit_0_1(&t);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             break; 
        default: return(1); break;
        }
      }
    }
  }
  return(0);
}


/* And finally, a function that calculates the partial volume vectors based on exact 
   Maximum likelihood criterion. Uses grid search, but should not be too time taking. 
   However, is more computationally intensive than function above.

 */

int Compute_partial_volume_vectors3_ml(Volume volume_inT1, Volume volume_inT2, 
                                    Volume volume_inPD, Volume volume_classified,
                                    Volume volume_pve[PURE_CLASSES],
                                    pVector means[PURE_CLASSES + 1], 
                                    pMatrix vars[PURE_CLASSES + 1], pMatrix var_measurement )
{
  int sizes[MAX_DIMENSIONS];
  int i,j,k,n;
  char c;
  double t;  /* mixing proportion to be estimated */
  double prob;
  Vector3D value;
  Vector3D mean_tmp1,mean_tmp2;
  double maxvalue = - 1.0;
  int maxindex = 0;

  Matrix3D var_tmp1,var_tmp2;

  get_volume_sizes(volume_inT1,sizes);
   
 
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
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               ScalarMultiplyVector(means[WMLABEL],((double) n)/100,mean_tmp1);
               ScalarMultiplyVector(means[GMLABEL],1 - ((double) n)/100,mean_tmp2);
               AddVectors(mean_tmp1,mean_tmp2,mean_tmp1);

               ScalarMultiply(vars[WMLABEL],pow(((double) n)/100,2),var_tmp1);
               ScalarMultiply(vars[GMLABEL],pow(1 - ((double) n)/100,2),var_tmp2);
               AddMatrices(var_tmp1,var_tmp2,var_tmp1);
               AddMatrices(var_tmp1,var_measurement,var_tmp1);
               prob = Compute_Gaussian_likelihood3(value,mean_tmp1,var_tmp1);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;  
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);  
	     set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
             break; 
        case GMCSFLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             maxvalue = -1.0;
             maxindex = 0;             
             for(n = 0;n < 101;n++) {
               ScalarMultiplyVector(means[GMLABEL],((double) n)/100,mean_tmp1);
               ScalarMultiplyVector(means[CSFLABEL],1 - ((double) n)/100,mean_tmp2);
               AddVectors(mean_tmp1,mean_tmp2,mean_tmp1);

               ScalarMultiply(vars[GMLABEL],pow(((double) n)/100,2),var_tmp1);
               ScalarMultiply(vars[CSFLABEL],pow(1 - ((double) n)/100,2),var_tmp2);
               AddMatrices(var_tmp1,var_tmp2,var_tmp1);
               AddMatrices(var_tmp1,var_measurement,var_tmp1);
               prob = Compute_Gaussian_likelihood3(value,mean_tmp1,var_tmp1);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100; 
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0, t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,1 - t);      
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);         
             break; 
        case CSFBGLABEL:
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             maxvalue = -1.0;
             maxindex = 0;             
             for(n = 0;n < 101;n++) {
               ScalarMultiplyVector(means[CSFLABEL],((double) n)/100,mean_tmp1);
               ScalarMultiplyVector(means[BGLABEL],1 - ((double) n)/100,mean_tmp2);
               AddVectors(mean_tmp1,mean_tmp2,mean_tmp1);

               ScalarMultiply(vars[CSFLABEL],pow(((double) n)/100,2),var_tmp1);
               ScalarMultiply(vars[BGLABEL],pow(1 - ((double) n)/100,2),var_tmp2);
               AddMatrices(var_tmp1,var_tmp2,var_tmp1);
               AddMatrices(var_tmp1,var_measurement,var_tmp1);
               prob = Compute_Gaussian_likelihood3(value,mean_tmp1,var_tmp1);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,t); 
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,0.0);
	     break;
        case WMSCLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               ScalarMultiplyVector(means[WMLABEL],((double) n)/100,mean_tmp1);
               ScalarMultiplyVector(means[SCLABEL],1 - ((double) n)/100,mean_tmp2);
               AddVectors(mean_tmp1,mean_tmp2,mean_tmp1);

               ScalarMultiply(vars[WMLABEL],pow(((double) n)/100,2),var_tmp1);
               ScalarMultiply(vars[SCLABEL],pow(1 - ((double) n)/100,2),var_tmp2);
               AddMatrices(var_tmp1,var_tmp2,var_tmp1);
               AddMatrices(var_tmp1,var_measurement,var_tmp1);
               prob = Compute_Gaussian_likelihood3(value,mean_tmp1,var_tmp1);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;  
             set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);  
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,0.0);              
             break;
        case SCGMLABEL: 
             InitializeVector(value,get_volume_real_value(volume_inT1,i,j,k,0,0),
                              get_volume_real_value(volume_inT2,i,j,k,0,0),
                              get_volume_real_value(volume_inPD,i,j,k,0,0));
             maxvalue = -1.0;
             maxindex = 0;
             for(n = 0;n < 101;n++) {
               ScalarMultiplyVector(means[SCLABEL],((double) n)/100,mean_tmp1);
               ScalarMultiplyVector(means[GMLABEL],1 - ((double) n)/100,mean_tmp2);
               AddVectors(mean_tmp1,mean_tmp2,mean_tmp1);

               ScalarMultiply(vars[SCLABEL],pow(((double) n)/100,2),var_tmp1);
               ScalarMultiply(vars[GMLABEL],pow(1 - ((double) n)/100,2),var_tmp2);
               AddMatrices(var_tmp1,var_tmp2,var_tmp1);
               AddMatrices(var_tmp1,var_measurement,var_tmp1);
               prob = Compute_Gaussian_likelihood3(value,mean_tmp1,var_tmp1);
               if(prob > maxvalue) {
                 maxvalue = prob;
                 maxindex = n;
               }
	     }
             
             t = ((double) maxindex)/100;  
             set_volume_real_value(volume_pve[SCLABEL - 1],i,j,k,0,0,t);
             set_volume_real_value(volume_pve[GMLABEL - 1],i,j,k,0,0,1 - t);
             set_volume_real_value(volume_pve[CSFLABEL - 1],i,j,k,0,0,0.0);  
	     set_volume_real_value(volume_pve[WMLABEL - 1],i,j,k,0,0,0.0);   
             break;
        default: return(1); break;
        }
      }
    }
  }
  return(0);
}
