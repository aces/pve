#include "matrix3.h"
#include <stdio.h>
#include <math.h>

/* Some basic routines for 3x3 matrices for dealing 
   with  multinormal distributions */

/* double GetElement(pMatrix pm,char i, char j)
{
  return(pm[(i - 1) * 3 + j - 1]);

}
*/
void SetElement(pMatrix pm,char i, char j, double value) 
{
  pm[(i -1) * 3 + j - 1] = value;

}
void PrintMatrix(pMatrix pm)
{
  char i,j;
  for(i = 1;i < 4;i++) {
    printf("\n");
    for(j = 1; j < 4;j++) {
      printf("%f ", GetElement(pm,i,j));
    }
  }
  printf("\n");
}

void PrintVector(pVector pv)
{
  printf("\n %f %f %f \n",pv[0],pv[1],pv[2]);
}


void SetZero(pMatrix pm)
{
  char i,j;
  for(i = 1;i < 4;i++) {
    for(j = 1; j < 4;j++) {
      SetElement(pm,i,j,0);
    }
  }
}

void CopyMatrix(pMatrix sourcem,pMatrix targetm)
{
  char i,j;
  for(i = 1;i < 4;i++) {
    for(j = 1; j < 4;j++) {
      SetElement(targetm,i,j,GetElement(sourcem,i,j));
    }
  }
}

void CopyVector(pVector sourcev, pVector targetv)
{
  targetv[0] = sourcev[0];
  targetv[1] = sourcev[1];
  targetv[2] = sourcev[2];
}

void SetZeroVector(pVector pv)
{
  pv[0] = 0; pv[1] = 0 ; pv[2] = 0;
}


void InitializeVector(pVector pv,double x1,double x2, double x3)
{
  pv[0] = x1;
  pv[1] = x2;
  pv[2] = x3;

}

void ScalarMultiplyVector(pVector pv1,double scalar, pVector result)
{
 
  result[0] = pv1[0] * scalar;
  result[1] = pv1[1] * scalar;
  result[2] = pv1[2] * scalar;

}

void ScalarMultiply(pMatrix pm1,double scalar, pMatrix result)
{
  char i1,i2;

  for(i1 = 1;i1 < 4;i1++) {
    for(i2 = 1;i2 < 4;i2++) {
      SetElement(result,i1,i2,scalar * GetElement(pm1,i1,i2));
    }
  }
}

void AddVectors(pVector v1, pVector v2, pVector result) 
{
  result[0] = v1[0] + v2[0];
  result[1] = v1[1] + v2[1];
  result[2] = v1[2] + v2[2];
 
}

void SubtractVectors(pVector v1, pVector v2, pVector result) 
{
  result[0] = v1[0] - v2[0];
  result[1] = v1[1] - v2[1];
  result[2] = v1[2] - v2[2];

}

void AddMatrices(pMatrix pm1, pMatrix pm2, pMatrix result)
{
  char i1,i2;

  for(i1 = 1;i1 < 4;i1++) {
    for(i2 = 1;i2 < 4;i2++) {
      SetElement(result,i1,i2,GetElement(pm1,i1,i2) + GetElement(pm2,i1,i2));
    }
  }

}

void SubtractMatrices(pMatrix pm1, pMatrix pm2, pMatrix result)
{
  char i1,i2;

  for(i1 = 1;i1 < 4;i1++) {
    for(i2 = 1;i2 < 4;i2++) {
      SetElement(result,i1,i2,GetElement(pm1,i1,i2) - GetElement(pm2,i1,i2));
    }
  }

}

double ScalarProduct(pVector pv1, pVector pv2)
{
  return(pv1[0] * pv2[0] + pv1[1] * pv2[1] + pv1[2] * pv2[2]);
}

void VectorProduct(pVector pv1, pVector pv2, pMatrix result)
{
  char i,j;
  for(i = 1; i < 4;i++) {
    for(j = 1; j < 4;j++) {
      SetElement(result,i,j,pv1[i - 1] * pv2[j - 1]);
    }
  }
}


double QuadraticForm(pMatrix pm, pVector v)
{
  double d = 0;
  char i1,i2;

  for(i1 = 1;i1 < 4;i1++) {
    for(i2 = 1;i2 < 4;i2++) {
      d = d + v[i1 - 1] * GetElement(pm,i1,i2) * v[i2 - 1];
    }
  }
  return(d);
}

double QuadraticForm2(pMatrix pm, pVector v1, pVector v2)
{
  double d = 0;
  char i1,i2;

  for(i1 = 1;i1 < 4;i1++) {
    for(i2 = 1;i2 < 4;i2++) {
      d = d + v1[i1 - 1] * GetElement(pm,i1,i2) * v2[i2 - 1];
    }
  }
  return(d);
}

int IsZero(pMatrix pm)
{
  int return_value  = 0;
  char i,j;

  for(i = 1;i < 4;i++) {
    for( j = 1;j < 4;j++) {
      if(fabs(GetElement(pm,i,j)) < VERY_SMALL) return_value = 1;
    }
  }
  return(return_value);
}


int IsSymmetric(pMatrix pm) 
{
 
  if(( fabs(GetElement(pm,1,2) - GetElement(pm,2,1)) < VERY_SMALL) &&
     ( fabs(GetElement(pm,1,3) - GetElement(pm,3,1)) < VERY_SMALL) &&
     ( fabs(GetElement(pm,2,3) - GetElement(pm,3,2)) < VERY_SMALL))
    return(1);
  else return(0);

}

int IsPositiveDefinite(pMatrix pm) 
{
  if(!(GetElement(pm,1,1) > 0)) return(0);
  else if(!(Determinant2x2(GetElement(pm,1,1), GetElement(pm,1,2),
                           GetElement(pm,2,1), GetElement(pm,2,2)) > 0)) return(0);
  else if(!(Determinant(pm) > 0)) return(0);
  else return(1);
}


double Determinant2x2(double m11, double m12, double m21, double m22)
{
  return(m11 * m22 - m12 * m21); 

} 


double Determinant(pMatrix pm)
{
  double d;

  d = GetElement(pm,1,1) * Determinant2x2(GetElement(pm,2,2),GetElement(pm,2,3),
                                          GetElement(pm,3,2),GetElement(pm,3,3));
  d = d + GetElement(pm,1,2) * Determinant2x2(GetElement(pm,2,3),GetElement(pm,2,1),
                                          GetElement(pm,3,3),GetElement(pm,3,1));
  d = d + GetElement(pm,1,3) * Determinant2x2(GetElement(pm,2,1),GetElement(pm,2,2),
                                          GetElement(pm,3,1),GetElement(pm,3,2));
  
  return(d);

}

int Invert(pMatrix pm, pMatrix pinvm) 
{
  double d;
  double tmp;

  d = Determinant(pm);
  if( fabs(d) <  VERY_SMALL) {
    return(1);
  }
  else {
    tmp = Determinant2x2(GetElement(pm,2,2),GetElement(pm,2,3),
                         GetElement(pm,3,2),GetElement(pm,3,3));
    SetElement(pinvm,1,1,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,3),GetElement(pm,1,2),
                         GetElement(pm,3,3),GetElement(pm,3,2));
    SetElement(pinvm,1,2,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,2),GetElement(pm,1,3),
                         GetElement(pm,2,2),GetElement(pm,2,3));
    SetElement(pinvm,1,3,tmp / d);
    tmp = Determinant2x2(GetElement(pm,2,3),GetElement(pm,2,1),
                         GetElement(pm,3,3),GetElement(pm,3,1));
    SetElement(pinvm,2,1,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,1),GetElement(pm,1,3),
                         GetElement(pm,3,1),GetElement(pm,3,3));
    SetElement(pinvm,2,2,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,3),GetElement(pm,1,1),
                         GetElement(pm,2,3),GetElement(pm,2,1));
    SetElement(pinvm,2,3,tmp / d);
    tmp = Determinant2x2(GetElement(pm,2,1),GetElement(pm,2,2),
                         GetElement(pm,3,1),GetElement(pm,3,2));
    SetElement(pinvm,3,1,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,2),GetElement(pm,1,1),
                         GetElement(pm,3,2),GetElement(pm,3,1));
    SetElement(pinvm,3,2,tmp / d);
    tmp = Determinant2x2(GetElement(pm,1,1),GetElement(pm,1,2),
                         GetElement(pm,2,1),GetElement(pm,2,2));
    SetElement(pinvm,3,3,tmp / d);
    return(0);
  }
}






