/* Some basic routines for 3x3 matrices for dealing 
   with  multinormal distributions */
#include <math.h>

#define VERY_SMALL 1E-16

typedef double  Matrix[9];
typedef double  Vector[3] ;
typedef double* pMatrix;
typedef double* pVector;

#define GetElement(pm, i, j) pm[(i - 1) * 3 + j - 1]
/* #define SetElement(pm, i, j, value) pm[(i -1) * 3 + j - 1] = value */

void SetElement(pMatrix pm,char i, char j, double value);
void PrintMatrix(pMatrix pm);
void PrintVector(pVector pv);
void SetZero(pMatrix pm);
void CopyMatrix(pMatrix sourcem,pMatrix targetm);
void CopyVector(pVector sourcev, pVector targetv);
void InitializeVector(pVector pv,double x1,double x2, double x3);
void SetZeroVector(pVector pv);

void ScalarMultiplyVector(pVector pv1,double scalar, pVector result);
void ScalarMultiply(pMatrix pm1,double scalar, pMatrix result);
void AddVectors(pVector v1, pVector v2, pVector result);
void SubtractVectors(pVector v1, pVector v2, pVector result);
void AddMatrices(pMatrix pm1, pMatrix pm2, pMatrix result);
void SubtractMatrices(pMatrix pm1, pMatrix pm2, pMatrix result);
double ScalarProduct(pVector pv1, pVector pv2);
void VectorProduct(pVector pv1, pVector pv2, pMatrix result);
double QuadraticForm(pMatrix pm, pVector v);
double QuadraticForm2(pMatrix pm, pVector v1, pVector v2);

int IsZero(pMatrix pm);
int IsSymmetric(pMatrix pm);
int IsPositiveDefinite(pMatrix pm);

double Determinant2x2(double m11, double m12, double m21, double m22);
double Determinant(pMatrix pm);
int Invert(pMatrix pm, pMatrix pinvm);


