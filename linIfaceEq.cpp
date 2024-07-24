#include <math.h>
#include "nrutil.h"
#include "nr.h"
#include <iostream>
#include <string.h>
#include <stdlib.h>
#define TINY 1.0e-20     /* a small number */
#define PI 3.1415920     /* err... what's this??? */




double getGammaD(void);
double getEpsilon(void);
double bessk0(double inArgument);
double bessi0(double inArgument);
double getLatticeSpacing(void);
int getGridLength(void);

double** fillKoeffizientenMatrix(int inNumberGridPoints, double inLatticeSpacing, double inGammaD, double inEpsilon);

double* initialisiereInhomogenitaet(int inNumberGridPoints, double inLatticeSpacing);
void solveLes(double** koeffizientenMatrix, int inNumberGridPoints, int* indx ,double* inhomogenitaet);
void decomposeMatrixToLAndUTriangularMatrices(double** koeffizientenMatrix, int inNumberGridPoints, int* indx, double* d);

int main(int argc, char* argv)
{
  double gammaD = PI / 2.0; 
 
  double** koeffizientenMatrix;

  double* inhomogenitaet;
  double d;
  double latticeSpacing = 1.0;
  int gridLength = 0;
  int* indx;  
  double epsilon = 10.0;
  

  gridLength = getGridLength();
  latticeSpacing = getLatticeSpacing();
  epsilon = getEpsilon();
  gammaD = getGammaD();
  int numberGridPoints = gridLength/latticeSpacing;

  inhomogenitaet = initialisiereInhomogenitaet(numberGridPoints, latticeSpacing);
  koeffizientenMatrix = fillKoeffizientenMatrix(numberGridPoints, latticeSpacing, gammaD, epsilon);
 
 
  indx = ivector(1,numberGridPoints);
  decomposeMatrixToLAndUTriangularMatrices(koeffizientenMatrix, numberGridPoints, indx, &d);
  solveLes(koeffizientenMatrix, numberGridPoints, indx, inhomogenitaet);


 
  FILE* fp;
  fp = fopen("loesung.dat", "w");     
  for (int i=1 ; i <= numberGridPoints; i++)
    {
      
      fprintf(fp, "%d  %g \n", i, inhomogenitaet[i]);
    }
  fclose(fp);

   
      
      
  return 0;
}

int getGridLength()
{
  int outGridLength = 0;
  
  std::cout << "please enter the gridlength:" << std::endl; 
 
  std::cin >> outGridLength ;
  
  return outGridLength;
} 

double getLatticeSpacing()
{
  double outLatticeSpacing = 0.0;
  
  std::cout << "please enter the latticespacing:" << std::endl;
 
  std::cin >> outLatticeSpacing;
 
  return outLatticeSpacing;
}

double getGammaD()
{
  double outGammaD = 1.0;

  std::cout << "please enter gammaD :" << std::endl;

  std::cin >> outGammaD;

  return outGammaD;
}

double getEpsilon()
{
  double outEpsilon = 10.0;

  std::cout << "please enter epsilon :" << std::endl;

  std::cin >> outEpsilon;

  return outEpsilon;
}


 double* initialisiereInhomogenitaet(int inNumberGridPoints, double inLatticeSpacing)
  {
  
    double* inhomogenitaet;
  
    inhomogenitaet = dvector(1,inNumberGridPoints);
  
    for (int i=1; i<=inNumberGridPoints; i++)
      {
	inhomogenitaet[i]=1.0;
      }

    FILE* fp;
    fp = fopen("inh.dat", "w");     
    for (int i=0; i <= inNumberGridPoints-1; i++)
      {
	fprintf(fp, "%g  %g\n", i*inLatticeSpacing, inhomogenitaet[i]);
      }
    fclose(fp);
   
    return inhomogenitaet;
  }


 double** fillKoeffizientenMatrix(int inNumberGridPoints, double inLatticeSpacing, double inGammaD, double inEpsilon)
  {
    double** koeffizientenMatrix;

    koeffizientenMatrix = dmatrix(1, inNumberGridPoints, 1, inNumberGridPoints);

      {
	for(int i=1; i <= inNumberGridPoints; i++)
	{
	  for(int j=1; j <= inNumberGridPoints; j++)
	  {
	    if (j==i)                                 //main diagonal elements      
	      {                                                                     
		koeffizientenMatrix[i][j] = (inGammaD / inLatticeSpacing) - inLatticeSpacing *2/PI * sqrt(inEpsilon/PI) * (log(inEpsilon * inLatticeSpacing / 4) - 1); 
		if (i == inNumberGridPoints)
		  {
		    koeffizientenMatrix[i][j] = inGammaD/inLatticeSpacing - inLatticeSpacing * 2/PI * sqrt(inEpsilon/PI) * (log(inEpsilon * inLatticeSpacing/4) - 1);
		  }
	      }
	    else if (j <=(i-1))                       //generic elements below the main diagonal
	    {
	    koeffizientenMatrix[i][j] = (2/PI)*(sqrt(inEpsilon/PI))*inLatticeSpacing*exp((i-j)*inEpsilon*inLatticeSpacing)*bessk0((i-j)*inEpsilon*inLatticeSpacing) ;
	    }  
	    else if (j == (i+1))                      //elements due to derivative on first upper beside diagonal
	    {
	    koeffizientenMatrix[i][j]= -(inGammaD/inLatticeSpacing)+(2/PI)*sqrt(inEpsilon/PI)*inLatticeSpacing*exp(-inEpsilon*inLatticeSpacing)*bessk0(inEpsilon*inLatticeSpacing);
	    }
            else if (j > (i+1))                       //elements on second or later upper beside diagonals
	    {
	    koeffizientenMatrix[i][j] = inLatticeSpacing*(2/PI)*(sqrt(inEpsilon/PI))*exp(-(j-i)*inEpsilon*inLatticeSpacing)*bessk0((j-i)*inEpsilon*inLatticeSpacing);
	    }
	  }
	}
  
       
      }
      
      
     
 return koeffizientenMatrix;
}

//double** fillOldKoeffizientenMatrix(int inNumberGridPoints, double inLatticeSpacing, double gammaD)
//{
//  double** oldKoeffizientenMatrix;

//  oldKoeffizientenMatrix = dmatrix(1, inNumberGridPoints, 1, inNumberGridPoints);

//    {
//	for(int i=1; i <= inNumberGridPoints; i++)
//{
//  for(int j=1; j <= inNumberGridPoints; j++)
//  {
//    if (j==i)                                 //main diagonal elements      //control matrix output shows that
//    {                                                                     //the i counts rows,j columns, so
//    oldKoeffizientenMatrix[i][j]                                           //allocating elements below the main
//      = (gammaD / inLatticeSpacing) + 2 * sqrt(2*inLatticeSpacing)/PI;  //diagonal to case j<=(i-1) is correct
//    if (i == inNumberGridPoints)
//      {
//        oldKoeffizientenMatrix[i][j] = 2 * sqrt(2 * inLatticeSpacing)/PI;
//      }
//    }
//    else if (j <=(i-1))                       //generic elements below the main diagonal
//    {
//    oldKoeffizientenMatrix[i][j] = sqrt(2 * inLatticeSpacing) / (PI  * sqrt(i-j));
//    }  
//    else if (j == (i+1))                      //elements due to derivative on first upper beside diagonal
//    {
//    oldKoeffizientenMatrix[i][j]= -gammaD/inLatticeSpacing;
///    }
//         else if (j > (i+1))                       //elements on second or later upper beside diagonals
//    {
//    oldKoeffizientenMatrix[i][j] = 0;
//    }
//  }
//}
  
	//std::cout << "after oldKoeffizientenMatrix" << std::endl; 
//}   
// return oldKoeffizientenMatrix;
// }

//double** getDifferenceOfMatrices(double** inOldKoeffizientenMatrix,double** inKoeffizientenMatrix,int inNumberGridPoints)
//{
// double** outDifferenceMatrix;
//outDifferenceMatrix = dmatrix(1,inNumberGridPoints,1,inNumberGridPoints);
  
//int i,j;
//for(i=1; i <= inNumberGridPoints; i++)
//  {
//    for(j=1; j <= inNumberGridPoints; j++)
//{
//  //std::cout << "in getDifferenceOfMatrices" << std::endl;
//  outDifferenceMatrix[i][j] = inKoeffizientenMatrix[i][j] - inOldKoeffizientenMatrix[i][j];
//  
//      }
      
      //std::cout << "after getDifferenceOfMatrices" << std::endl;
//  }
//return outDifferenceMatrix;
//}

 


void decomposeMatrixToLAndUTriangularMatrices(double** koeffizientenMatrix, int inNumberGridPoints, int* indx,double* d)
{
  int i,j,imax,k;
  double big,dum,sum,temp;
  double *vv;

  vv = dvector(1,inNumberGridPoints);
  //std::cout << "richtich"  << endl;
  *d=1.0;
      
  for (i = 1; i <= inNumberGridPoints; i++)
    {
	  
      big = 0.0;
      for (j = 1; j <= inNumberGridPoints; j++)
	
	if ((temp=fabs(koeffizientenMatrix[i][j])) > big ) big=temp;
      if (big == 0.0) nrerror("Singular matrix in routine ludcomp");
      vv[i]=1.0/big;
    } 
    
  
  for (j = 1; j <= inNumberGridPoints; j++)
    {
      for (i = 1; i < j; i++ )
	{ 
	  sum = koeffizientenMatrix[i][j];
	  for (k = 1; k < i; k++)
	    {
	      sum -= koeffizientenMatrix[i][k]*koeffizientenMatrix[k][j];
	    }
	  koeffizientenMatrix[i][j]=sum;
	}
      

      big=0.0;
      for (i = j; i <= inNumberGridPoints; i++)
	{
	  sum=koeffizientenMatrix[i][j];
	  for (k=1; k<j; k++)
	    {
	      sum -= koeffizientenMatrix[i][k]*koeffizientenMatrix[k][j];
	    }

	  koeffizientenMatrix[i][j]=sum;
	  
	  if ( (dum=vv[i]*fabs(sum)) >= big) 
	    {
	      big=dum;
	      imax=i;
	    }
	}

	  
      if (j != imax) 
	{
	  //cout << "Vor for" << endl;
	  for (k=1; k<= inNumberGridPoints; k++)
	    {
	      //cout << "Schleifendurchlauf " << k << "angefangen." << endl;
	      dum = koeffizientenMatrix[imax][k];
	      koeffizientenMatrix[imax][k] = koeffizientenMatrix[j][k];
	      koeffizientenMatrix[j][k] = dum;
	      //cout << "Schleifendurchlauf " << k << "beendet." << endl;
	    }
	  
	  *d = -(*d);
	  vv[imax]=vv[j];
	}
      
      //cout << "imax = " << imax << " , j = " << j << endl;
      
      //cout << "indx[j] = " << indx[j] << endl;
      
      
      indx[j] = imax;
      //cout << "indx[j] = " << indx[j] << endl;
      
      if (koeffizientenMatrix[j][j] == 0.0) 
	{
	  koeffizientenMatrix[j][j]=TINY;
	}
	  
	  
      if (j != inNumberGridPoints)
	{
	  dum=1.0/(koeffizientenMatrix[j][j]);
	  
	  for (i = j+1; i <= inNumberGridPoints; i++) 
	    {
	      koeffizientenMatrix[i][j] *= dum;
	    }
	}
    }
  //cout << "na?" << endl;
  free_dvector(vv,1, inNumberGridPoints);
}
   
 

void solveLes(double** koeffizientenMatrix, int numberGridPoints, int* indx ,double* inhomogenitaet)
  
{
  int i,ii=0,ip,j;
  double sum;
  
  
  for (i=1;i<=numberGridPoints;i++)
    {
      ip=indx[i];
      sum=inhomogenitaet[ip];
      inhomogenitaet[ip]=inhomogenitaet[i];
      if (ii)
	for (j = ii; j <= i-1; j++) sum -= koeffizientenMatrix[i][j]*inhomogenitaet[j];
      else if (sum) ii=i;
      inhomogenitaet[i]=sum;
    }
  
  for (i = numberGridPoints; i >= 1; i--)
    {
      sum=inhomogenitaet[i];
      for (j = i+1; j <= numberGridPoints; j++) sum -= koeffizientenMatrix[i][j]*inhomogenitaet[j];
      inhomogenitaet[i]=sum/koeffizientenMatrix[i][i];
    }
}

double bessk0(double inArgument) //aus nr-code auskopiert
{
  double bessi0(double inArgument);
  double y,ans;
  
  if (inArgument <= 2.0) {
    y=inArgument*inArgument/4.0;
    ans=(-log(inArgument/2.0)*bessi0(inArgument))+(-0.57721566+y*(0.42278420
			+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
			+y*(0.10750e-3+y*0.74e-5))))));
  } else {
    y=2.0/inArgument;
    ans=(exp(-inArgument)/sqrt(inArgument))*(1.25331414+y*(-0.7832358e-1
			+y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
			+y*(-0.251540e-2+y*0.53208e-3))))));
  }
  return ans;
}

double bessi0(double inArgument) //aus nr-code 
{
  double ax,ans;
  double y;
  
  if ((ax=fabs(inArgument)) < 3.75) {
    y=inArgument/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
  } else {
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
  }
  return ans;
}
