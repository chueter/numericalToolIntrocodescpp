#include <math.h>
#include "nrutil.h"
#include <iostream>
#include <stdlib.h>
#define TINY 1.0e-20     /* a small number */
#define PI 3.1415920     /* err... what's this??? */



using namespace std;

double getLatticeSpacing(void);
int getGridLength(void);
double** fillkoeffizientenMatrix(int inNumberGridPoints, double inLatticeSpacing, double inGammaD);
double* initialisiereInhomogenitaet(int inNumberGridPoints, double inLatticeSpacing);
void solveLes(double** koeffizientenMatrix, int inNumberGridPoints, int* indx ,double* inhomogenitaet);
void decomposeMatrixToLAndUTriangularMatrices(double** koeffizientenMatrix, int inNumberGridPoints, int* indx, double* d);

int main(int argc, char* argv)
{
  double gammaD = 0.00; 
  double** koeffizientenMatrix;
  double* inhomogenitaet;
  double d;
  double latticeSpacing = 1;
  int gridLength = 0;
  int* indx;  

  gridLength = getGridLength();
  latticeSpacing = getLatticeSpacing();
  int numberGridPoints = gridLength/latticeSpacing;

  inhomogenitaet = initialisiereInhomogenitaet(numberGridPoints,latticeSpacing);
  koeffizientenMatrix = fillkoeffizientenMatrix(numberGridPoints,latticeSpacing,gammaD);
  indx = ivector(1,numberGridPoints);
  decomposeMatrixToLAndUTriangularMatrices(koeffizientenMatrix, numberGridPoints, indx, &d);
  solveLes(koeffizientenMatrix, numberGridPoints, indx, inhomogenitaet);



  


   FILE* fp;
   fp = fopen("loesung.dat", "w");     
   for (int i=1 ; i <= numberGridPoints; i++)
     {
       //cout << i  << inhomogenitaet[i] << endl;
       fprintf(fp, "%d  %g \n", i, inhomogenitaet[i]);
     }
   fclose(fp);
      
      
   return 0;
}

int getGridLength()
{
  int outGridLength = 0;
  
  cout << "please enter the gridlength:" << endl; 
 
  cin >> outGridLength ;
  
  return outGridLength;
} 

double getLatticeSpacing()
{
  double outLatticeSpacing = 0.0;
  
  cout << "please enter the latticespacing:" << endl;
 
  cin >> outLatticeSpacing;
 
  return outLatticeSpacing;
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
    for (int i=1; i <= inNumberGridPoints; i++)
      {
	fprintf(fp, "%g  %g\n", i*inLatticeSpacing, inhomogenitaet[i]);
      }
    fclose(fp);
   
    return inhomogenitaet;
  }


 double** fillkoeffizientenMatrix(int inNumberGridPoints, double inLatticeSpacing, double gammaD)
  {
    double** koeffizientenMatrix;

    koeffizientenMatrix = dmatrix(1, inNumberGridPoints, 1, inNumberGridPoints);

      {
	for(int i=1; i <= inNumberGridPoints; i++)
	{
	  for(int j=1; j <= inNumberGridPoints; j++)
	  {
           
            //if (j <= i)
	    //{
	    //koeffizientenMatrix [i][j] = i*j;
	    //}
	    //if (j > i)
	    //{
	    //koeffizientenMatrix [i][j] = 0;
	    //}




	    if (j==i)                                 //main diagonal elements      //control matrix output shows that
	    {                                                                     //the i counts rows,j columns, so
	    koeffizientenMatrix[i][j]                                           //allocating elements below the main
	      = (gammaD / inLatticeSpacing) + 2 * sqrt(2*inLatticeSpacing)/PI;  //diagonal to case j<=(i-1) is correct
	    if (i == inNumberGridPoints)
	      {
	        koeffizientenMatrix[i][j] = 2 * sqrt(2 * inLatticeSpacing)/PI;
	      }
	    }
	    else if (j <=(i-1))                       //generic elements below the main diagonal
	    {
	    koeffizientenMatrix[i][j] = sqrt(2 * inLatticeSpacing) / (PI  * sqrt(i-j));
	    }  
	    else if (j == (i+1))                      //elements due to derivative on first upper beside diagonal
	    {
	    koeffizientenMatrix[i][j]= -gammaD/inLatticeSpacing;
	    }
            else if (j > (i+1))                       //elements on second or later upper beside diagonals
	    {
	    koeffizientenMatrix[i][j] = 0;
	    }
	  }
	}
  
       
      }
      
      //cout << "the matrix is: " << endl;                        // for correction purposes
      //for (int i = 1; i <= inNumberGridPoints; i++)
      //{
      //  for (int j = 1; j <= inNumberGridPoints; j++)
      //  cout << koeffizientenMatrix [i] [j] << " ";
      //cout << endl;
      //}

      //cout << "control matrix is: " << endl;                        // for correction purposes :
      //for (int i = 1; i <= inNumberGridPoints; i++)                 // i counts the rows,
      //{                                                           // j the columns
      //  for (int j = 1; j <= inNumberGridPoints; j++)
      //    cout << i << j << " ";
      //  cout << endl;
      //}
      //for (int i=1; i <= inNumberGridPoints; i++)
      //{
      //for (int j=1; j <= inNumberGridPoints; j++)
      //	
      //	  cout << koeffizientenMatrix[i] [j] << " " ;
      //    cout << endl; 
	
      //}
  return koeffizientenMatrix;
  }


void decomposeMatrixToLAndUTriangularMatrices(double** koeffizientenMatrix, int inNumberGridPoints, int* indx,double* d)
{
  int i,j,imax,k;
  double big,dum,sum,temp;
  double *vv;

  vv = dvector(1,inNumberGridPoints);
  cout << "richtich"  << endl;
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

