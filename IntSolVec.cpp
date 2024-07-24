#include <iostream>
#include <fstream>
#include <cmath>

double readVector(ifstream&, char filename[20], int inDim);
int getDim();
int getLength();
double integrateVector(double disc, double* solVector, int inDim);
double integrate(double* solVector, int index, double disc);

int main(int argc, char* argv)
{
  int dim=1;
  int length=1;
  double disc=1.0;

  dim = getDim();
  length = getLength();
  disc = length/dim;

  ///// control //// begin
  cout << disc << " disc "  
       << length << " length " 
       << dim << " dim " << endl;
  ///// control //// end 

  double* solVector;
  solVector = new double[dim];

  

  
  

  return 0;
}
