#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <string>
using namespace std;
#include "vectorXYZ.h" 
#include "constants.h" 

double frand(double fmin, double fmax);

int main(int argc, char **argv)
{
  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/
  string ofname;
  vectorXYZ domainmin;
  vectorXYZ domainmax;
  int ni, nj, nk, n;

  if(argc!= 11)
    {
      cout << "incorrect number of parameters" <<endl;
      return 1;
    }

  domainmin.x = atof(argv[1]);
  domainmin.y = atof(argv[2]);
  domainmin.z = atof(argv[3]);

  domainmax.x = atof(argv[4]);
  domainmax.y = atof(argv[5]);
  domainmax.z = atof(argv[6]);

  ni = atof(argv[7]);
  nj = atof(argv[8]);
  nk = atof(argv[9]);

  n =ni*nj*nk;

  ofname = argv[10];
  ofstream ofile(ofname.c_str());

  cout << " domain min: "<< domainmin<<endl;
  cout << " domain max: "<< domainmax<<endl;
  cout << " ni, nj, nk: "<< ni<<" "<<nj<<" "<<nk<<endl;
  cout << " n: "<< n<<endl;
  cout << "ofname:" << ofname <<endl;

  // RANDOM Positions in the domain limits
  srand(time(NULL));
  for(int q=0; q<n; q++)
    {
      ofile<< frand(domainmin.x, domainmax.x)<<" ";
      ofile<< degrees*asin(frand(sin(rads*domainmin.y), sin(rads*domainmax.y)))<<" ";
      ofile<< frand(domainmin.z, domainmax.z)<<endl;
    }
  return 0;
}

double frand(double fmin, double fmax)
{
    double f = (double)rand() / RAND_MAX;
    return fmin + f * (fmax - fmin);
}
