#include <string>
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace std;
#include "date.h"
#include "vectorXYZ.h"
#include "velocity.h"
#include "ioutil.h"
#include "constants.h"
#include "lagrangian_engine.h"

double vsink;

int main(int argc, char **argv)
{

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/
  int opt;
  int optflag=0;

  string fnameitracers;
  ifstream fileitracers;
  date startdate;
  double tau;
  double intstep;
  int eqvelocity;

  char velocitydir[] = "/scratch/pmonroy/";
  //char velocitydir[] = "/data/geo/escola/roms_benguela/";

  string wdir="/home/pmonroy/ESCOLA/RUNS/SLtrajectory/";

  date reference_date = {8,  //year
			 1,  //month
			 1,  //day
  };

  int (*velocity)(double ,vectorXYZ , vectorXYZ* );
 

  while ((opt = getopt(argc, (char **)argv, "f:t:v:e:")) != -1) 
    {
      cout << "Option: " << (char)opt;
      switch(opt)
	{
	case 'f':
	  optflag++;
	  fnameitracers=optarg;
	  fileitracers.open(fnameitracers.c_str());
	  if (optarg)
	    cout << ", argument: " << optarg;
	  else
	    cout << ",no argument! ";
	  cout << '\n';
	  break;
	case 't':
	  optflag++;
	  if(sscanf(optarg,"%u-%u-%u:%lf:%lf",&startdate.day,&startdate.month,&startdate.year,&tau,&intstep)!=5)
	    {
	      cout << ": Date format in startdate is incorrect" << endl;
	      return 1;
	    }
	  else
	    cout << ", arguments: " <<startdate.day<<" "<<startdate.month<<" "<<startdate.year <<endl;
	  break;
	case 'v':
	  optflag++;
	  vsink=atof(optarg); 
	  if (optarg)
	    cout << ", argument: " << vsink;
	  else
	    cout << ",no argument! ";
	  cout << '\n';
	  break;
	case 'e':
	  optflag++;
	  eqvelocity=atoi(optarg); 
	  if (optarg)
	    cout << ", argument: " << eqvelocity;
	  else
	    cout << ",no argument! ";
	  cout << '\n';
	  break;
	case '?':
	  cout << "Unknown option "<<endl ;
	}
    }

  cout << "number of options =" << optflag<<endl;
  
  switch(eqvelocity)
    {
    case 0:
      velocity = GetVelocity;  
      break;
    case 1:
      velocity = FlowVplusSinkV;
      break;
    case 2:
      velocity = vmesoscale;
      break;
    case 3:
      velocity = vsubmesoscale;
      break;  
    case 4:
      velocity = vfullinertia;
      break;
    }

  /**********************************************
   * READ TRACER INITIAL POSTIONS 
   **********************************************/
  double x, y, z;
  vector<vectorXYZ> itracer;
  unsigned int numtracers;

  while( fileitracers >> x >> y >> z )
    {
      itracer.push_back(vectorXYZ(x,y,z));
    }  
  fileitracers.close();

  numtracers = itracer.size();
  

  /**********************************************
   * LAGRANGIAN ENGINE 
   **********************************************/

  // Setup velocity field
  if(LoadVelocityGrid(reference_date, velocitydir)!=0)
    {
      cout << "Error in reading reference date netcdf"<< endl;
      return 1;
    }

  if(LoadVelocities(startdate, (int) tau+1, velocitydir)!=0)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }

  int indexextension = fnameitracers.find_last_of("."); 
  int indexpath = fnameitracers.find_last_of("/")+1; 
  string rawfnameitracers = fnameitracers.substr(indexpath, indexextension-indexpath); 

  vectorXYZ **tracer;

  unsigned int ntau;

  ntau = ((unsigned int)(tau/intstep)+1);

  tracer = new vectorXYZ *[itracer.size()];
  for (unsigned int i = 0; i < itracer.size(); i++)
    {
      tracer[i] = new vectorXYZ [ntau];
      tracer[i][0] = itracer[i];
    }

   int outsider[itracer.size()];  
   unsigned int j;

   for(unsigned int i = 0; i <itracer.size() ; i++)
     outsider[i]=0;


   string nameotracers;  
   ofstream fileotracers;

   for(double t=0.0; t<=tau-intstep; t+=intstep)
     {
       
       j = ((unsigned int)(t/intstep)+1);
       
       nameotracers = wdir + rawfnameitracers + 
	 "_t"+ DoubletoString(2, 0, startdate.day) + 
	 DoubletoString(2, 0, startdate.month) + 
	 DoubletoString(4, 0, startdate.year) +
	 DoubletoString(3, 0, tau) + 
	 DoubletoString(3, 2, intstep) + 
	 "_v"+ DoubletoString(3, 0,vsink) +
	 "_e"+ DoubletoString(1, 0,eqvelocity) + 
	 "-"+ DoubletoString(5, 2, t+intstep) + 
	 ".trac";
       fileotracers.open(nameotracers.c_str());	

       for (unsigned int i = 0; i <itracer.size() ; i++)
	 {
	   tracer[i][j] = tracer[i][j-1];
	   // Semi-implicit 4th order Runge-Kutta
	   if(outsider[i]==0)
	     outsider[i]=RK4(t, intstep, &tracer[i][j], velocity);

	   fileotracers << tracer[i][j] << endl;
	 }
       fileotracers.close();
     }
  
   /****************************************************************************************************
    * STATISTICS
    ****************************************************************************************************/

   vectorXYZ delta;
   vectorXYZ h;
   vectorXYZ sumdelta,sum2delta;
   vectorXYZ meandisplacement,meandisplacement2, dispersion;

   string fnameodispersion;

   fnameodispersion = wdir + rawfnameitracers + 
	 "_t"+ DoubletoString(2, 0, startdate.day) + 
	 DoubletoString(2, 0, startdate.month) + 
	 DoubletoString(4, 0, startdate.year) +
	 DoubletoString(3, 0, tau) + 
	 DoubletoString(3, 2, intstep) + 
	 "_v"+ DoubletoString(3, 0,vsink) +
	 "_e"+ DoubletoString(1, 0,eqvelocity) + 
	 ".disp";


   ofstream fileodispersion(fnameodispersion.c_str());

   int ninsiders;

   for(double t=0.0; t<=tau-intstep; t+=intstep)
     {       
       j = (unsigned int)(t/intstep);
       
       sumdelta.x = 0.0;
       sumdelta.y = 0.0;
       sumdelta.z = 0.0; 

       sum2delta.x = 0.0;
       sum2delta.y = 0.0;
       sum2delta.z = 0.0; 

       ninsiders = 0;

       for (unsigned int i = 0; i <numtracers ; i++)
	 {
	   if(outsider[i]==0)
	     {
	       ninsiders++;
	       delta = tracer[i][j] - tracer[i][0];
	       
	       delta.x = rads * delta.x;
	       delta.y = rads * delta.y;
	       	       
	       h.x = rearth*cos(rads*tracer[i][0].y);
	       h.y = rearth;
	       h.z = 1.0;

	       delta = h*delta;
	       sumdelta += delta;

	       delta *=delta;
	       sum2delta += delta;
	     }
	 }


       meandisplacement = sumdelta*(1.0/(double) ninsiders);
       meandisplacement2 = sum2delta*(1.0/(double) ninsiders);

       dispersion = meandisplacement2 - meandisplacement*meandisplacement;

       fileodispersion << t<<" "<<  meandisplacement.x <<" "<< meandisplacement.y <<" " << meandisplacement.z <<" ";
       fileodispersion << setprecision(16) << sqrt(meandisplacement2.x) <<" "<< sqrt(meandisplacement2.y) <<" " << sqrt(meandisplacement2.z) <<" ";
       fileodispersion << setprecision(16) << sqrt(meandisplacement2.x+meandisplacement2.y+meandisplacement2.z) <<" ";
       fileodispersion << setprecision(16) << sqrt(dispersion.x) <<" "<< sqrt(dispersion.y) <<" "<< sqrt(dispersion.z) <<endl;
     }
   fileodispersion.close();  


   /****************************************************************************************************
    * VTK FILE
    ****************************************************************************************************/
 
   string namevtktracers;

    namevtktracers = wdir + rawfnameitracers + 
	 "_t"+ DoubletoString(2, 0, startdate.day) + 
	 DoubletoString(2, 0, startdate.month) + 
	 DoubletoString(4, 0, startdate.year) +
	 DoubletoString(3, 0, tau) + 
	 DoubletoString(3, 2, intstep) + 
	 "_v"+ DoubletoString(3, 0,vsink) +
	 "_e"+ DoubletoString(1, 0,eqvelocity) + 
	 ".vtk";

    ofstream vtktracers(namevtktracers.c_str());

   vtktracers << "# vtk DataFile Version 3.0" <<endl;
   vtktracers << "yeah!!!!!!!!!!!!!!!!!!!!!!" <<endl;
   vtktracers << "ASCII" <<endl;
   vtktracers << "DATASET POLYDATA" <<endl;
   vtktracers << "POINTS "<< itracer.size()*((int)(tau/intstep)+1) <<" float"<<endl;
   for (unsigned int i = 0; i <itracer.size() ; i++)
	 {
	   vtktracers << tracer[i][0].x <<" "<< tracer[i][0].y<<" "<< tracer[i][0].z*0.001 <<endl;
	 }
   for(double t=0.0; t<=tau-intstep; t+=intstep)
     {

       j = ((unsigned int)(t/intstep)+1);

       for (unsigned int i = 0; i <itracer.size() ; i++)
	 {

	   vtktracers << tracer[i][j].x <<" "<< tracer[i][j].y<<" "<< tracer[i][j].z*0.001 <<endl;
	 }
     }

   vtktracers << "LINES "<< itracer.size() <<" "<< (itracer.size())*((int)(tau/intstep) + 2) <<endl;
   for(unsigned int i=0; i< itracer.size(); i++)
     {
       vtktracers << ((int)(tau/intstep) + 1) <<" ";
       for(unsigned int j=0; j<itracer.size()*((int)(tau/intstep)+1) ; j+=itracer.size())
	 { 
	   vtktracers << i+j <<" ";
	 }
       vtktracers << endl;
     }

   vtktracers << "CELL_DATA "<< itracer.size() <<endl;
   vtktracers <<"SCALARS IDtracer float"<<endl;
   vtktracers <<"LOOKUP_TABLE default"<<endl; 
   for(unsigned int i=0; i< itracer.size(); i++)
     {
       vtktracers << i<<" ";
     }
   vtktracers << "POINT_DATA "<< itracer.size()*((int)(tau/intstep)+1) <<endl;
   vtktracers <<"SCALARS lines_color float"<<endl;
   vtktracers <<"LOOKUP_TABLE default"<<endl;
   for(double t=0.0; t<=tau; t+=intstep)
     {
       for(unsigned int i=0; i< itracer.size(); i++)
	 {
	   vtktracers << t<<" ";
	 }
     }
   vtktracers.close();
   
   FreeMemoryVelocityGrid();
   FreeMemoryVelocities(tau);


   
  for (unsigned int i = 0; i < itracer.size(); i++)
    {
       delete[] tracer[i];
    }

  return 0;
}
