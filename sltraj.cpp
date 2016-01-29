#include <string>
#include <iostream>
#include <iomanip>  
#include <fstream>
#include <sstream>
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
#include "optionparser.h"
double vsink;


struct Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "%s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }

  static option::ArgStatus Unknown(const option::Option& option, bool msg)
  {
    if (msg) printError("Unknown option '", option, "'\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus NonEmpty(const option::Option& option, bool msg)
  {
    if (option.arg != 0 && option.arg[0] != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a non-empty argument\n");
    return option::ARG_ILLEGAL;
  }

  static option::ArgStatus Numeric(const option::Option& option, bool msg)
  {
    char* endptr = 0;
    if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
    if (endptr != option.arg && *endptr == 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires a numeric argument\n");
    return option::ARG_ILLEGAL;
  }
};


enum optionIndex { UNKNOWN, HELP, INPUTFILE, TIMEPARAMS, VEQUATION, SVELOCITY, RANDOM};
enum optionType {OPTIONAL, REQUIRED};
const option::Descriptor usage[] = 
{
  {UNKNOWN, OPTIONAL,"", "",       Arg::Unknown, "USAGE: ./sltraj [options]\n\n"
                                          "Options:" },
  {HELP,    OPTIONAL,"h", "h",  Arg::None,      "  -h \t \t Print usage and exit." },
  {INPUTFILE,REQUIRED,"f","f",  Arg::Required,   "  -f <arg>,\t \t Input file." },
  {TIMEPARAMS,REQUIRED,"t","t", Arg::Required, "  -t <arg=dd/mm/yyyy:tau:intstep>,\t\t Time parameters. " },
  {VEQUATION,REQUIRED,"e","e", Arg::Required, "  -e <arg>,\t \t Velocity equation." },
  {SVELOCITY, REQUIRED,"s","s", Arg::Required,"  -s <arg>,\t \t Requires a number as argument." },
  {RANDOM,OPTIONAL,"r","r", Arg::None,        "  -r, \t \t Enable random displacement equation." },
  { UNKNOWN, OPTIONAL,"", "", Arg::None, "\nExamples:\n"
  " ./sltraj -t =01/02/1989:10:0.01 -f yeah.dat -s 12 -e 0 \n" 
  " ./sltraj -t =01/02/1989:10:0.01 -f yeah.dat -s 12 -e 0 -r \n" 
  },
  { 0, 0, 0, 0, 0, 0} 
};



int main(int argc, char* argv[])
{

  /**********************************************
   * READ COMAND LINE PARAMETERS
   **********************************************/

  string fnameitracers;
  date startdate;
  double tau;
  double intstep;
  int eqvelocity;
  int random=0;
  int verbose=0;
  /* Default parameters */
  string wdir="/home/pmonroy/ESCOLA/RUNS/SLtrajectories/";

  /* VELOCITY FLOW PARAMS */
  char velocitydir[] = "/scratch/pmonroy/";
  //char velocitydir[] = "/data/geo/escola/roms_benguela/";
  date reference_date = {8,  //year
			 1,  //month
			 1,  //day
  };
  int (*velocity)(double ,vectorXYZ , vectorXYZ* );
 
  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present
  option::Stats  stats(usage, argc, argv);
  vector<option::Option> options(stats.options_max);
  vector<option::Option> buffer(stats.buffer_max);
  option::Parser parse(usage, argc, argv, &options[0], &buffer[0]); // Parsing options

  /**************************************************************
   * CHECK COMAND LINE OPTIONS
   *************************************************************/
  if (parse.error())// Error parsing
    return 1;

  if (options[HELP] || argc == 0) // Help options or no options given
    {
      option::printUsage(std::cout, usage);
      return 0;
    }

  for(unsigned int i=0; i<stats.options_max; i++)
    {
      if(options[i].type()!=usage[i].type)// Options that are required 
	{
	  if(options[i].count()==0)
	    cout<<"Option '"<< usage[i].shortopt << "' is required."<<endl;
	  return 0;
	}
      if(options[i].count()>1) // Repeated options are not permited
	{
	  cout<<"Option '"<< usage[i].shortopt << "' is repeated."<<endl;
	  return 0;
	}
    }


  /**************************************************************
   * PARSING ARGUMENTS
   *************************************************************/ 

  if (options[INPUTFILE])
    fnameitracers = options[INPUTFILE].arg;

  if (options[TIMEPARAMS])
    {
     stringstream ss(options[TIMEPARAMS].arg); 
     if (ss.peek() == '=')
       ss.ignore();

     ss >> startdate.day;
     
     if (ss.peek() == '/')
       ss.ignore();
     
     ss >> startdate.month;
     
     if (ss.peek() == '/')
       ss.ignore();
     
     ss >> startdate.year;

    if (ss.peek() == ':')
       ss.ignore();
     
    ss >> tau;

    if (ss.peek() == ':')
       ss.ignore();
     
    ss >> intstep;
    }
  if (options[VEQUATION])
    {
      stringstream ss(options[VEQUATION].arg); 
      ss >> eqvelocity;
    }
  
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
    default:
      cout<<"Invalid argument of vequation"<<endl;
      return 0;
    }
  
  if (options[SVELOCITY])
    {
      stringstream ss(options[SVELOCITY].arg); 
      ss >> vsink;
    }

  
  if (options[RANDOM])
     random = 1;

  /* VERBOSE OPTIONS */

  cout << "stardate: " <<startdate.day<<" "<<startdate.month<<" "<<startdate.year <<endl;
  cout << "tau: " << tau <<endl;
  cout << "intstep: " << intstep <<endl;
  if (random)
    cout << "random: enabled "<<endl; 
  else
    cout << "random: disabled"<<endl;
  cout << "vsink: " << vsink <<endl;
  cout << "vequation: " << eqvelocity <<endl;

  
  /**********************************************
   * READ TRACER INITIAL POSTIONS 
   **********************************************/
  double x, y, z;
  vector<vectorXYZ> itracer;
  unsigned int numtracers;
  ifstream fileitracers(fnameitracers.c_str());  
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

  if(LoadVelocities(startdate, (int) (tau), velocitydir)!=0)
    {
      cout << "Error in reading velocities"<< endl;
      return 1;
    }
  cout << "Reading vflow"<< endl;

  int indexextension = fnameitracers.find_last_of("."); 
  int indexpath = fnameitracers.find_last_of("/")+1; 
  string rawfnameitracers = fnameitracers.substr(indexpath, indexextension-indexpath); 

  vectorXYZ **tracer;

  unsigned int ntau;

  ntau = ((unsigned int)(abs(tau)/intstep)+1);

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

   string nameotracersvtk;  
   ofstream fileotracersvtk;

   string nameoutsider;  
   ofstream fileoutsider;

   double tstart;
   double tend;
   double h;

   int ascnd;


   ascnd = tau > 0;

   if(ascnd)
    {
      tstart = 0.0;
      tend = tau-intstep;
      h=intstep;
    }
  else
    {
      tend = 0.0+intstep;
      tstart = abs(tau);
      h=-1.0*intstep;
    }

   double t;
   int count;
   if(random==0)
     {
       for(t=tstart,count=0; ((t<tend)==ascnd) || (t==tend); t+=h,count++)
	 {
	   cout << t+h << " "<< count+1 <<endl;
	   nameotracers = wdir + rawfnameitracers + 
	     "_t"+ DoubletoString(2, 0, startdate.day) + 
	     DoubletoString(2, 0, startdate.month) + 
	     DoubletoString(4, 0, startdate.year) +
	     DoubletoString(3, 0, tau) + 
	     DoubletoString(3, 2, intstep) + 
	     "_v"+ DoubletoString(3, 0,vsink) +
	     "_e"+ DoubletoString(1, 0,eqvelocity) +
	     "_r"+ DoubletoString(1, 0,random) + 
	     "-"+ DoubletoString(5, 2, t+h) + 
	     ".trac";
	   fileotracers.open(nameotracers.c_str());
	
	   nameotracersvtk = wdir + rawfnameitracers + "_xyzpos" +
	     "_t"+ DoubletoString(2, 0, startdate.day) + 
	     DoubletoString(2, 0, startdate.month) + 
	     DoubletoString(4, 0, startdate.year) +
	     DoubletoString(3, 0, tau) + 
	     DoubletoString(3, 2, intstep) + 
	     "_v"+ DoubletoString(3, 0,vsink) +
	     "_e"+ DoubletoString(1, 0,eqvelocity) + 
	     "_r"+ DoubletoString(1, 0,random) + 
	     "-"+ DoubletoString(5, 2, t+h) + 
	     ".vtk";
	   fileotracersvtk.open(nameotracersvtk.c_str());
	   fileotracersvtk << "# vtk DataFile Version 3.0" << endl;
	   fileotracersvtk <<  "vtk output" << endl;
	   fileotracersvtk <<  "ASCII " << endl;
	   fileotracersvtk << "DATASET POLYDATA" << endl;
	   fileotracersvtk << "POINTS "<< itracer.size() <<" float" << endl;

	   for (unsigned int i = 0; i <itracer.size() ; i++)
	     {
	       tracer[i][count+1] = tracer[i][count];
	       // Semi-implicit 4th order Runge-Kutta
	       if(outsider[i]==0)
		 outsider[i]=RK4(t, h, &tracer[i][count+1], velocity);
	       
	       fileotracers << setprecision(20) << tracer[i][count+1] << endl;
	       fileotracersvtk << tracer[i][count+1].x<<" "<<tracer[i][count+1].y<<" "<<tracer[i][count+1].z*0.001 << endl;
	     }
	   fileotracersvtk << "VERTICES "<< itracer.size() <<" " << 2*itracer.size() << endl;
	   for(unsigned int i = 0; i<itracer.size(); i++ )
	     {
	       fileotracersvtk  << "1" << " " << i << endl;
	     }
	   fileotracers.close();
	   fileotracersvtk.close();
	 }
     }
   else
     {
       for(t=tstart,count=0; ((t<tend)==ascnd) || (t==tend); t+=h,count++)
	 {
	   cout << t+h << " "<< count+1 <<endl;
	   nameotracers = wdir + rawfnameitracers + 
	     "_t"+ DoubletoString(2, 0, startdate.day) + 
	     DoubletoString(2, 0, startdate.month) + 
	     DoubletoString(4, 0, startdate.year) +
	     DoubletoString(3, 0, tau) + 
	     DoubletoString(3, 2, intstep) + 
	     "_v"+ DoubletoString(3, 0,vsink) +
	     "_e"+ DoubletoString(1, 0,eqvelocity) + 
	     "_r"+ DoubletoString(1, 0,random) + 
	     "-"+ DoubletoString(5, 2, t+h) + 
	     ".trac";
	   fileotracers.open(nameotracers.c_str());	

	   nameotracersvtk = wdir + rawfnameitracers + "_xyzpos" +
	     "_t"+ DoubletoString(2, 0, startdate.day) + 
	     DoubletoString(2, 0, startdate.month) + 
	     DoubletoString(4, 0, startdate.year) +
	     DoubletoString(3, 0, tau) + 
	     DoubletoString(3, 2, intstep) + 
	     "_v"+ DoubletoString(3, 0,vsink) +
	     "_e"+ DoubletoString(1, 0,eqvelocity) +
	     "_r"+ DoubletoString(1, 0,random) +  
	     "-"+ DoubletoString(5, 2, t+h) + 
	     ".vtk";
	   fileotracersvtk.open(nameotracersvtk.c_str());
	   fileotracersvtk << "# vtk DataFile Version 3.0" << endl;
	   fileotracersvtk <<  "vtk output" << endl;
	   fileotracersvtk <<  "ASCII " << endl;
	   fileotracersvtk << "DATASET POLYDATA" << endl;
	   fileotracersvtk << "POINTS "<< itracer.size() <<" float" << endl;

	   for (unsigned int i = 0; i <itracer.size() ; i++)
	     {
	       tracer[i][count+1] = tracer[i][count];
	       // 2th order Heun 
	       if(outsider[i]==0)
		 outsider[i]=heun(t, h, &tracer[i][count+1], velocity, *stddev_consteddydiff);
	       else
		 randomtrial();

	       fileotracers << setprecision(20) << tracer[i][count+1] << endl;
	       fileotracersvtk << tracer[i][count+1].x<<" "<<tracer[i][count+1].y<<" "<<tracer[i][count+1].z*0.001 << endl;
	     }

	   fileotracersvtk << "VERTICES "<< itracer.size() <<" " << 2*itracer.size() << endl;
	   for(unsigned int i = 0; i<itracer.size(); i++ )
	     {
	       fileotracersvtk  << "1" << " " << i << endl;
	     }
	   fileotracers.close();
	   fileotracersvtk.close();
	 }
     } 


       
   nameoutsider = wdir + rawfnameitracers + 
     "_t"+ DoubletoString(2, 0, startdate.day) + 
     DoubletoString(2, 0, startdate.month) + 
     DoubletoString(4, 0, startdate.year) +
     DoubletoString(3, 0, tau) + 
     DoubletoString(3, 2, intstep) + 
     "_v"+ DoubletoString(3, 0,vsink) +
	 "_e"+ DoubletoString(1, 0,eqvelocity) +
     "_r"+ DoubletoString(1, 0,random) +  
     ".outsider";
   fileoutsider.open(nameoutsider.c_str());
   for(unsigned int i = 0; i<itracer.size(); i++ )
     {
       fileoutsider  << outsider[i]<< endl;
     }
   fileoutsider.close();   

   /****************************************************************************************************
    * STATISTICS
    ****************************************************************************************************/
   vectorXYZ delta;
   vectorXYZ scalefactor;
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
         "_r"+ DoubletoString(1, 0,random) + 
	 ".disp";


   ofstream fileodispersion(fnameodispersion.c_str());

   int ninsiders;
   for(t=tstart,j=0; ((t<tend)==ascnd) || (t==tend); t+=h,j++)
     {       
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
	       	       
	       scalefactor.x = rearth*cos(rads*tracer[i][0].y);
	       scalefactor.y = rearth;
	       scalefactor.z = 1.0;

	       delta = scalefactor*delta;
	       sumdelta += delta;

	       delta *=delta;
	       sum2delta += delta;
	     }
	 }


       meandisplacement = sumdelta*(1.0/(double) ninsiders);
       meandisplacement2 = sum2delta*(1.0/(double) ninsiders);

       dispersion = meandisplacement2 - meandisplacement*meandisplacement;

       fileodispersion << t << " ";
       fileodispersion << setprecision(20) << meandisplacement.x <<" "<< meandisplacement.y <<" " << meandisplacement.z <<" ";
       fileodispersion << setprecision(20) << meandisplacement2.x <<" "<< meandisplacement2.y <<" " << meandisplacement2.z <<" ";
       fileodispersion << setprecision(20) << sqrt(dispersion.x+dispersion.y) <<" "<< sqrt(dispersion.z) <<endl;
     }
   fileodispersion.close();  

   /* VERBOSE */

   cout << "num. part. insiders = "<< ninsiders <<endl;
   cout << "num. part. seeding = "<< numtracers <<endl;

   /****************************************************************************************************
    * VTK FILES
    ****************************************************************************************************/
 
   string namevtktracers;

   namevtktracers = wdir + rawfnameitracers + "_path" +
     "_t"+ DoubletoString(2, 0, startdate.day) + 
     DoubletoString(2, 0, startdate.month) + 
     DoubletoString(4, 0, startdate.year) +
     DoubletoString(3, 0, tau) + 
     DoubletoString(3, 2, intstep) + 
     "_v"+ DoubletoString(3, 0,vsink) +
     "_e"+ DoubletoString(1, 0,eqvelocity) + 
     "_r"+ DoubletoString(1, 0,random) + 
     ".vtk";

    ofstream vtktracers(namevtktracers.c_str());

   vtktracers << "# vtk DataFile Version 3.0" <<endl;
   vtktracers << "yeah!!!!!!!!!!!!!!!!!!!!!!" <<endl;
   vtktracers << "ASCII" <<endl;
   vtktracers << "DATASET POLYDATA" <<endl;
   vtktracers << "POINTS "<< itracer.size()*ntau <<" float"<<endl;
   for (unsigned int i = 0; i <itracer.size() ; i++)
	 {
	   vtktracers << tracer[i][0].x <<" "<< tracer[i][0].y<<" "<< tracer[i][0].z*0.001 <<endl;
	 }
   for(t=tstart,j=0; ((t<tend)==ascnd) || (t==tend); t+=h,j++) 
     {
       for (unsigned int i = 0; i <itracer.size() ; i++)
	 {
	   vtktracers << tracer[i][j+1].x <<" "<< tracer[i][j+1].y<<" "<< tracer[i][j+1].z*0.001 <<endl;
	 }
     }

   vtktracers << "LINES "<< itracer.size() <<" "<< (itracer.size())*(ntau + 1) <<endl;
   for(unsigned int i=0; i< itracer.size(); i++)
     {
       vtktracers << ntau <<" ";
       for(unsigned int j=0; j<itracer.size()*ntau ; j+=itracer.size())
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
   vtktracers << "POINT_DATA "<< itracer.size()*ntau <<endl;
   vtktracers <<"SCALARS lines_color float"<<endl;
   vtktracers <<"LOOKUP_TABLE default"<<endl;
   for(unsigned int i=0; i< itracer.size(); i++)
	 {
	   vtktracers << tstart <<" ";
	 }
   for(t=tstart; ((t<tend)==ascnd) || (t==tend); t+=h) 
     {
       for(unsigned int i=0; i< itracer.size(); i++)
	 {
	   vtktracers << t+h <<" ";
	 }
     }
   vtktracers.close();


   /**************************************
    * FREE MEMORY 
    *************************************/
 
   FreeMemoryVelocityGrid();
   FreeMemoryVelocities((int) (tau));   
   cout << "free vectors "<<endl;   
   for (unsigned int i = 0; i < itracer.size(); i++)
    {
      delete[] tracer[i];
    }
   delete[] tracer;

   return 0;
}
