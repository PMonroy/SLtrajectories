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
#include "ioutil.h"
#include "constants.h"
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
  {VEQUATION,REQUIRED,"e","e", Arg::Required, "  -e <arg,arg>,\t \t Velocity equation comparation." },
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
  int eqvelocity[2];
  int random=0;
  int verbose=0;
  /* Default parameters */
  string wdir="/home/pmonroy/ESCOLA/RUNS/SLtrajectories/";

  /* VELOCITY FLOW PARAMS */
  //char velocitydir[] = "/scratch/pmonroy/";
  //char velocitydir[] = "/data/geo/escola/roms_benguela/";
  //date reference_date = {8,  //year
  //			 1,  //month
  //			 1,  //day
  //};
 
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

  int indexextension = fnameitracers.find_last_of("."); 
  int indexpath = fnameitracers.find_last_of("/")+1; 
  string rawfnameitracers = fnameitracers.substr(indexpath, indexextension-indexpath); 

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
      ss >> eqvelocity[0];

      if (ss.peek() == ',')
       ss.ignore();
     
      ss >> eqvelocity[1];

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
  cout << "vequation: " << eqvelocity[0]<< " "<< eqvelocity[1]<<endl;


  /**********************************************
   * READ outsider files
   **********************************************/
 
  vector<int> outsider1;
  string nameoutsider1;  
  ifstream fileoutsider1;
  int i;

  nameoutsider1 = wdir + rawfnameitracers + 
    "_t"+ DoubletoString(2, 0, startdate.day) + 
    DoubletoString(2, 0, startdate.month) + 
    DoubletoString(4, 0, startdate.year) +
    DoubletoString(3, 0, tau) + 
    DoubletoString(3, 2, intstep) + 
    "_v"+ DoubletoString(3, 0,vsink) +
    "_e"+ DoubletoString(1, 0,eqvelocity[0]) +
    "_r"+ DoubletoString(1, 0,random) +  
    ".outsider";
  fileoutsider1.open(nameoutsider1.c_str());
  
  while( fileoutsider1 >> i )
    {
     outsider1.push_back(i);
    }  
  fileoutsider1.close();



  cout << "num. tracer=" << outsider1.size() << endl;


  vector<int> outsider2;
  string nameoutsider2;  
  ifstream fileoutsider2;

  nameoutsider2 = wdir + rawfnameitracers + 
    "_t"+ DoubletoString(2, 0, startdate.day) + 
    DoubletoString(2, 0, startdate.month) + 
    DoubletoString(4, 0, startdate.year) +
    DoubletoString(3, 0, tau) + 
    DoubletoString(3, 2, intstep) + 
    "_v"+ DoubletoString(3, 0,vsink) +
    "_e"+ DoubletoString(1, 0,eqvelocity[1]) +
    "_r"+ DoubletoString(1, 0,random) +  
    ".outsider";
  fileoutsider2.open(nameoutsider2.c_str());
  
  while( fileoutsider2 >> i )
    {
     outsider2.push_back(i);
    }  
  fileoutsider2.close();



  cout << "num. tracer=" << outsider2.size() << endl;

  vector<int> outsider;

  for(i=0; i<outsider1.size(); i++)
    {
      if(outsider1[i]==0 && outsider2[i]==0)
	outsider.push_back(0);
      else
	outsider.push_back(1);
    }

  cout << "num. tracer=" << outsider.size() << endl;



  /**********************************************
   * STATISTICS
   **********************************************/
 
  string nametracers1;  
  ifstream filetracers1;
  
  string nametracers2;  
  ifstream filetracers2;

  string nameoutput;  
  ofstream fileoutput;

  double t,h;

  h=intstep;

  vector<vectorXYZ> tracer;
  vector<vectorXYZ> itracer;

  vectorXYZ vect0,vect1;


  vectorXYZ delta;
  vectorXYZ scalefactor;
  vectorXYZ sumPos, sum2Pos, sumdelta,sum2delta;
  vectorXYZ meanPos, mean2Pos, varPos, meandisplacement,meandisplacement2;

  double dist2;
  double sumDist, sum2Dist;
  double meanDist, mean2Dist, varDist;

  double dist2Hor;
  double sumDistHor, sum2DistHor;
  double meanDistHor, mean2DistHor, varDistHor;

  double dist2Ver;
  double sumDistVer, sum2DistVer;
  double meanDistVer, mean2DistVer, varDistVer;

  // NAME FILE OUTPUT

  nameoutput = wdir + rawfnameitracers + 
    "_t"+ DoubletoString(2, 0, startdate.day) + 
    DoubletoString(2, 0, startdate.month) + 
    DoubletoString(4, 0, startdate.year) +
    DoubletoString(3, 0, tau) + 
    DoubletoString(3, 2, intstep) + 
    "_v"+ DoubletoString(3, 0,vsink) +
    "_e"+ DoubletoString(1, 0,eqvelocity[0]) + DoubletoString(1, 0,eqvelocity[1]) +
    "_r"+ DoubletoString(1, 0,random) +
    ".rdist";

  fileoutput.open(nameoutput.c_str());

  for(t=0; t<tau; t+=h)
    {
      nametracers1 = wdir + rawfnameitracers + 
	"_t"+ DoubletoString(2, 0, startdate.day) + 
	DoubletoString(2, 0, startdate.month) + 
	DoubletoString(4, 0, startdate.year) +
	DoubletoString(3, 0, tau) + 
	DoubletoString(3, 2, intstep) + 
	"_v"+ DoubletoString(3, 0,vsink) +
	"_e"+ DoubletoString(1, 0,eqvelocity[0]) + 
	"_r"+ DoubletoString(1, 0,random) + 
	"-"+ DoubletoString(5, 2, t+h) + 
	".trac";
      filetracers1.open(nametracers1.c_str());	
      
      nametracers2 = wdir + rawfnameitracers + 
	"_t"+ DoubletoString(2, 0, startdate.day) + 
	DoubletoString(2, 0, startdate.month) + 
	DoubletoString(4, 0, startdate.year) +
	DoubletoString(3, 0, tau) + 
	DoubletoString(3, 2, intstep) + 
	"_v"+ DoubletoString(3, 0,vsink) +
	"_e"+ DoubletoString(1, 0,eqvelocity[1]) + 
	"_r"+ DoubletoString(1, 0,random) + 
	"-"+ DoubletoString(5, 2, t+h) + 
	".trac";
      filetracers2.open(nametracers2.c_str());	
      
      
      for(i=0; i<outsider.size(); i++)
	{
	  filetracers1 >> vect0;
	  filetracers2 >> vect1;
      
	  if(outsider[i]==0)
	    {
	      tracer.push_back(vect0);
	      itracer.push_back(vect1);
	    }
	}

      sumPos.x = 0.0;
      sumPos.y = 0.0;
      sumPos.z = 0.0; 

      sum2Pos.x = 0.0;
      sum2Pos.y = 0.0;
      sum2Pos.z = 0.0;

      sumdelta.x = 0.0;
      sumdelta.y = 0.0;
      sumdelta.z = 0.0; 
      
      sum2delta.x = 0.0;
      sum2delta.y = 0.0;
      sum2delta.z = 0.0; 

      sumDist = 0.0; 
      sum2Dist = 0.0;

      sumDistHor = 0.0; 
      sum2DistHor = 0.0;

      sumDistVer = 0.0; 
      sum2DistVer = 0.0;

      for (i=0; i<tracer.size(); i++)
	{
	  sumPos += tracer[i];
	  sum2Pos += (tracer[i]*tracer[i]);

	  delta = itracer[i]-tracer[i];
	  
	  delta.x = rads * delta.x;
	  delta.y = rads * delta.y;
	  
	  scalefactor.x = rearth*cos(rads*tracer[i].y);
	  scalefactor.y = rearth;
	  scalefactor.z = 1.0;
	  
	  delta = scalefactor*delta;
	  //sumdelta += delta;
	  
	  delta *=delta;
	  sum2delta += delta;

	  dist2 = delta.x + delta.y + delta.z; 
	  dist2Hor = delta.x + delta.y; 
	  dist2Ver = delta.z;

	  sum2Dist += dist2;
	  sumDist += sqrt(dist2);

	  sum2DistHor += dist2Hor;
	  sumDistHor += sqrt(dist2Hor);

	  sum2DistVer += dist2Ver;
	  sumDistVer += sqrt(dist2Ver);
	}
      

      meanPos = sumPos*(1.0/(double) tracer.size());
      mean2Pos = sum2Pos*(1.0/(double) tracer.size());
      varPos = mean2Pos - (meanPos*meanPos);

      meanDist = sumDist*(1.0/(double) tracer.size());
      mean2Dist = sum2Dist*(1.0/(double) tracer.size());
      varDist = mean2Dist - (meanDist*meanDist);

      meanDistHor = sumDistHor*(1.0/(double) tracer.size());
      mean2DistHor = sum2DistHor*(1.0/(double) tracer.size());
      varDistHor = mean2DistHor - (meanDistHor*meanDistHor);
      
      meanDistVer = sumDistVer*(1.0/(double) tracer.size());
      mean2DistVer = sum2DistVer*(1.0/(double) tracer.size());
      varDistVer = mean2DistVer - (meanDistVer*meanDistVer);

      //meandisplacement = sumdelta*(1.0/(double) tracer.size());
      meandisplacement2 = sum2delta*(1.0/(double) tracer.size());
      
      fileoutput << t+h <<" ";
      fileoutput << meanPos <<" ";
      fileoutput << sqrt(varPos.x/(tracer.size() - 1.0)) <<" ";
      fileoutput << sqrt(varPos.y/(tracer.size() - 1.0)) <<" ";
      fileoutput << sqrt(varPos.z/(tracer.size() - 1.0)) <<" ";
      fileoutput << meanDist <<" ";
      fileoutput << sqrt(varDist/(tracer.size() - 1.0)) <<" ";
      fileoutput << meanDistHor <<" ";
      fileoutput << sqrt(varDistHor/(tracer.size() - 1.0)) << " ";
      fileoutput << meanDistVer <<" ";
      fileoutput << sqrt(varDistVer/(tracer.size() - 1.0)) << " ";
      fileoutput << tracer.size();
      fileoutput <<endl;

      itracer.clear();
      tracer.clear();
      
      filetracers1.close();
      filetracers2.close();
    }
  fileoutput.close();
  return 0;
}
