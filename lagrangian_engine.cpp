#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <omp.h>
#include <ctime>
using namespace std;

#include "velocity.h" // Function to read velocities 
#include "ioutil.h"
#include "constants.h"
#include "random.h"

extern double vsink;

/* Global constant: Eddy diffusivity*/
double Dh = 10.0*secondsday;
double Dz = 1.0e-5*secondsday;
 
int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));


int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ))
{
  vectorXYZ point2,point3,point4;
  vectorXYZ v1,v2,v3,v4;
  
  double tstep2,tstep6;
  double h; // scale factor of equally spaced cordinate system
  double t;

  /* Time increments */
  tstep2 = intstep*0.5;
  tstep6 = intstep/6.0;
  
  /* Calculate V1: */
  if(velocity(t0,*point, &v1))
    return 1;
  h = rearth * cos(rads*(point->y));
  v1.x = degrees*(v1.x / h ); // rads velocity
  v1.y = degrees*(v1.y / rearth); // rads velocity


  /* Calculate V2: */
  t = t0 + tstep2;
  point2 = *point + (tstep2 * v1);

  if(velocity(t,point2, &v2))
    return 1;

  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity


  /* Calculate V3: */
  point3 = *point + (tstep2 * v2);

  if(velocity(t,point3, &v3))
    return 1;

  h = rearth * cos(rads*(point3.y));
  v3.x = degrees*(v3.x / h);
  v3.y = degrees*(v3.y / rearth);

  
  /* Calculate V4: */
  t = t0 + intstep;
  point4 = *point + (intstep * v3);
  
  if(velocity(t,point4, &v4))
    return 1;

  h = rearth * cos(rads*(point4.y));
  v4.x = degrees*(v4.x / h);
  v4.y = degrees*(v4.y / rearth);

  /* Calculate Final point */  
  *point += (tstep6 * (v1 + v4 + 2.0*v2 + 2.0*v3));

  return 0;
}

int heun(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ),vectorXYZ (*stdevnoise)(double ,vectorXYZ))
{
  vectorXYZ point1,point2;
  vectorXYZ v1,v2;
  vectorXYZ g1,g2;
  vectorXYZ u, uh, aux;
  vectorXYZ k,l;
  double h;
  static int rset=0;
  static unsigned long seed=329571112321;

  if(rset == 0)
    {
      rset = 1;
      srand64(seed);
      init_genrand(seed);
    }

  u.x = gasdev2(genrand_real2); 
  u.y = gasdev2(genrand_real2); 
  u.z = gasdev2(genrand_real2); 
  
  /* Calculate V1: */
  point1 = *point;
  if(velocity(t0, point1, &v1))
    return 1;
  h = rearth * cos(rads*(point1.y));
  v1.x = degrees*(v1.x / h ); 
  v1.y = degrees*(v1.y / rearth);

  /* Calculate G1*/
  g1 = stdevnoise(t0, point1);
  g1.x = degrees*(g1.x / h );
  g1.y = degrees*(g1.y / rearth);
 

  /* Calculate V2: */
  uh=sqrt(intstep)*u;
  aux = intstep * v1 + uh * g1;
  point2 = point1 + aux;

  if(velocity(t0+intstep,point2, &v2))
    return 1;
  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h ); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity

  /* Calculate G2*/
  g2 = stdevnoise(t0+intstep, point2);
  g2.x = degrees*(g2.x / h ); // rads velocity
  g2.y = degrees*(g2.y / rearth); // rads velocity

  /* Compute the new position */
  *point += 0.5*(aux + intstep*v2+ uh*g2);

  return 0;
}

void randomtrial(void)
{
  vectorXYZ u;
  
  u.x = gasdev2(genrand_real2); 
  u.y = gasdev2(genrand_real2); 
  u.z = gasdev2(genrand_real2); 
}

int FlowVplusSinkV(double t,vectorXYZ point, vectorXYZ *vint)
{
  if(GetVelocity( t, point, vint))
    return 1;
  
  vint->z = vint->z - vsink;
  return 0;
}

int vmesoscale(double t,vectorXYZ point, vectorXYZ *vint)
{
  
  if(GetVelocity( t, point, vint))
    return 1;
  
  vectorXYZ vsinking;
  vsinking=vectorXYZ(0.0,0.0,-vsink);

  double factor=vsink/gravity;

  vectorXYZ omega;
  omega=vectorXYZ(0.0,2.0*fearth*cos(rads*point.y),2.0*fearth*sin(rads*point.y));

  *vint = *vint + vsinking - factor * cross(omega,*vint);
  
  return 0;
}

int vsubmesoscale(double t,vectorXYZ point, vectorXYZ *vint)
{
  double h=0.1;

  if(GetVelocity( t, point, vint))
    return 1;

  vectorXYZ aflow;

  if(t<h)
    {
      vectorXYZ p1;
      p1 = point;
      if(RK4(t, 0.1, &p1, GetVelocity))
	return 1;
      
      vectorXYZ v1;
      if(GetVelocity(t+0.1, p1, &v1))
	return 1;
      
      vectorXYZ sphomega;
      
      sphomega.x = -(*vint).y/rearth;
      sphomega.y = (*vint).x/rearth;
      sphomega.z = ((*vint).x*tan(rads*point.y))/rearth;

      aflow = (v1-(*vint))*(10.0) + cross(sphomega,*vint);
    }
  else
    {
      vectorXYZ p1;
      p1 = point;
      if(RK4(t, 0.1, &p1, GetVelocity))
	return 1;

      vectorXYZ v1;
      if(GetVelocity(t+0.1, p1, &v1))
	return 1;

      vectorXYZ p2;
      p2 = point;
      if(RK4(t, 0.1, &p1, GetVelocity))
	return 1;

      vectorXYZ v2;
      if(GetVelocity(t-0.1, p2, &v2))
	return 1;     

      vectorXYZ sphomega;
      
      sphomega.x = -(*vint).y/rearth;
      sphomega.y = (*vint).x/rearth;
      sphomega.z = ((*vint).x*tan(rads*point.y))/rearth;
      
      aflow = (v1-v2)*(5.0) + cross(sphomega,*vint);
    }

  vectorXYZ vsinking;
  vsinking=vectorXYZ(0.0,0.0,-vsink);
  
  double factor=vsink/gravity;

  *vint = *vint + vsinking - factor * aflow;
  
  return 0;
}
int vfullinertia(double t,vectorXYZ point, vectorXYZ *vint)
{
  
  if(GetVelocity( t, point, vint))
    return 1;
  

  vectorXYZ p1;
  p1 = point;

  if(RK4(t, 1.0, &p1, GetVelocity))
    return 1;

  vectorXYZ v1;
  if(GetVelocity( t+1.0, p1, &v1))
    return 1;
  

  vectorXYZ aflow,sphomega;

  sphomega.x = -(*vint).y/rearth;
  sphomega.y = (*vint).x/rearth;
  sphomega.z = ((*vint).x*tan(rads*point.y))/rearth;

  aflow = (v1-(*vint)) + cross(sphomega,*vint);

  vectorXYZ vsinking;
  vsinking=vectorXYZ(0.0,0.0,-vsink);

  double factor=vsink/gravity;

  vectorXYZ omega;
  omega=vectorXYZ(0.0,2.0*fearth*cos(rads*point.y),2.0*fearth*sin(rads*point.y));

  
  *vint = *vint + vsinking - factor * (cross(omega,*vint) + aflow);
  return 0;
}

vectorXYZ stddev_consteddydiff(double t,vectorXYZ point)
{

  vectorXYZ stddev;

  stddev.x = sqrt(2.0*Dh);
  stddev.y = stddev.x;
  stddev.z = sqrt(2.0*Dz);

  return stddev;
}
