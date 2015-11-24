
int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));
int FlowVplusSinkV(double t,vectorXYZ point, vectorXYZ *vint);
int vmesoscale(double t,vectorXYZ point, vectorXYZ *vint);
int vsubmesoscale(double t,vectorXYZ point, vectorXYZ *vint);
int vfullinertia(double t,vectorXYZ point, vectorXYZ *vint);
