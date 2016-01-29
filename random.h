#ifndef RANDOM
#define RANDOM
double drand64(void) ;
void srand64(unsigned long long int seed);
double gasdev(double (*rand)(void ));
double gasdev2(double (*rand)(void ));
void init_genrand(unsigned long s);
double genrand_real2(void);
#endif
