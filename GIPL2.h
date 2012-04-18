#ifndef GIPL2_H
#define GIPL2_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iomanip>
using namespace std;


#define TRUE 1
#define FALSE 0

extern int numTimeSteps;

class GIPL
{ 
public:
  double *snowDepth;    
  double *snowDensity;
  double *airTemp;
  double *depth;
  double *T_init;    
  double *D_init;
  double *Q1;    
  double *P1;
//  double **T;    
  
GIPL();
~GIPL();
void init(void);
void run(void);
double interp_1(double x, double xi[], double yi[], int imax);
double alfaSnow (int j0, double rho_snow, double d_snow);
double soilThermalConductivity(double depth,double temper);
int numLayers(double depth);
double heatCapacity (double depth0,double temper);
double unfrWater(double temper, double ac, double bc, double cc);
};
#endif
