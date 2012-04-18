#include "GIPL2.h"

int numTimeSteps=1500;
int numSoilLayers=5;
int TStep =24;
//soil properties
double Thick[]= {0.5, 0.5 ,9 ,20, 70};
int Num[]={5, 5, 10, 19, 14};
double Tfr[]={0,0,0,0,0};
double Wvol[]={0.48,0.42,0.31,0.26,0.22};
double Wunf[]={0.02, 0.02, 0.01,0.02,0.02};
double aclv[]={0.005,       0.035,  0.061,  0.018,  0.064};
double bclv[]={-0.1,        -0.32,  -0.35,  -0.17,  -0.340};
double cclv[]={0.,          0,              0.0,    0.0,    0.};;
double Cond_Th[]= {0.34,0.73,1.23,2.12,2.16};
double Cond_Fr[]={0.81,1.12,1.52,2.54,2.51};
double Cvol[]={2.1e6,       2.4e6,  2.2e6,  1.8e6,  2.7e6};
double FIT=0.08;
double ALFA0 = 20.14;
int Dn = 230;
int N=230; //depth nodes //N=number of grid points
int Dn_out = 15; //Depth for output
int Dn_init = 14;

int m;

float S;
double L2;
double A1, B1;
double W1, W2;
double ALFA;
double H0,H1,H2;
double C0,C1;
double G1,G2,G3,G4;
double L0,L1;

// variables in the main function

double maxABS;
double E0 = 1e-6;
int M=numTimeSteps;
int iter;
int iter0=21;
double dinit, tinit;
double dayNo, junk;

//below extent of array default
int in_out_max = 15;
int maxim = 230;
int Time_max = 1500;


GIPL::GIPL () {
  snowDepth=NULL;
  snowDensity=NULL;
  airTemp=NULL;
  depth=NULL;
  T_init=NULL;    
  D_init=NULL; 
  Q1=NULL;    
  P1=NULL;    


//  T = new double*[Time_max];
  snowDensity = new double[Time_max];
  airTemp = new double[Time_max];
  depth = new double[maxim];
  snowDepth = new double[Time_max];
  T_init = new double[in_out_max];
  D_init = new double[in_out_max];
  Q1 = new double[maxim];
  P1 = new double[maxim];

}

GIPL::~GIPL () {
//  delete T;
  delete snowDensity;
  delete airTemp;
  delete snowDepth;
  delete depth;
  delete T_init;
  delete D_init;
  delete Q1;    
  delete P1;    

}



void GIPL::init(void)
{	
FILE *e_fptr;
  int j, i;
  std::string line;
  std::ifstream myfile;
	//---------- 1 --------------------------
    myfile.open ("snowdepth.txt");
	if (myfile.is_open())
	{
		for (i = 0; i < numTimeSteps; i++)
		{
		    getline	(myfile, line);
		    std::istringstream stm;
		    stm.str(line);
		    double tmpdayNo;
		    stm >> tmpdayNo;
			stm >> snowDepth[i];
		}
    }
    myfile.close();
	// --------------  2  -------------------
    myfile.open ("Snowdensity.txt");
    if (myfile.is_open()){
	    for (i = 0; i < numTimeSteps; i++){
		    getline	(myfile, line);
		    std::istringstream stm;
		    stm.str(line);
		    double tmpdayNo;
		    stm >> tmpdayNo;
		    stm >> snowDensity[i];}	
    }
    myfile.close();
	// ------------- 3 -----------------------
    myfile.open ("AirTemperature.txt");
	if (myfile.is_open())
	{
		for (i = 0; i < numTimeSteps; i++)
		{
		    getline	(myfile, line);
		    std::istringstream stm;
		    stm.str(line);
		    double tmpdayNo;
		    stm >> tmpdayNo;
			stm >> airTemp[i];
		}
    }
	myfile.close();
    //------------------  4  -----------------
    myfile.open("depthGrid.txt");
	if (myfile.is_open())
	{
		for (i = 0; i < N; i++)
		{
		    getline	(myfile, line);
		    std::istringstream stm;
		    stm.str(line);
			stm >> depth[i];
		}
    }
    myfile.close();
	// ---------------  5  --------------------
	e_fptr=fopen( "Initial.txt", "r");
	for(i=1;i<=14;i++)
        {
		    double junk;
		    fscanf(e_fptr, "%lf %lf \n", &dinit ,&tinit);
		    D_init[i]=dinit;    //depth
		    T_init[i]=tinit;    // initemper
		}
	fclose(e_fptr);

}  	// void GIPL::init(void)

//====================
void GIPL::run(void)
{
double G0[Time_max];
int i,j,k = 1;
  double *U1=NULL;    
  U1 = new double[maxim];
  double Timex[Time_max];
  double T[maxim][Time_max];
//  double* T = new double[Time_max];
  FILE *out_fptr;

    for (j=1;j<=numTimeSteps;j++){G0[j]=0.015;}

    cout.setf(ios::fixed | ios::showpoint);
    cout.precision(5); 	

//	for (i = 1; i <= 14; i++){//do i=1,Dn
//            cout << "inid = " << D_init[i] << " iniT = " <<  T_init[i] << endl;
//    	}

//   interpolation of the initial conditions
     for (i = 0; i < Dn; i++){	//do i=1,Dn
	     U1[i] = interp_1(depth[i], D_init, T_init, 14);   // call for linear interpolation
	     T[i][0]=U1[i];// gives a segmentation fault on that
//        cout << "depth = " << depth[i] << " U1 = " << U1[i] << endl;
     }
//    cout << " Thick = " << Thick[0] << "numSoilLayers = " << numSoilLayers <<endl;
     for (i = 1;i< numSoilLayers ;i++){
	     Thick[i]=Thick[i]+Thick[i-1];
//    		cout  << " Thick = " << Thick[i] << endl;
     }


     for (i=0;i< M ;i++){ //	--- Time grid
	 Timex[i]=i*TStep;
//	    cout << " i = " << i << "Timex = " << Timex[i] <<endl;
     }
    
    S=TStep*60*60;//   ! 24 hours time step in seconds
    
    for (j=1;j< 2;j++){	// to M DO j=2,M  ! begining of the time cycle
	iter = 0.0;
	maxABS = 1.41e-6;
	while (iter < iter0 && maxABS > E0){
	    //NOTE: Please, see the comments in soilThermalConductivity function
	    L0=soilThermalConductivity(depth[0],U1[0]);   
	    L1=soilThermalConductivity(depth[1],U1[1]);   
	    printf ("ITER = %i L0 %f L1 %f\n",iter, L0, L1,"\n");

	    H0=depth[1]-depth[0];
	    printf ("depth1 %f depth0 %f\n", depth[1], depth[0],"\n");
     	    printf ("snowDepth[j] %f depth0 %f\n", snowDepth[j], depth[0],"\n");
	    printf ("T0 %f airTemp[j] %f\n", T[0][j-1], airTemp[j],"\n");
	    printf ("depth[1] %f T[][] %f\n", depth[1], T[0][j-1],"\n");
  	    if (snowDepth[j] == 0.0){
 		G1=0.0;
		G2=airTemp[j];	    }
            //	else if (snowDepth[j] > 0.0)
	    else { //! if (snowDepth(j) > 0.) then call alfaSnow(j)! heat transfer within the snow cover
		ALFA =  alfaSnow(j, snowDensity[j], snowDepth[j]);
		ALFA=(double) 1/ALFA;       //alfa is another coeficient
		C1=heatCapacity(depth[0],T[0][j-1]);   //main solution para
		W1=(double) 0.5* (L0+L1);
		W2=H0*ALFA/W1;
		W1=(double) 0.5*pow(H0,2)*C1/W1/S;
		G1=(double)1.0 + W1 + W2; //G1 G2 are coefficient of upper boundary
		G2=(W2*airTemp[j]+W1*T[0][j-1])/G1;
		G1=(double) 1.0 / G1;
	    }//endif
  //	----- Permutation and forward elimination
	    P1[1]=G1;
	    Q1[1]=G2;
            for (i=1;i<N-1;i++){//do i=2,N-1 // depth and temperature main where temperature takes previous time step
		printf ("depth[i] %f T[][] %f\n", depth[i+1], T[i+1][j-1],"\n");
		C1=heatCapacity(depth[i],T[i][j-1]);
		L2=soilThermalConductivity(depth[i+1],T[i+1][j-1]);
		printf ("C1 %f L2 %f\n", C1, L2,"\n");
		H1=depth[i+1]-depth[i];
		H2=(double) 0.5*(H0+H1);
		A1=(double) 0.5 * (L0+L1)*S/C1/(H0*H2);
		B1=(double) 0.5 *(L1+L2)*S/C1/(H1*H2);
		C0=(double) 1.0+A1+B1;
		P1[i+1]=B1/(C0-A1*P1[i]);
		Q1[i+1]=(A1*Q1[i]+T[i][j-1])*P1[i+1]/B1;
		H0=H1 ;
		L0=L1 ;
		L1=L2;}
  //	  computation of the Lower boundary koef. G3 & G4
	    C1=heatCapacity(depth[N-1],T[N-1][j-1]);
	    G3=(double) 0.5*pow(H1,2)*C1/L2/S ;
	    G4=H1*G0[j]+G3*T[N-1][j-1];
	    G3=(double) 1.0/(1.0+G3);
	    G4=G4*G3;
  //	 Temperature computation in the last (deepest) grid node
	    W1=(G3*Q1[N-1]+G4)/(1.0-G3*P1[N-1]);
	    maxABS=fabs(W1-U1[N-1]);
	    U1[N-1]=W1;
  //---- Back substitution
	    while (i>=0){	//DO WHILE (I>=1)
		W1=P1[i+1]*U1[i+1]+Q1[i+1];
		if (fabs(W1-U1[i])>maxABS) maxABS=fabs(W1-U1[i]);  //check for the iterative convergence
		U1[i]=W1;
		i=i-1;}//ENDDO !WHILE
	    iter=iter+1;
	   }//enddo ! while ((ITER < ITER0).AND.(maxABS > E0))
	   for (i=0;i<N;i++){ //do i=1,N
		T[i][j]=U1[i]; //printf ("ITER = %i j%i TEMP %f\n",i,j, T[i][j]);
	   }//enddo
    } //ENDDO	! j=2,M time cycle

   out_fptr = fopen("Results.txt", "w");{
      for (i=0;i<numTimeSteps;i++){
	   for(j=0;j<=N;j++){
	      fprintf(out_fptr, "%d %lf\n", i, T[j][i]);}
      }
   }
   fclose(out_fptr);
    
}// end of the run function

double GIPL::interp_1(double x, double xi[], double yi[], int imax)
/*
====================================================================
 Linear interpolation  
 Alex Godunov (demo for Phys420)
--------------------------------------------------------------------
 input ...
 x   - the abscissa at which the linear interpolation is to be evaluated
 xi[]- the arrays of data abscissas
 yi[]- the arrays of data ordinates
 imax- size of the arrays xi[] and yi[]

 output ...
 y    - interpolated value
====================================================================
*/
{
    double y;
    int j;
//    printf ("x %f xi[0] %f yi[0] %f\n", x, xi[imax], yi[imax]);
// if x is ouside the xi[] interval take a boundary value (left or right)
    if (x <= xi[0])      return y = yi[1];
    if (x >= xi[imax]) return y = yi[imax];    
// loop to find j so that x[j-1] < x < x[j]
    j = 0;
    while (j <= imax)
    {
     if (xi[j] >= x) break;
     j = j + 1;
	}
    y = yi[j-1] + (yi[j] - yi[j-1])*(x - xi[j-1])/(xi[j]-xi[j-1]);

    return y;
}

//--------------------------------------------------------------------
double GIPL::alfaSnow (int j0, double rho_snow, double d_snow)
{
	double ALFA;
		ALFA=1.0/ALFA0 +1./((0.018+0.00087*rho_snow)/d_snow);
	return ALFA;
}
//---------------------------------------------------
double GIPL::heatCapacity (double depth0,double temper){// Volumetric Heat Capacity
   double unw, CAP;

   int i=numLayers(depth0); // soil layer number
   if (temper < Tfr[i]-FIT){
      unw=unfrWater(temper,aclv[i],bclv[i],cclv[i]);
      Wunf[i]=unw;
      CAP=Cvol[i]+5e2*(Wvol[i]-Wunf[i])+1e6*Wunf[i];}
   else if (temper > Tfr[i]+FIT){
      CAP=Cvol[i]+1e6*Wvol[i];}
   else{
      unw=unfrWater(temper,aclv[i],bclv[i],cclv[i]);
      Wunf[i]=unw;
      CAP=0.5*((Cvol[i]+1e6*Wvol[i])+(Cvol[i]+5e2*(Wvol[i]-Wunf[i])+1e6*Wunf[i]))+0.5*334.2*1e6*(Wvol[i]-Wunf[i] )/FIT;}
   return CAP;
}

//--------------------------------------------------------------------------------------------
double GIPL::soilThermalConductivity(double depth,double temper)
//      Soil Thermal Conductivity
   {
   int i;
   double ice,ice_max,unw,fr_t,COND;
   fr_t = -30.0;
	   i=numLayers(depth);
	   
	   if (temper < Tfr[i]-FIT) // frozen
	   {
	   //COND=Cond_Fr[i];
	   //ice=Wvol(i)-unw
	   unw=unfrWater(temper,aclv[i],bclv[i],cclv[i]);
	   ice=Wvol[i]-unw;
	   unw=unfrWater(fr_t,aclv[i],bclv[i],cclv[i]);
	   ice_max=Wvol[i] -unw;
	   COND=Cond_Th[i]+( Cond_Fr[i]-Cond_Th[i] )*(ice/ice_max);
	   }
	   else if (temper > FIT+Tfr[i]) // thawed
	   {
	   unw=unfrWater(fr_t,aclv[i],bclv[i],cclv[i]);
	   ice_max=Wvol[i]; //- unw
	   COND=Cond_Th[i]+( Cond_Th[i]-Cond_Th[i] )*(Wvol[i]/ice_max);
	   // COND=Cond_Th[i];
	   }
	   else
	   { // within the smoothing interval = very close to zero C
	   COND=(double) 0.5 *(Cond_Th[i]+Cond_Fr[i]);
	   }
	   return COND;
//    double Cond_Th[]= {0.34,0.73,1.23,2.12,2.16};
 //   double Cond_Fr[]={0.81,1.12,1.52,2.54,2.51};
   }

int GIPL::numLayers(double depth){
   int Numl,i=1;

   while (depth > Thick[i]+1e-3){ i=i+1;}
   Numl=i;
   return Numl;}

   
//----------------------------------------------------------------------------------------------------
double GIPL::unfrWater(double temper, double ac, double bc, double cc)
{
   double unWater, diff;
   diff = fabs(cc-temper);
   unWater=ac*(pow(diff,bc));

   return unWater;
}
//----------------------------------------------

