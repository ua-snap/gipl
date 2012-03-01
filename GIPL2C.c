//---------------------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#define TRUE 1
#define FALSE 0



int numTimeSteps=1500;
int m;
	int numSoilLayers=5;
	int TStep =24;
	// soil properties
	 double Thick[]= {0.5, 0.5 ,9 ,20, 70};
	 int Num[]={5, 5, 10, 19, 14};
	 double Tfr[]={0,0,0,0,0};
	 double Wvol[]={0.48,0.42,0.31,0.26,0.22};
	 double Wunf[]={0.02, 0.02, 0.01,0.02,0.02};
	 double aclv[]={0.005,	0.035,	0.061,	0.018,	0.064};
	 double bclv[]={-0.1,	-0.32,	-0.35,	-0.17,	-0.340};
	 double cclv[]={0.,		0,		0.0,	0.0,	0.};;
	 double Cond_Th[]= {0.34,0.73,1.23,2.12,2.16};
	 double Cond_Fr[]={0.81,1.12,1.52,2.54,2.51};
   double Cvol[]={2.1e6,	2.4e6,	2.2e6,	1.8e6,	2.7e6};
   double FIT=0.08;
   double ALFA0 = 20.14;
		int Dn = 230;
		int N=230; //depth nodes
		int Dn_out = 15; //Depth for output
  	int Dn_init = 14;

	float S;
	double L2;
	double A1, B1;
	double W1, W2;
 	double ALFA;
 	double H0,H1,H2;
	double C0,C1;
  double G1,G2,G3,G4;
	double L0,L1;

//for memory

//for subroutines
  extern int numLayers(double depth);
	extern double soilThermalConductivity(double depth,double temper);
  extern double snowProperties (int j0, double *snowDensity, double *snowDepth);
	extern double heatCapacity (double depth0,double temper);
	extern double unfrWater(double temp, double ac, double bc, double cc);
	
int main(int argc, char* argv[])
{
FILE *a_fptr;
FILE *b_fptr;
FILE *c_fptr;
FILE *d_fptr;
FILE *e_fptr;
FILE *f_fptr;
FILE *out_fptr;

  int j;
	int i, i0, k;
	k = 1;
  double maxABS;
	double E0 = 1e-6;
	int M=numTimeSteps;
	int iter;
	int iter0=21;
  double dinit, tinit;
  double dayNo, junk;

//all dimensional parameters starts
 	double *T_init=NULL;    
  double *D_init=NULL; 
  double *D_out=NULL;  
	double *U1=NULL;    
	double *Q1=NULL;    
	double *P1=NULL;    
  //double **T=NULL;
  double *snowDepth=NULL;
  double *snowDensity=NULL;
  double *airTemp=NULL;
  double *depth=NULL;
  //double *G0=NULL;
  //double *Timex=NULL;

	
	//below extent of array default
  int Time_max; 
  int in_out_max; 
  int maxim;
  //default values below
  in_out_max = 15;
  maxim = 230;
  Time_max = 1500;

  //allocate memory
T_init = new double[in_out_max];
D_init = new double[in_out_max];
D_out = new double[in_out_max];
U1 = new double[maxim];
Q1 = new double[maxim];
P1 = new double[maxim];
snowDepth = new double[Time_max];
snowDensity = new double[Time_max];
airTemp = new double[Time_max];
depth = new double[maxim];


double **T;
T = new double*[Time_max];
for (int i = 0; i < Time_max; i++){
	T[i] = new double[Time_max];
}


double G0[Time_max];
double Timex[Time_max];
std::string line;
	std::ifstream myfile;
	myfile.open ("snowdepth.txt");
	if (myfile.is_open()){
		for (i = 0; i < numTimeSteps; i++){
			getline	(myfile, line);
			std::istringstream stm;
			stm.str(line);
			double tmpdayNo;
			stm >> tmpdayNo;
			stm >> snowDepth[i];
		}	
	}
	myfile.close();
 // --------------  2  ----------------------
	myfile.open ("Snowdensity.txt");
	if (myfile.is_open()){
		for (i = 0; i < numTimeSteps; i++){
			getline	(myfile, line);
			std::istringstream stm;
			stm.str(line);
			double tmpdayNo;
			stm >> tmpdayNo;
			stm >> snowDensity[i];
		}	
	}
	myfile.close();

//-----------------  3 ----------
	myfile.open ("AirTemperature.txt");
	if (myfile.is_open()){
		for (i = 0; i < numTimeSteps; i++){
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
	if (myfile.is_open()){
		for (i = 0; i < N; i++){
			getline	(myfile, line);
			std::istringstream stm;
			stm.str(line);
			stm >> depth[i];
		}	
	}
	myfile.close();


 //-----------------  5  --------------------------

		e_fptr=fopen( "Initial.txt", "r");
	
	for(i=1;i<=14;i++) 
            {
			double junk;
			fscanf(e_fptr, "%lf %lf \n", &dinit ,&tinit);
			D_init[i]=dinit;    //depth
			T_init[i]=tinit;    // initemper
						 }
		fclose(e_fptr);


 //----------------- 6 -----------------------------

		f_fptr=fopen( "Initial.txt", "r");
	
	for(i=1;i<=14;i++) 
            {
			double junk;
			fscanf(f_fptr, "%lf %lf \n", &dinit ,&tinit);
			D_init[i]=dinit;    //depth
			T_init[i]=tinit;    // initemper
						 }
		fclose(f_fptr);


//--------------------------------------------------


    for (j=1;j<=numTimeSteps;j++){G0[j]=0.015;}

//---------------------------------------------
  	//!  interpolation of the initial conditions
	for (i = 0; i < Dn; i++){//do i=1,Dn
		double dividend = T_init[k]-T_init[k-1];
		double divisor = (D_init[k]-D_init[k-1]) * (depth[i]-D_init[k-1]);
		//T[i][0]=T_init[k-1]+ ((T_init[k]-T_init[k-1]) / (D_init[k]-D_init[k-1])  *  (depth[i]-D_init[k-1])) ;
		if (divisor != 0){
			T[i][0] = T_init[k-1] + dividend / divisor;
		} else {
			T[i][0] = 0;
		}
		  U1[i]=T[i][0];
		if (depth[i]==D_init[k]) k=k+1;
		printf ("i = %i k = %i Init T %f Init D %f\n", i,k,T_init[i], D_init[i]);
	}//enddo


 	//for (j=0;j< 14;j++){printf ("Depth %f Temp %f\n",D_init[j], T_init[j]);}
    
    for (i0 = 1;i0< numSoilLayers ;i0++){
		Thick[i0]=Thick[i0]+Thick[i0-1];
    }// enddo
//	!--- Time grid
		for (i0=0;i0< M ;i0++){
			Timex[i0]=i0*TStep;
    }//	enddo


//	!--------------
	 S=TStep*60*60;//   ! 24 hours time step in seconds


	 for (j=1;j< M+1;j++){// to M //DO j=2,M	! begining of the time cycle

			iter = 0;
			maxABS = 1.41e-6;
			while (iter < iter0 && maxABS > E0){
				//printf ("boundary coeff\n" );
//			! computation of boundary coefficients  G1,G2
				L0=soilThermalConductivity(depth[0],U1[0]);
				L1=soilThermalConductivity(depth[1],U1[1]);   ///main solution para/lamda 

				H0=depth[1]-depth[0];
				if (snowDepth[j] < E0 || T[0][j-1] > E0){
					G1=0.0;
					G2=airTemp[j];
				}else if (T[0][j-1] <= 0.0 &&  snowDepth[j] < E0) {
					     G1=0.0;
						 G2=airTemp[j];
				}

				else{ //! if (snowDepth(j) > 0.) then
//				 call snowProp(j)	! heat transfer within the snow cover
						ALFA = 	snowProperties(j, snowDensity, snowDepth);
				    ALFA=(double) 1/ALFA;       //alfa is another coeficient
					  C1=heatCapacity(depth[0],T[0][j-1]);   //main solution para
					  W1=(double) 0.5* (L0+L1);              
					  W2=H0*ALFA/W1;
					  W1=(double) 0.5*pow(H0,2)*C1/W1/S;
					  G1=(double)1.0 +W1+W2;                   //G1 G2 are coefficient of boundary
					  G2=(W2*airTemp[j]+W1*T[0][j-1])/G1;
					  G1=(double) 1/G1;
				}//endif

//				!----- Permutation and forward elimination
				P1[1]=G1;
				Q1[1]=G2;
			//printf ("N = = %i\n",N);
		 		for (i=1;i<N-1;i++){//do i=2,N-1

				  C1=heatCapacity(depth[i],T[i][j-1]);     // depth and temperature main where temperature takes previous time step
				  L2=soilThermalConductivity(depth[i+1],T[i+1][j-1]);
	               H1=depth[i+1]-depth[i];
				         H2=(double) 0.5*(H0+H1);
	               A1=(double) 0.5 * (L0+L1)*S/C1/(H0*H2);
	               B1=(double) 0.5 *(L1+L2)*S/C1/(H1*H2);
	               C0=(double) 1.0+A1+B1;
	               P1[i+1]=B1/(C0-A1*P1[i]);
	               Q1[i+1]=(A1*Q1[i]+T[i][j-1])*P1[i+1]/B1;
	               H0=H1 ;
				   L0=L1 ;
				   L1=L2;

				}
//				! computation of the Lower boundary koef. G3 & G4
				C1=heatCapacity(depth[N-1],T[N-1][j-1]);
	            G3=(double) 0.5*pow(H1,2)*C1/L2/S ;
	            G4=H1*G0[j]+G3*T[N-1][j-1];

	            G3=(double) 1.0/(1.0+G3);
	            G4=G4*G3;
//					! Temperature computation in the last (deepest) grid node
					W1=(G3*Q1[N-1]+G4)/(1.0-G3*P1[N-1]);
					maxABS=fabs(W1-U1[N-1]);
					U1[N-1]=W1;
//				!---- Back substitution
				i=N-2;
				while (i>=0){//DO WHILE (I>=1)
					W1=P1[i+1]*U1[i+1]+Q1[i+1];
//					! check for the iterative convergence
					if (fabs(W1-U1[i])>maxABS) maxABS=fabs(W1-U1[i]);
						U1[i]=W1;
						i=i-1;
				}//ENDDO !WHILE
				iter=iter+1;
			}//enddo ! while ((ITER < ITER0).AND.(maxABS > E0))
			//			write(201,*) iter
			for (i=0;i<N;i++){//do i=1,N

			   T[i][j]=U1[i];
		//	   printf ("ITER = %i j%i TEMP %f\n",i,j, T[i][j]);
			}//enddo
//
	 } //ENDDO	! j=2,M
//		close(201)
//	return
	
	out_fptr = fopen("Results.txt", "w");
  {
				for (i=0;i<numTimeSteps;i++)
      {
        for(j=0;j<=N;j++)  //N=number of grid points
        {
        fprintf(out_fptr, "%d %lf\n", i, T[j][i]);
        }
      }
				
	}

	fclose(out_fptr);
	
	return 0;
}

//---------------------------------------------------------------------------
	int numLayers(double depth){
	// Number of the Soil Layer Function
	//use GLOBVAR
	//implicit none
	//	real*8 depth
	int Numl;
		int i;
	      i=1;
	      while (depth > Thick[i]+1e-3){
			i=i+1;
	      }
	      Numl=i;

		return Numl;
		 }
	double soilThermalConductivity(double depth,double temper)
//	!/* Soil Thermal Conductivity
	{
	double COND;

	int i;
	        i=numLayers(depth);
	        if (temper < Tfr[i]-FIT){
					COND=Cond_Fr[i];
	        }
			else if (temper > FIT+Tfr[i]){
			    COND=Cond_Th[i];
	        }else{
			    COND=(double) 0.5 *(Cond_Th[i]+Cond_Fr[i]);
			}
			return COND;
	}

double snowProperties (int j0, double *snowDensity, double *snowDepth)
{

	double ALFA;
		ALFA=1.0/ALFA0 +1./((0.018+0.00087*snowDensity[j0])/snowDepth[j0]);
	return ALFA;
}
	double heatCapacity (double depth0,double temper)
	{
	// Volumetric Heat Capacity
	double CAP;
	        int i=numLayers(depth0);
			if (temper < Tfr[i]-FIT){
				CAP=Cvol[i]+5e2*(Wvol[i]-Wunf[i])+1e6*Wunf[i];
	        }
	        else if (temper > Tfr[i]+FIT){
				CAP=Cvol[i]+1e6*Wvol[i];
	        }
	        else{
				CAP=0.5*( (Cvol[i]+1e6*Wvol[i] )+( Cvol[i]+5e2*
				( Wvol[i]-Wunf[i] )+1e6*Wunf[i] ) )+0.5*334.2*1e6*
				 (Wvol[i]-Wunf[i] )/FIT;
			}
			return CAP;
	   }

	double unfrWater(double temper, double ac, double bc, double cc)
	{
		double unWater;
    double diff;
    diff = fabs(cc-temper);
	 unWater=ac*(pow(diff,bc));

return unWater;
	}

