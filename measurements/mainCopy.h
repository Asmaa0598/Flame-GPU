#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define nsp 7
#define nr 3
#define nasaC 9
#define Ru 8.134
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Arrhenius coefficient
void  k_coeff(double* k, double* As, float* betas, double* Tas, double T) {

  for (int r=0; r<nr; ++r) k[r] = As[r]*pow(T,betas[r])*exp(-Tas[r]/T);
  //for (int r=0; r<nr; ++r) printf("%g\n", k[r]);
}

// Returning the adequate NASA coefficient: // Ignore the case where the temp is lower or higher than the range
void nasa_coeff(double nasaCoeff[][nasaC], const float nasaTemp[][3], const double low[][nasaC], const double up[][nasaC],  double T) {

  for (int k=0; k<nsp; ++k)
    if ( T >= nasaTemp[k][0] && T <= nasaTemp[k][1] ) {
      for (int i=0; i<nasaC; ++i)
	nasaCoeff[k][i] = low[k][i];
    } else if (T >= nasaTemp[k][1] && T <= nasaTemp[k][2]) {
      for (int i=0; i<nasaC; ++i)
	nasaCoeff[k][i] = up[k][i];
    } else {
      printf("ERROR!! OUT OF RANGE");
      abort();
    }
}

// Calculating the equilibrium constant for each reaction:
void Kceq(double* Kc, int vnet[][nr], double nasaCoeff[][nasaC], double T) {

  for (int r=0; r<nr; ++r) {
    int vsum = 0.0;
    double as[nasaC] = {0.0};
    
    for (int i=0; i<nasaC; ++i) {
      for(int k=0; k<nsp; ++k){
	vsum += vnet[k][r];
	as[i] += vnet[k][r] * nasaCoeff[k][i];
      }
    }

    // NASA 7:
    //Kc[r] = pow(Ru*T, -vsum)* exp(as[0]*(log(T)-1) + as[1]*T/2 + as[2]*T*T/6 + as[3]*pow(T,3)/12 + as[4]*pow(T, 4)/20 - as[5]/T + as[6]);

    // NASA 9:
    
    Kc[r] = exp(as[8]);
    Kc[r] *=  exp(0.5*as[0]*pow(T, -2) - as[1]*(log(T)+1.0)/T - as[2]*(1.0 - log(T)) + 0.5*as[3]*T + as[4]*T*T/6.0 + as[5]*T*T*T/12.0 + as[6]*pow(T,4)/20.0 - as[7]/T );
    Kc[r] *= pow(2.3025851*Ru*T, -vsum);
    
    /*
    double gibbs = -as[8] - 0.5*as[0]/T/T + as[1]*(log(T)+1.0)/T + as[2]*(1.0 - log(T)) - 0.5*as[3]*T - as[4]*T*T/6.0 - as[5]*T*T*T/12.0 - as[6]*pow(T,4)/20.0 + as[7]/T;
    gibbs -= vsum*log(1.0E+05/Ru/T);
    Kc[r] = exp(-gibbs);
    */
  }
}




// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void zeldovich(double* omega, double Ts, double Ps, double Yk[nsp]) {

  //-------------------------------- Species Details ----------------------------------//
  // The array of species' names:
  char names[nsp][3] = {"N", "NO", "N2", "O", "O2", "H", "OH"};

  // The array of species Molecular weight:
  const double wk[nsp] =  { 0.014007, 0.030006, 0.028013, 0.015999, 0.031999, 0.001008, 0.017007 };

  // Array of the NASA temperature range:
  const float Tnasa[nsp][3] = {
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000700, 6000.000700 }, 
    { 200.000000, 1000.000000, 6000.000000 } 
 };

  // Array of the NASA  polynomials lower range coefficients:
  const double nasaLower[nsp][nasaC] = { 
 { 0.000000, 0.000000, 2.500000, 0.000000, 0.000000, 0.000000, 0.000000, 56104.637800, 4.193909 }, 
 { -11439.165030, 153.646759, 3.431469, -0.002669, 0.000008, -0.000000, 0.000000, 9098.214410, 6.728727 }, 
 { 22103.714970, -381.846182, 6.082738, -0.008531, 0.000014, -0.000000, 0.000000, 710.846086, -10.760033 }, 
 { -7953.611300, 160.717779, 1.966226, 0.001014, -0.000001, 0.000000, -0.000000, 28403.624370, 8.404242 }, 
 { -34255.634200, 484.700097, 1.119011, 0.004294, -0.000001, -0.000000, 0.000000, -3391.454870, 18.496995 }, 
 { 0.000000, 0.000000, 2.500000, 0.000000, 0.000000, 0.000000, 0.000000, 25473.708010, -0.446683 }, 
 { -1998.858990, 93.001362, 3.050854, 0.001530, -0.000003, 0.000000, -0.000000, 3239.683480, 4.674111 } 
 }; 

  const double nasaUpper[nsp][nasaC] = { 
 { 88765.0138, -107.12315, 2.362188287, 0.0002916720081, -1.7295151e-07, 4.01265788e-11, -2.677227571e-15, 56973.5133, 4.86523579 }, 
 { 223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.88096645e-08, -1.416076856e-11, 9.38018462e-16, 17503.17656, -8.50166709 }, 
 { 587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86639599 }, 
 { 261902.0262, -729.872203, 3.31717727, -0.000428133436, 1.036104594e-07, -9.43830433e-12, 2.725038297e-16, 33924.2806, -0.667958535 }, 
 { -1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16, -16890.10929, 17.38716506 }, 
 { 60.7877425, -0.1819354417, 2.500211817, -1.226512864e-07, 3.73287633e-11, -5.68774456e-15, 3.410210197e-19, 25474.86398, -0.448191777 }, 
 { 1017393.379, -2509.957276, 5.11654786, 0.000130529993, -8.28432226e-08, 2.006475941e-11, -1.556993656e-15, 20444.8713, 0 } 
 };
  
  //------------------------------ Mechanism Details ---------------------------------//
  //const int nsp=7; // Number of species
  //const int nr=6; // 3 reversible reactions = 6 irreversible reactions.

  // The matrix of reactants coefficients:
  int vreact[nsp][nr] = {
    { 1, 1, 1 }, //N
    { 1, 0, 0 }, //NO
    { 0, 0, 0 }, //N2
    { 0, 0, 0 }, //O
    { 0, 1, 0 }, //O2
    { 0, 0, 0 }, //H
    { 0, 0, 1 } }; //OH

  // The matrix of products coefficients:
  int vprod[nsp][nr] = {
    { 0, 0, 0 },
    { 0, 1, 1 },
    { 1, 0, 0 },
    { 1, 1, 0 },
    { 0, 0, 0 },
    { 0, 0, 1 },
    { 0, 0, 0 } };
  
  // The array of the pre-exponent factor:
  double A[nr] = { 1.8E+14, 1.8E+10, 7.1E+13 };
  
  // The array of the powers beta:
  float beta[nr] = {0.0, 0.0, 0.0};

  // The array of activation temperature:
  double Ta[nr] = { 38370.0, 4680.0, 450.0 };

  // Array of third body efficiencies:
  double M[nr] = {1.0, 1.0 , 1.0 };

  // Reversibility:
  float reversibility[nr] = {1.0};
  
  
  //---------------------------------------------------------------------------------//

  //---------------------------- Initial conditions ---------------------------------//
  double T = Ts;
  double P = Ps;

  //double Yk[nsp] = { 2.31141, 2.31141, 0.0, 2.38145, 0.0, 0.0, 0.0 };
  //double Yk[nsp] =;
  //double Yk[nsp] = { 0.5, 0.3, 0.0, 0.2, 0.0, 0.0, 0.0 };
  double R = 0.0;
  for (int k=0; k<nsp; ++k) R += Yk[k]/wk[k];
  R *= Ru;
  double rho = P/R/T;
  //---------------------------------------------------------------------------------//

  //------------------------------ Arrhenuis Law ------------------------------------//

  // The matrix of net stoichiometric coefficients:
  int vnet[nsp][nr];
  for (int k=0; k<nsp; ++k) {
    for (int r=0; r<nr; ++r) {
      vnet[k][r] = vprod[k][r] - vreact[k][r];
      //printf("%d ", vnet[k][r]);
    }
    //printf("\n");
  }

  // Calculating concentrations:
  // [Xk] = \rho*Yk/Wk
  double concentrations_IC[nsp];
  for (int k=0; k<nsp; ++k) concentrations_IC[k] = Yk[k]*rho/wk[k];

  double kcoeff[nr];
  
  k_coeff(kcoeff, &A[0], &beta[0], &Ta[0], T);

  double nasa9[nsp][nasaC];
  nasa_coeff(nasa9, &Tnasa[0], &nasaLower[0], &nasaUpper[0], T);

  double KcEq[nr];
  Kceq(KcEq, &vnet[0], &nasa9[0], T);

  // Displaying Kc:
  //printf("\n\nEquilibrium Constant Kc:\n");
  //for (int r=0; r<nr; ++r)
  //printf("Kc,%d is %g \n", r, KcEq[r]);
  
  // Calculating kb,r:
  double kback[nr];
  //printf("\n Backward coefficients: \n");
  for (int r=0; r<nr; ++r) {
    kback[r] = kcoeff[r]/KcEq[r];
    //printf("kb,%d: %g \n", r, kback[r]); 
  }

  //---------------------------- Production rate ------------------------------------//
  //double omega[nsp] = {0.0};
  double tmp[nsp] = {0.0};

  double RR[nr] = {1.0};
  
  for (int r=0; r<nr; ++r) {
   double dummyForward = 1.0;
   double dummyBackward = 1.0;
   //printf("\nFor r=%d: \n", r);
    for (int k=0; k<nsp; ++k){
      //printf("[%3s]=%g, with v%d%d=%d  gives:%g  \n", names[k], concentrations_IC[k], k, r, vreact[k][r], pow(concentrations_IC[k], vreact[k][r]));
      dummyForward *= pow(concentrations_IC[k], vreact[k][r]);
      dummyBackward *= pow(concentrations_IC[k], vprod[k][r]);
    }
 
    RR[r] = kcoeff[r]*dummyForward - reversibility[r]*kback[r]*dummyBackward;
    RR[r] *= M[r]; 
    //printf("\nReaction rate of %d is %g and kf=%g and kb=%g \n", r, RR[r], kcoeff[r], kback[r]);
  }

   for (int k=0; k<nsp; ++k) {
    for (int r=0; r<nr; ++r) {
      tmp[k] += vnet[k][r]*RR[r];
    }
    omega[k] = wk[k]*tmp[k];
  }
  //--------------------------------------------------------------------------------//

  // Printing results:
  //for (int k=0; k<nsp; ++k) printf("Production rate of %3s is %g \n", names[k], omega[k]);


}



/*
  //---------------------------- Production rate ------------------------------------//

  // Displaying NASA coefficients:
  printf("Coefficients:");
  for (int k=0; k<nsp; ++k){
    printf("\n\n%3s:\n ", names[k]);
    for (int i=0; i<7; ++i)
      printf("%g   ", nasa7[k][i]);
  }

 


 double omega[nsp] = {0.0};
  double tmp[nsp] = {0.0};

  double RR[nr] = {1.0};
  for (int r=0; r<nr; ++r) {
    RR[r] = 1.0;
    printf("\nFor r=%d: \n", r);
    for (int k=0; k<nsp; ++k){
      printf("[%3s]=%g, with v%d%d=%d  gives:%g  \n", names[k], concentrations_IC[k], k, r, vreact[k][r], pow(concentrations_IC[k], vreact[k][r]));
      RR[r] *= pow(concentrations_IC[k], vreact[k][r]);
    }
    if (r ==3 ) printf("before: %g \n", RR[r]);
    RR[r] *= kcoeff[r];
    printf("\nReaction rate of %d is %g and k=%g \n", r, RR[r], kcoeff[r]);
  }
  
  for (int k=0; k<nsp; ++k) {
    for (int r=0; r<nr; ++r) {
  universal gas constant    tmp[k] += vnet[k][r]*RR[r];
    }
    omega[k] = wk[k]*tmp[k];
  }
  //--------------------------------------------------------------------------------//

  // Printing results:
   for (int k=0; k<nsp; ++k) printf("Production rate of %3s is %g \n", names[k], omega[k]);
  

 */
