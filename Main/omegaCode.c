
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define nsp 7
#define nr 3
#define nasaC 7
#define speChara 8
#define Ru 8.134
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Arrhenius coefficient
void  k_coeff(double* k, double* As, float* betas, double* Tas, double T);

// Returning the adequate NASA coefficient: // Ignore the case where the temp is lower or higher than the range
void nasa_coeff(double nasaCoeff[][nasaC], const float nasaTemp[][3], const double low[][nasaC], const double up[][nasaC],  double T);

// Calculating the equilibrium constant for each reaction:
void Kceq(double* Kc, int vnet[][nr], double nasaCoeff[][nasaC], double T);

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main() {

  //-------------------------------- Species Details ----------------------------------//
  // The array of species' names:
  char names[nsp][speChara] =  {"N", "NO", "N2", "O", "O2", "H", "OH"}; 

  // The vector of molecular weight:
  char wk[nsp] =  { 0.001, 0.0002, 0.02 }; 

  // Array of the NASA temperature range:
  const float Tnasa[nsp][3] =  { 200.0, 1000.0, 6000.0 }; 

  // Array of the NASA  polynomials lower range coefficients:
  const double nasaLower[nsp][nasaC] =  { 200.0, 1000.0, 6000.0, 2.0, 6.0, 8.0, 7.0 }; 

  // Array of the NASA polynomials upper range coefficients:
  const double nasaUpper[nsp][nasaC] =  { 200.0, 1000.0, 6000.0, 2.0, 6.0, 8.0, 7.0 }; 

  // Vector of enthalpy of formation:
  const double Hf[nsp] =  { 200.0, 1000.0, 6000.0, 2.0, 6.0, 8.0, 7.0 }; 

  
  //------------------------------ Mechanism Details ---------------------------------//
  //const int nsp=7; // Number of species
  //const int nr=6; // 3 reversible reactions = 6 irreversible reactions.

  // The matrix of reactants coefficients:
  int vreact[nsp][nr] =  {
    {vr, vr},
    {vr, vr} 
   };

  // The matrix of products coefficients:
  int vprod[nsp][nr] =  vp;
  
  // The array of the pre-exponent factor:
  double A[nr] =  Arr;
  
  // The array of the powers beta:
  float beta[nr] =  beta;

  // The array of activation temperature:
  double Ta[nr] =  Ta;
  
  //---------------------------------------------------------------------------------//

  //---------------------------- Initial conditions ---------------------------------//
  double T =  1900.0;
  double P =  1E+05;

  double Yk[nsp] = @Yk

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

  double nasa7[nsp][nasaC];
  nasa_coeff(nasa7, &Tnasa[0], &nasaLower[0], &nasaUpper[0], T);

  double KcEq[nr];
  Kceq(KcEq, &vnet[0], &nasa7[0], T);
  
  // Calculating kb,r:
  double kback[nr];
  for (int r=0; r<nr; ++r)
    kback[r] = kcoeff[r]/KcEq[r];

  //---------------------------- Production rate ------------------------------------//
  double omega[nsp] = {0.0};
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
 
    RR[r] = kcoeff[r]*dummyForward - kback[r]*dummyBackward;
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
  for (int k=0; k<nsp; ++k) printf("Production rate of %8s is %g \n", names[k], omega[k]);
  

  
  return 0;


}

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

    Kc[r] = pow(Ru*T, -vsum)* exp(as[0]*(log(T)-1) + as[1]*T/2 + as[2]*T*T/6 + as[3]*pow(T,3)/12 + as[4]*pow(T, 4)/20 - as[5]/T + as[6]);
  }
}
