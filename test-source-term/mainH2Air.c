// Here is the main code:

/*

The test case is the Zeldovich mechanism provided in the paper:
"Computational analysis of the Extended Zeldovich Mechanism" by Lucky Anetor, Christopher Odetunde, Edward E. Osakue

nb of Species: 13
Species: H2 H OH H2O HO2 H2O2 O O2 N N2 NO NO2 HNO2

nb of reactions


====> To find kb,r we use the formula:  kb,r = kf,r/ Kc,r
====> Kc,r = (Ru*T)^(- sum^nsp v_sr) exp(Delta_a1(lnT-1) + Delta_a2*T/2 + Delta_a3*T^2/6 + Delta_a4*T^3/12 +  Delta_a5*T^4/20 + Delta_a6/T + Delta_a7)
====> Delta_a = sum_(s=1)^nsp v_sr a_s

IC: [N]=[NO]=[O]=1 mol/m^3, P=1.106e05 Pa, and T=1900K.


[Xk] = \rho*Yk/Wk ==> P = \sum [Xk] Ru T

To compile:
gcc main.c -o main -lm
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

#define nsp 13
#define nr 25
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

int main() {
   //---------------------------- Initial conditions ---------------------------------//
  double T = 1900.0;
  double P = 1.1065E+05;

  //double Yk[nsp] = { 2.31141, 2.31141, 0.0, 2.38145, 0.0, 0.0, 0.0 };
  double Yk[nsp] = { 0.33, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.33, 0.0, 0.34, 0.0, 0.0, 0.0 };
  //double Yk[nsp] = { 0.5, 0.3, 0.0, 0.2, 0.0, 0.0, 0.0 };
  //---------------------------------------------------------------------------------//


  //-------------------------------- Species Details ----------------------------------//
  // The array of species' names:
  char names[nsp][5] = {"H2","H","OH","H2O","HO2","H2O2","O","O2","N","N2","NO","NO2","HNO2"};

  // The array of species Molecular weight:
  const double wk[nsp] =   { 0.002016, 0.001008, 0.017007, 0.018015, 0.033007, 0.034015, 0.015999, 0.031999, 0.014007, 0.028013, 0.030006, 0.046006, 0.047013 }; 

  // Array of the NASA temperature range:
  const float Tnasa[nsp][3] = { 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000000, 6000.000000 }, 
  { 200.000000, 1000.000000, 6000.000000 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000000, 6000.000000 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 }, 
  { 200.000000, 1000.000700, 6000.000700 } 
 }; 


  // Array of the NASA  polynomials lower range coefficients:
  const double nasaLower[nsp][nasaC] ={ 
 { 40783.2281, -800.918545, 8.21470167, -0.0126971436, 1.75360493e-05, -1.20286016e-08, 3.36809316e-12, 2682.48438, -2.20276405 }, 
 { 0, 0, 2.5, 0, 0, 0, 0, 25473.70801, -0.448191777 }, 
 { -1998.85899, 93.0013616, 3.050854229, 0.001529529288, -3.157890998e-06, 3.31544618e-09, -1.138762683e-12, 3239.68348, -11.01282337 }, 
 { -39479.6083, 575.573102, 0.931782653, 0.00722271286, -7.34255737e-06, 4.95504349e-09, -1.336933246e-12, -33039.7431, -7.97814851 }, 
 { -75988.8254, 1329.383918, -4.67738824, 0.02508308202, -3.006551588e-05, 1.895600056e-08, -4.82856739e-12, -5809.36643, 40.6685092 }, 
 { -92795.3358, 1564.748385, -5.97646014, 0.0327074452, -3.93219326e-05, 2.509255235e-08, -6.46504529e-12, -24940.04728, -46.5085566 }, 
 { -7953.6113, 160.7177787, 1.966226438, 0.00101367031, -1.110415423e-06, 6.5175075e-10, -1.584779251e-13, 28403.62437, -0.667958535 }, 
 { -34255.6342, 484.700097, 1.119010961, 0.00429388924, -6.83630052e-07, -2.0233727e-09, 1.039040018e-12, -3391.45487, 17.38716506 }, 
 { 0, 0, 2.5, 0, 0, 0, 0, 56104.6378, 4.86523579 }, 
 { 22103.71497, -381.846182, 6.08273836, -0.00853091441, 1.384646189e-05, -9.62579362e-09, 2.519705809e-12, 710.846086, -15.86639599 }, 
 { -11439.16503, 153.6467592, 3.43146873, -0.002668592368, 8.48139912e-06, -7.68511105e-09, 2.386797655e-12, 9098.21441, -8.50166709 }, 
 { -56420.3878, 963.308572, -2.434510974, 0.01927760886, -1.874559328e-05, 9.14549773e-09, -1.777647635e-12, -1547.925037, -43.0512991 }, 
 { 8591.98506, 120.3644046, 0.941297912, 0.01942891839, -2.253174194e-05, 1.384587594e-08, -3.47355046e-12, -11063.37202, 20.73967459 } 
 }; 


  const double nasaUpper[nsp][nasaC] = { 
 { 560812.338, -837.149134, 2.97536304, 0.00125224993, -3.74071842e-07, 5.93662825e-11, -3.60699573e-15, 5339.81585, -2.20276405 }, 
 { 60.7877425, -0.1819354417, 2.500211817, -1.226512864e-07, 3.73287633e-11, -5.68774456e-15, 3.410210197e-19, 25474.86398, -0.448191777 }, 
 { 1017393.379, -2509.957276, 5.11654786, 0.000130529993, -8.28432226e-08, 2.006475941e-11, -1.556993656e-15, 20444.8713, -11.01282337 }, 
 { 1034972.096, -2412.698562, 4.64611078, 0.002291998307, -6.83683048e-07, 9.42646893e-11, -4.82238053e-15, -13842.86509, -7.97814851 }, 
 { -1810669.724, 4963.19203, -1.039498992, 0.00456014853, -1.061859447e-06, 1.144567878e-10, -4.76306416e-15, -31944.1874, 40.6685092 }, 
 { 1489428.027, -5170.82178, 11.2820497, -8.04239779e-05, -1.818383769e-08, 6.94726559e-12, -4.8278319e-16, 14182.51038, -46.5085566 }, 
 { 261902.0262, -729.872203, 3.31717727, -0.000428133436, 1.036104594e-07, -9.43830433e-12, 2.725038297e-16, 33924.2806, -0.667958535 }, 
 { -1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16, -16890.10929, 17.38716506 }, 
 { 88765.0138, -107.12315, 2.362188287, 0.0002916720081, -1.7295151e-07, 4.01265788e-11, -2.677227571e-15, 56973.5133, 4.86523579 }, 
 { 587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86639599 }, 
 { 223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.88096645e-08, -1.416076856e-11, 9.38018462e-16, 17503.17656, -8.50166709 }, 
 { 721300.157, -3832.6152, 11.13963285, -0.002238062246, 6.54772343e-07, -7.6113359e-11, 3.32836105e-15, 25024.97403, -43.0512991 }, 
 { 878790.413, -3990.45503, 11.87349269, -0.000488190061, 7.13363679e-08, -5.37630334e-12, 1.581778986e-16, 12463.43241, 0 } 
 }; 


  // Vector of enthalpy of formation:
  const double Hf[nsp] = {0.56850013E+05, 1.09770882E+04, 0.00000000E+00, 2.99687009E+04, 0.00000000E+00, 0.26219035E+05, -1.74702052E+04};

  
  //------------------------------ Mechanism Details ---------------------------------//
  //const int nsp=7; // Number of species
  //const int nr=6; // 3 reversible reactions = 6 irreversible reactions.

  // The matrix of reactants coefficients:
  int vreact[nsp][nr] = {
    { 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, //H2
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, //H
    { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1 }, //OH
    { 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 }, //H2O
    { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //HO2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //H2O2
    { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0 }, //O
    { 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, //O2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //N
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, //N2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0 }, //NO
    { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1 }, //NO2
    { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }  //HNO2
  }; 

  // The matrix of products coefficients:
  int vprod[nsp][nr] = {
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //H2
    { 1, 2, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0 }, //H
    { 1, 0, 0, 0, 0, 1, 0, 2, 1, 1, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, //OH
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //H2O
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1 }, //HO2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //H2O2
    { 0, 0, 2, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, //O
    { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, //O2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0 }, //N
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, //N2
    { 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1 }, //NO
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0 }, //NO2
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }  //HNO2
  }; 
  
  // The array of the pre-exponent factor:
  double A[nr] = { 5.2E+15, 5.5E+12, 7.2E+12, 8.5E+12, 1.7E+10, 5.0E+12, 1.1E+10, 5.8E+07, 8.4E+07, 2.2E+08, 7.5E+07, 1.7E+05, 1.9E+07, 2.0E+05, 3.7E+05, 1.7+07, 5.8E+05, 1.2E+06, 5.0E+07, 1.7E+08, 2.4E+05, 2.0E+05, 1.0E+06, 2.4E+07, 1.0E+05 };
  
  // The array of the powers beta:
  float beta[nr] = {-1.5, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.64, 0.0, 0.5, 0.21, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5 };

  // The array of activation temperature:
  double Ta[nr] = { 59386.0, 51987.0, 59340.0, 50830.0, 32100.0, 25000.0, 32712.0, 9059.0, 10116.0, 8455.0, 5586.0, 21137.0, 24100.0, 362960.0, 27840.0, 24232.0, 28686.0, 39815.0, 37940.0, 24500.0, 19200.0, 15500.0, 22800.0, 14500.0, 6000.0 };

  // Array of third body efficiencies
  double M[nr] = {1.0, 1.0 , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  
  //---------------------------------------------------------------------------------//
  
  //------------------------------ Arrhenuis Law ------------------------------------//

  clock_t start_time = clock();
  
  double R = 0.0;
  for (int k=0; k<nsp; ++k) R += Yk[k]/wk[k];
  R *= Ru;
  double rho = P/R/T;
  

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

  /*
  // Displaying Kc:
  printf("\n\nEquilibrium Constant Kc:\n");
  for (int r=0; r<nr; ++r)
    printf("Kc,%d is %g \n", r, KcEq[r]);
  */
  
  // Calculating kb,r:
  double kback[nr];
  printf("\n Backward coefficients: \n");
  for (int r=0; r<nr; ++r) {
    kback[r] = kcoeff[r]/KcEq[r];
    //printf("kb,%d: %g \n", r, kback[r]); 
  }

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
   clock_t end_time = clock();

   // Calculate elapsed time in seconds
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

  // Printing results:
  for (int k=0; k<nsp; ++k) printf("Production rate of %3s is %e \n", names[k], omega[k]);
  printf("Execution time: %f seconds\n", elapsed_time);
  

  
  return 0;


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
