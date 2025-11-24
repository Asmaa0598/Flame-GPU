// Here is the main code:

/*

The test case is the Zeldovich mechanism provided in the paper:
"Computational analysis of the Extended Zeldovich Mechanism" by Lucky Anetor, Christopher Odetunde, Edward E. Osakue

Species: N NO N2 O O2 H OH

Global reaction: N2+O2->2NO

Version 1:
Reactions:
 N+NO -> N2+O
 N2+O -> N+NO
 N+O2 -> NO+O
 NO+O -> N+O2
 N+OH -> NO+H
 NO+H -> N+OH

k1 = 1.8e14 exp(-38370/T) m^3/(mol.s)
k2 = 3.8e13 exp(-425/T) m^3/(mol.s)
k3 = 1.8e10 exp(-4680/T) m^3/(mol.s)
k4 = 3.81e09 exp(-20820/T) m^3/(mol.s)
k5 = 7.1e13 exp(-450/T) m^3/(mol.s)
k6 = 1.7e14 exp(-24560/T) m^3/(mol.s)


Version 2:
Reactions:
 N+NO = N2+O
 N+O2 = NO+O
 N+OH = NO+H

k1 = 1.8e14 exp(-38370/T) m^3/(mol.s)
k2 = 1.8e10 exp(-4680/T) m^3/(mol.s)
k3 = 7.1e13 exp(-450/T) m^3/(mol.s)

H                 L 6/94H   1    0    0    0G   200.000  6000.000 1000.        1  
 0.25000000E+01 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
 0.25473660E+05-0.44668285E+00 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.25473660E+05-0.44668285E+00 0.26219035E+05    4

OH-               g 4/02O  1.H  1.E  1.   0.G   298.150  6000.000 1000.        1
 2.80023747E+00 1.13380509E-03-2.99666184E-07 4.01911483E-11-1.78988913E-15    2
-1.82535298E+04 4.69394620E+00 3.43126659E+00 6.31146866E-04-1.92914359E-06    3
 2.40618712E-09-8.66679361E-13-1.85085918E+04 1.07990541E+00-1.74702052E+04    4

N                 L 6/88N   1    0    0    0G   200.000  6000.000 1000.        1  
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226244E-10-0.20360983E-14    2
 0.56133775E+05 0.46496095E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104638E+05 0.41939088E+01 0.56850013E+05    4

NO                RUS 89N   1O   1    0    0G   200.000  6000.000 1000.        1  
 3.26071234E+00 1.19101135E-03-4.29122646E-07 6.94481463E-11-4.03295681E-15    2
 9.92143132E+03 6.36900518E+00 4.21859896E+00-4.63988124E-03 1.10443049E-05    3
-9.34055507E-09 2.80554874E-12 9.84509964E+03 2.28061001E+00 1.09770882E+04    4

N2                G 8/02N  2.   0.   0.   0.G   200.000  6000.000 1000.        1
 2.95257637E+00 1.39690040E-03-4.92631603E-07 7.86010195E-11-4.60755204E-15    2
-9.23948688E+02 5.87188762E+00 3.53100528E+00-1.23660988E-04-5.02999433E-07    3
 2.43530612E-09-1.40881235E-12-1.04697628E+03 2.96747038E+00 0.00000000E+00    4

O                 L 1/90O   1    0    0    0G   200.000  6000.000 1000.        1  
 2.54363697E+00-2.73162486E-05-4.19029520E-09 4.95481845E-12-4.79553694E-16    2
 2.92260120E+04 4.92229457E+00 3.16826710E+00-3.27931884E-03 6.64306396E-06    3
-6.12806624E-09 2.11265971E-12 2.91222592E+04 2.05193346E+00 2.99687009E+04    4

O2                RUS 89O   2    0    0    0G   200.000  6000.000 1000.        1  
 3.66096065E+00 6.56365811E-04-1.41149627E-07 2.05797935E-11-1.29913436E-15    2
-1.21597718E+03 3.41536279E+00 3.78245636E+00-2.99673416E-03 9.84730201E-06    3
-9.68129509E-09 3.24372837E-12-1.06394356E+03 3.65767573E+00 0.00000000E+00    4

HO2               T 1/09H  1.O  2.   0.   0.G   200.000  5000.000 1000.        1
 4.17228741E+00 1.88117627E-03-3.46277286E-07 1.94657549E-11 1.76256905E-16    2
 3.10206839E+01 2.95767672E+00 4.30179807E+00-4.74912097E-03 2.11582905E-05    3
-2.42763914E-08 9.29225225E-12 2.64018485E+02 3.71666220E+00 1.47886045E+03    4

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

int main() {

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

  // Array of the NASA polynomials upper range coefficients:
  /*
  const double nasaUpper[nsp][nasaC] = { 
 { 88765.013800, -107.123150, 2.362188, 0.000292, -0.000000, 0.000000, -0.000000, 56973.513300, 4.865236 }, 
 { 223901.871600, -1289.651623, 5.433936, -0.000366, 0.000000, -0.000000, 0.000000, 17503.176560, -8.501667 }, 
 { 587712.406000, -2239.249073, 6.066949, -0.000614, 0.000000, -0.000000, 0.000000, 12832.104150, -15.866396 }, 
 { 261902.026200, -729.872203, 3.317177, -0.000428, 0.000000, -0.000000, 0.000000, 33924.280600, -0.667959 }, 
 { -1037939.022000, 2344.830282, 1.819732, 0.001268, -0.000000, 0.000000, -0.000000, -16890.109290, 17.387165 }, 
 { 60.787743, -0.181935, 2.500212, -0.000000, 0.000000, -0.000000, 0.000000, 25474.863980, -0.448192 }, 
 { 1017393.379000, -2509.957276, 5.116548, 0.000131, -0.000000, 0.000000, -0.000000, 20444.871300, -11.012823 } 
 };
  */
  const double nasaUpper[nsp][nasaC] = { 
 { 88765.0138, -107.12315, 2.362188287, 0.0002916720081, -1.7295151e-07, 4.01265788e-11, -2.677227571e-15, 56973.5133, 4.86523579 }, 
 { 223901.8716, -1289.651623, 5.43393603, -0.00036560349, 9.88096645e-08, -1.416076856e-11, 9.38018462e-16, 17503.17656, -8.50166709 }, 
 { 587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679e-07, -1.923105485e-11, 1.061954386e-15, 12832.10415, -15.86639599 }, 
 { 261902.0262, -729.872203, 3.31717727, -0.000428133436, 1.036104594e-07, -9.43830433e-12, 2.725038297e-16, 33924.2806, -0.667958535 }, 
 { -1037939.022, 2344.830282, 1.819732036, 0.001267847582, -2.188067988e-07, 2.053719572e-11, -8.19346705e-16, -16890.10929, 17.38716506 }, 
 { 60.7877425, -0.1819354417, 2.500211817, -1.226512864e-07, 3.73287633e-11, -5.68774456e-15, 3.410210197e-19, 25474.86398, -0.448191777 }, 
 { 1017393.379, -2509.957276, 5.11654786, 0.000130529993, -8.28432226e-08, 2.006475941e-11, -1.556993656e-15, 20444.8713, 0 } 
 }; 

  // Vector of enthalpy of formation:
  const double Hf[nsp] = {0.56850013E+05, 1.09770882E+04, 0.00000000E+00, 2.99687009E+04, 0.00000000E+00, 0.26219035E+05, -1.74702052E+04};

  
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
  
  //---------------------------------------------------------------------------------//

  //---------------------------- Initial conditions ---------------------------------//
  double T = 1900.0;
  double P = 1.1065E+05;

  //double Yk[nsp] = { 2.31141, 2.31141, 0.0, 2.38145, 0.0, 0.0, 0.0 };
  double Yk[nsp] = { 0.231527, 0.495993, 0.0, 0.27248, 0.0, 0.0, 0.0 };
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
  printf("\n\nEquilibrium Constant Kc:\n");
  for (int r=0; r<nr; ++r)
    printf("Kc,%d is %g \n", r, KcEq[r]);
  
  // Calculating kb,r:
  double kback[nr];
  printf("\n Backward coefficients: \n");
  for (int r=0; r<nr; ++r) {
    kback[r] = kcoeff[r]/KcEq[r];
    printf("kb,%d: %g \n", r, kback[r]); 
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
  for (int k=0; k<nsp; ++k) printf("Production rate of %3s is %g \n", names[k], omega[k]);
  

  
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
