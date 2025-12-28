#include "mainCopy.h"
#include <time.h>

int main() {
  srand((unsigned int)time(NULL));

    // Output files
    FILE* output;
    FILE* measurements;

    // Opening the files in write mode
    output = fopen("outputZeldovich.txt", "a");
    measurements = fopen("costZeldovich.txt", "a");

    // --------------------------------------------------------------------------------
    // number of iterations
    const int nIter = 1e+06;
    // ranges limits
    const double minT = 300.0;
    const double maxT = 3500.0;
    const double minP = 1.0e+04;
    const double maxP = 1.0e+08;
    // creating the temperature and pressure arrays:
    double T[nIter];
    double P[nIter];
    for (int i = 0; i < nIter; i++){
      T[i] = (double)rand() / (double)RAND_MAX; 
      T[i] = minT + T[i] * (maxT - minT); 
    }
    for (int i = 0; i < nIter; i++) {
      P[i] = (double)rand() / (double)RAND_MAX; 
      P[i] = minP + P[i] * (maxP - minP); 
    }
    // The mass fractions array:
    double Y[nsp] = { 0.231527, 0.495993, 0.0, 0.27248, 0.0, 0.0, 0.0 };
    // --------------------------------------------------------------------------------
    // Mechanism specific information 


    // --------------------------------------------------------------------------------
    double omegas[nsp] = {0.0};
    // Start time
    clock_t start_time = clock();
    // main loop
    for (int i=0; i<nIter; ++i) {
      zeldovich(omegas, T[i], P[i], &Y[0]);
    }
    clock_t end_time = clock();
    // Calculate elapsed time in seconds
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    // --------------------------------------------------------------------------------
    // writing in files:
    // writing cost:
    fprintf(measurements, "%d %g \n", nIter, elapsed_time);
    // writing omegas:
    for (int i = 0; i < nsp; i++) {
        fprintf(output, "%g ", omegas[i]);
    }
   fprintf(output, "\n");

    
  return 0; 
}
