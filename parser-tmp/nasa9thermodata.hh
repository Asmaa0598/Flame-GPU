#ifndef NASA9_THERMO_DATA_HH
#define NASA9_THERMO_DATA_HH

#include <cstdio>

struct NASA9ThermoData {
  struct ChemicalFormula {
    int numElements;
    char element[5][3];
    double valency[5];
  };

  char name[16-0];
  char description[80-18];
  int numIntervals;
  char referenceDateCode[9-3];
  ChemicalFormula chemicalFormula;
  int isGas;
  double molecularWeight;
  double heatOfFormation;
  double temperatureRange[5][2];
  double H298;
  double coefficients[5][9];
};

int readNASA9ThermoDataFromFile(
  FILE * fp, NASA9ThermoData ** alldata, int * speciesCount,
  char ** message
);

void print(FILE * fp, NASA9ThermoData *alldata, int nSpecies);

#endif
