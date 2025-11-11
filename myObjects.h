#include<iostream>
#include<fstream>
#include<sstream>
#include<filesystem>
#include<string>
#include<cctype>
#include<vector>
#include<tuple>
#include <bits/stdc++.h>

using namespace std;

class Species {
  string name;
  
  // Stores Element index and how many of it in the species
  vector<tuple<int,int>> numElements;
 
  double moleweight; // in kg/mol
  double hf; // Enthalpy of formation at 298.75K in J/mol
  double hf0; // Enthalpy of formation at 0K: H(298.15) - H(0) in J/mol
  std::array<double, 9> nasaUp;
  std::array<double, 9> nasaLow;
  double Tlow;
  double Tmid;
  double Tup;

public:
  // Default construtor, needs to be fixed, this is a memory disaster:
  Species() { }

  // Main constructor:
  /* This constructor takes for now only the name of the species, once the name
     and the number of elements has been established, it will also take the an
     ifstream object to parse the thermo file to load molecular weight, NASA
     coefficients, and enthalpy of formations.
  */
  Species(string& speMechanism, vector<string> elements, int numEL) {
    name = speMechanism;
    // Test:
    cout << name << endl;
    // Parsing name
      for (int e=0; e<numEL; ++e) {
	size_t found = name.find(elements[e]);
	if ( found != string::npos) {
	  if (isdigit(name[++found])){
	    char number = name[found];
	    numElements.push_back(make_tuple((int)e, atoi(&number)));
	  }
	  else
	    numElements.push_back(make_tuple((int)e, (int)1));
	}
    }
    
  }

  void printingElements() {
     // Verification:
    for (const auto& tuple : numElements) {
        // Access tuple elements using get<index>(tuple)
      cout << "It contains " << get<1>(tuple) << " of element " << get<0>(tuple) << endl;
    }
  }
  // -------------------- Getters -------------------- //
  //
  string getname() {  // Getter for specie's name
    return name;
  }
  double getwk() { // Getter for molecular weight
    return moleweight;
  }
  double getHf() { // Getter for Enthalpy of formation
    return hf;
  }
  double getH0f() { // Getter for Enthalpy of formation difference with respect to 0K
    return hf0;
  }
  double getTlow() { // Getter for NASA Polynomials lowest temperature
    return Tlow;
  }
  double getTmid() { // Getter for NASA Polynomials mid temperature
    return Tmid;
  }
  double getTup() { // Getter for NASA Polynomials highest temperature
    return Tup;
  }
  std::array<double, 9> getNASAlow() { // Getter to get low range NASA coefficients
    return nasaLow;
  }
  double getNASAlowID(int index) {
    // Getter to get low coefficients to index
    return nasaLow[index];
  }
  std::array<double, 9> getNASAup() { // Getter to get upper range NASA coefficients
    return nasaUp;
  }
  double getNASAupID(int index) {
    // Getter to get upper NASA coefficients to index
    return nasaUp[index];
  }

  // -------------------- Setters -------------------- //
  // Set the species molecular weight:
  void setwk(double wK) {
    moleweight = 0.001*wK; // To covert to kg/mol
  }
  // Set the species enthalpy of formation:
  void setHf(double hF) { hf = hF; }
  // Set the species enthalpy of formation difference with respect to 0K:
  void setHf0(double hF) { hf0 = hF; }
  // Set NASA Polynomial lowest temperature
  void setTlow(double Tl) { Tlow = Tl; }
  // Set NASA Polynomial mid temperature
  void setTmid(double Tm) { Tmid = Tm; }
  // Set NASA Polynomial highest temperature
  void setTup(double Th) { Tup = Th; }
  // Set the species Lower range of NASA coefficients
  void setNASAlow(std::array<double,9> low) {
    for(int i=0; i<9; ++i) nasaLow[i] = low[i];
  }
  // Set the species  range of NASA coefficients
  void setNASAup(std::array<double,9> up) {
    for(int i=0; i<9; ++i) nasaUp[i] = up[i];
  }
  
};

