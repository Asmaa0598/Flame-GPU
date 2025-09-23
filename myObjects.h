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
 
  double moleweight;
  double hf; // Enthalpy of formation at 298.75K in J/mol
  double nasaUp[9];
  double nasaLow[9]; 

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
  // -------------------- Getters --------------------//
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

  // -------------------- Setters --------------------
  // Set the species molecular weight:
  void setwk(double wK) {
    moleweight = 0.001*wK; // To covert to kg/mol
  }
  // Set the species enthalpy of formation:
  void setHf(double hF) { hf = hF; }
};

