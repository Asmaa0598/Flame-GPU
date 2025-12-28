#include "myObjects.h"

// Here the reaction objects are stored.
// ********************************************************************************
// Things to be added:
//  - Unit conversion depending on reaction order
//  - set enums for reaction types: 0-Unimolecular, 1-bimolecular, 2-Trimolecular, 3-THird body, 4-Falloff
// ********************************************************************************

class Reaction {

  int index;

  enum ReactionType {
        Unimolecular,
	Bimolecular,
        Trimolecular,
	ThirdBody,
	Lindemann, // Falloff formulation
	TROE, // Falloff formulation
	SRI // Falloff formulation
    };
  ReactionType rtype;
  
  bool reversibility = true; // 0 for irreversible reactions;
  double A;
  double beta;
  double Ta; // Activation temperature, either obtained from file or from Activation
  // energy in K

  // Array storing reactants index and their reactant coeff, initially set to include number higher than 3.
  std::array<int, 2> reactCoeff[3] = {
        {1024, 0},
        {1024, 0},
        {1024, 0}
  };
  // Array storing products index and their product coeff
  std::array<std::array<int, 2>, 3> prodCoeff; 
  // Vector of M coefficients
  vector<array<int,2>> mefficiency;

 public:
  // Default constructor
  Reaction() { }
  

  // Main constructor
 Reaction(int id, int nsp, vector<Species> mySpecies, bool rever, string reactClause, string prodClause) : index(id), reversibility(rever) {
    cout << "Successfull construction!" << endl;
    cout << "reaction ID is " << index << endl; 
    // Parsing reactants:
    stringstream ss(reactClause);
    char delim = '+';
    string entry;
    int nbr =0;

    while (getline(ss, entry, delim) && !(nbr>2)) {
      cout << "I have: " << entry << endl;
      if (isdigit(entry[0])) {
	int reactcoeff = entry[0]-'0';
	entry.erase(0, 1);
	int k=0;
	while (k<nsp) {
	  if (entry.compare(mySpecies[k].getname())>=0) {
	    reactCoeff[nbr][0] = k;
	    reactCoeff[nbr][1] = reactcoeff;
	    nbr += reactcoeff;
	    break;
	  }
	  ++k;
	}
      } else {
	if (entry.compare("M")>=0) {
	  ++nbr;
	  rtype = ThirdBody;
	}
	else {
	  int k=0;
	  while (k<nsp) {
	    if (entry.compare(mySpecies[k].getname())>=0) {
	      reactCoeff[nbr][0] = k;
	      reactCoeff[nbr][1] = 1;
	      ++nbr;
	      break;
	    }
	    ++k;
	  }
	}
      }     
    }

    // Checking if reaction has more than 3 reactants
    if (nbr>3) {
      cout << "Reaction has more than 3 reactants!!! nbr=" << nbr << endl;
      abort();
    } else {
      cout << " We have " << nbr << " reactants" << endl;
    }
    
  }

  // Setters:
  // Setter for M efficiencies

  // Getter:
  
  
  
};
