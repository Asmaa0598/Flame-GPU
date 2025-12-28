// --------------------------------------------------------------------------------
// Here testing a code that reads a line that contains a reaction
//  1- test reading: H2+O2 = 2OH 0.170E+14 0.00 47780
//     Task to be performed:
//          - Identify whether it's a reaction or not
//          - Recognize reactant species and product species
//          - Recognize whether there is more than 3 reactants or products
//          - Recognize reversibility
//          - Recognize reaction order
//          - Recognize if there are theree Arrhenius law parameters
//          - Recognize if there is a foreign species
//          - Check elements balance
//          - Check charge balance
// --------------------------------------------------------------------------------


#include "reactionMech.h"
#include <iomanip>
#include <cmath>

using namespace std;

int main() {

  const int numEL = 3;
  const int numSP = 11;

  vector<string> elements;
  elements.push_back("H");
  elements.push_back("O");
  elements.push_back("N");

  string names[numSP] = { "H2", "H", "O2", "O", "OH", "HO2", "H2O2", "H2O", "N", "N2", "NO"};

  vector<Species> mySpecies;
  mySpecies.push_back(Species(names[0], elements, numEL));
  mySpecies.push_back(Species(names[1], elements, numEL));
  mySpecies.push_back(Species(names[2], elements, numEL));
  mySpecies.push_back(Species(names[3], elements, numEL));
  mySpecies.push_back(Species(names[4], elements, numEL));
  mySpecies.push_back(Species(names[5], elements, numEL));
  mySpecies.push_back(Species(names[6], elements, numEL));
  mySpecies.push_back(Species(names[7], elements, numEL));
  mySpecies.push_back(Species(names[8], elements, numEL));
  mySpecies.push_back(Species(names[9], elements, numEL));
  mySpecies.push_back(Species(names[10], elements, numEL));

  string dividers[5] = {"!", "<=>", "=>", "=", "/" };

  string line = "H2 +O2 = 2OH 0.170E+14 0.00 47780";

  string reactantClause;
  string productClause;
  string mClause; 
  bool reversibility = false;

  // Ignore the commented region
  size_t position = line.find(dividers[0]);
  line = line.substr(0, position);
  
  size_t position1 = line.find(dividers[1]);
  size_t position2 = line.find(dividers[2]);
  size_t position3 = line.find(dividers[3]);
  size_t position4 = line.find(dividers[4]);

 

  if (position1 != std::string::npos) {
    cout << "It's a reversible reaction" << endl;
    reversibility = true; 
    reactantClause = line.substr(0, position3-1);
    productClause =  line.substr(position3+dividers[4].length()+1);

    cout << "The reactant clause is: " << reactantClause << endl;
    cout << "The product clause is: " << productClause << endl;
  }
  else if (position2 != std::string::npos){
    cout << "It's a irreversible reaction" << endl;
    reactantClause = line.substr(0, position3);
    productClause =  line.substr(position3+dividers[4].length()+1);

    cout << "The reactant clause is: " << reactantClause << endl;
    cout << "The product clause is: " << productClause << endl;
  }
  else if (position3 != std::string::npos) {
    cout << "It's a reversible reaction" << endl;
    reversibility = true; 
    reactantClause = line.substr(0, position3);
    productClause =  line.substr(position3+dividers[4].length());

    cout << "The reactant clause is: " << reactantClause << endl;
    cout << "The product clause is: " << productClause << endl;
  }
  else if (position4 != std::string::npos) {
    cout << "It's the continuation of a reaction" << endl;
  }

  Reaction(1, numSP, mySpecies, reversibility, reactantClause, productClause); 
    
  return 0;
}
