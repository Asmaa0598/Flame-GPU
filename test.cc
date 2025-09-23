// --------------------------------------------------------------------------------
/* ac3804 Notes:
- Milestone 1:
--> Writing a code that takes the string representing the elements and the species, counts the number of species, count how many elements in each species, outputs that in terminal.
- Milestone 2:
--> Writing a code that will take the stored species and look into the thermo data to retrieve the molecular weight and the NASA polynomials.  
 */
// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------
/* Things to keep in mind:
  - The delimeter in each line labeled "1" is the space character.
  - Before "1" there is a number representing the name of the species.
  - The species maximum number of characters is 16.
  - The lines labeled at the end "2" or "3" or "4" store the NASA coefficients. The last coefficient in line "4" is enthalpy of formation.
  - There is no delimiter for the lines of NASA coefficients, so partitioning is done through number of characters, each 15 characters represents one coefficient 
 */
// --------------------------------------------------------------------------------

// --------------------------------------------------------------------------------
/* ac3804:
   There are 4 sections in the reaction mechanism file:
   UNITS
   ELEMENTS - SPECIES  - THERMO  - REACTIONS
 */
// --------------------------------------------------------------------------------

#include "myObjects.h"

using namespace std;


int main() {

  // The vector of variables to be filled:
  vector<string> variables;

  // My test input file:
  const string sample = {"sample.txt"};
  const string thermdata = {"nasa9.dat"};
  ifstream mechfile;
  ifstream thermofile;
  mechfile.open(sample);

  if (mechfile.is_open()) {
    cout << "File opened successfully." << endl;

    // Things to store:
    vector<string> myElements;
    vector<Species> mySpecies;
    int numEL = 0;
    int numSP = 0;
    int numR= 0; // PLACEHOLDER
    vector<int> speAvail; 

    // Reading the lines in the file:
    string line;
    while ( getline(mechfile, line) ) {

      // Extracting the chemical elements
      if (line == "ELEMENTS") {
	while (getline(mechfile, line) && line != "END" ){

	  string entry;
	  stringstream ss(line);
	  
	  while ( ss >> entry ){
	    myElements.push_back(entry);
	    ++numEL;
	  }
	}
	
      }
      
      else if (line == "SPECIES") {

	// Read the involved species
	while (getline(mechfile, line) && line != "END" ){

	  string entry;
	  stringstream ss(line);
	  
	  while ( ss >> entry ){
	    //entry = entry + ' ';
	    Species sp(entry, myElements, numEL);
	    mySpecies.push_back(sp);
	    speAvail.push_back(0.0);
	    ++numSP;
	  }  
	  }

	// Fetch the molecular weight and NASA coefficients and the enthalpy of 
	// formation of the species at 298.15 K from the thermal data
	thermofile.open(thermdata);
	string thermline;
	// counting lines:
	int linenum = 0;
	while (getline(thermofile, thermline)) {
	  stringstream ss(thermline);
	  string entry;
	  ++linenum;
	  ss >> entry;
	  for(int k=0; k<numSP ; ++k) {
	      if (entry == mySpecies[k].getname()) {
		cout << "Found " << entry << " at " << linenum << endl;
		// Move to the next line, the 1st and 2nd entry from the right
		// are the enthalpy of dormation and the molecular weight:
		getline(thermofile, thermline);
		++linenum;
		reverse(thermline.begin(), thermline.end());  
		stringstream moleW(thermline);
		moleW >> entry;
		reverse(entry.begin(), entry.end());
		double hFtest = stod(entry); 
		moleW >> entry;
		reverse(entry.begin(), entry.end());
		double wTest = stod(entry);
		
		mySpecies[k].setHf(hFtest);
		mySpecies[k].setwk(wTest);
		speAvail[k] = 1.0;
		break;
	      }
	    }
	   
	}
	thermofile.close(); 
      }
    }

    // Create the variable that are going to replace the token:
    // The number of species:
    variables.push_back(to_string(numSP));
    variables.push_back(to_string(numR));

    string names = {" {"};
    for(int k=0; k<numSP-1 ; ++k)
      names = names + "\"" + mySpecies[k].getname() + "\",";
    names = names + "\"" + mySpecies[numSP-1].getname()+ "\"}; ";

    variables.push_back(names);

    cout << "strings to replace the tokens: " << endl;
    for (const std::string& v : variables) {
            std::cout << v << std::endl;
    }

    cout << "\n\nAvailable Species in thermo data: " << endl;
    for (int k=0; k<numSP; ++k) {
      if (speAvail[k] == 1.0 ) {
	cout << " Species " << mySpecies[k].getname() << " was found in " << thermdata << endl;
	cout << "  The molecular weight is " << mySpecies[k].getwk() << " kg/mol"
	     << " and its enthalpy of formation is " << mySpecies[k].getHf() << " J/mol" << endl;
      }
      else
	cout << " Species " << mySpecies[k].getname() << " was not found in " << thermdata << endl;
    }
    

  mechfile.close();

  }
  
  // There is no file
  else {
    cout << "No file named " << sample << " is found!" << endl;
  }

  return 0;
}

/*
  // Displaying stored elements
	for (const std::string& s : myElements) {
            std::cout << s << std::endl;
        }
 cout << "Number of species is " << numSP << endl;
*/
