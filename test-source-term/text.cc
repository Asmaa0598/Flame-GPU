// --------------------------------------------------------------------------------
/* ac3804 Notes:
- Milestone 1: Completed
--> Writing a code that takes the string representing the elements and the species, counts the number of species, count how many elements in each species, outputs that in terminal.
- Milestone 2: Completed
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
#include <cmath>

using namespace std;

void getnasa9(std::array<double,9> &myarr, string coeff);
void copyTOoutput(string& theLine, vector<string> token, vector<string> toReplace, int index);

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
		// are the enthalpy of formation and the molecular weight:
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
		// Move to the next line, the first entry is the lower temperature
		// the second entry is the higher temperature
		getline(thermofile, thermline);
		++linenum;// !!! //
		moleW.str(thermline);
		moleW >> entry;
		double TL = stod(entry);
		moleW >> entry;
		double TMID = stod(entry);
		// Reverse the line to get the enthalpy of formation difference with
		// respect to 0K
		reverse(thermline.begin(), thermline.end());
		moleW.str(thermline);
		moleW >> entry;
		reverse(entry.begin(), entry.end());
		double hF0 = stod(entry); 
		// Move to the next line, each coefficient is composed of 15 characters
		// before the the exponent part there is a D (2 lines in total)
		getline(thermofile, thermline);
		++linenum;
		string coeffs = thermline;
		getline(thermofile, thermline);
		++linenum;
		coeffs += thermline;
		std::array<double, 9> nasaLOW;
		// Function to scan and store the 9 coefficients
		getnasa9(nasaLOW, coeffs);
		// Move to the next line to scan the upper range of temperature
		getline(thermofile, thermline);
		++linenum;
		moleW.str(thermline);
		moleW >> entry;
		moleW >> entry;
		double TU = stod(entry);
		// Move to the next line and scan the upper range of the NASA coefficient
		getline(thermofile, thermline);
		++linenum;
		coeffs = thermline;
		getline(thermofile, thermline);
		++linenum;
		coeffs += thermline;
		std::array<double, 9> nasaUP;
		
		mySpecies[k].setNASAup(nasaUP);
		mySpecies[k].setNASAlow(nasaLOW);
		mySpecies[k].setTup(TU);
		mySpecies[k].setTmid(TMID);
		mySpecies[k].setTlow(TL);
		mySpecies[k].setHf0(hF0);
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
    
    string names = {" {"};
    string wks = { " { " };
    string Tnasa = { "{ \n" };
    string nasa9low = { "{ \n" };
    string nasa9up = { "{ \n" };
    for (int k=0; k<numSP-1; ++k) {
      names = names + "\"" + mySpecies[k].getname() + "\",";
      wks = wks + to_string(mySpecies[k].getwk()) + ", ";
      Tnasa += "  { " + to_string(mySpecies[k].getTlow()) + ", " +
	to_string(mySpecies[k].getTmid()) + ", " +
	to_string(mySpecies[k].getTup()) +" }, \n";
      nasa9low += " { " ;
      for (int i=0; i<9; ++i)
	nasa9low += to_string(mySpecies[k].getNASAlowID(i)) + ", ";
      nasa9low += " }, \n";
    }
    names = names + "\"" + mySpecies[numSP-1].getname()+ "\"}; ";
    wks = wks + to_string(mySpecies[numSP-1].getwk()) + " }; ";
    Tnasa += "  { " + to_string(mySpecies[numSP-1].getTlow()) + ", " +
	to_string(mySpecies[numSP-1].getTmid()) + ", " +
	to_string(mySpecies[numSP-1].getTup()) +" } \n }; ";
    nasa9low += " { " ;
    for (int i=0; i<9; ++i)
	nasa9low += to_string(mySpecies[numSP-1].getNASAlowID(i)) + ", ";
    nasa9low += " } \n }; ";
    

    variables.push_back(to_string(numSP));
    variables.push_back(to_string(numR));
    variables.push_back(names);
    variables.push_back(wks);
    variables.push_back(Tnasa);
    variables.push_back(nasa9low);
    
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

// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// -------- Functions
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------

void getnasa9(std::array<double,9> &myarr, string coeff) {
  for(int i=0; i<10*16; i += 16) {
       string base{};
       string exponent{};
       bool divider = true;
	  if (i== 7*16)
	     continue;
		  
       for(int  j=i; j<i+16; ++j) {
	     if (coeff[j] == 'D') {
		 divider = false;
		 continue;
	     }
	     else if (coeff[j] == ' ')
		 continue;

	     if (divider)
		base += coeff[j];
	     else
		exponent += coeff[j];
	}
        myarr[i/16] = stod(base) * pow(10,stod(exponent));
   }

}

// Function to modify the line containg the token:
void copyTOoutput(string& theLine, vector<string> token, vector<string> toReplace, int index) {

  if (theLine.find(token[index]) != string::npos) {
    while (theLine.find(token[index]) != string::npos) 
	theLine.replace(theLine.find(token[index]), token[index].length(), toReplace[index]);
  }      
}

/*
  // Displaying stored elements
	for (const std::string& s : myElements) {
            std::cout << s << std::endl;
        }
 cout << "Number of species is " << numSP << endl;

 	for(int i=0; i<10*16; i += 16) {
		  string base{};
		  string exponent{};
		  bool divider = true;
		  if (i== 7*16)
		    continue;
		  
		  for(int  j=i; j<i+16; ++j) {
		    if (coeffs[j] == 'D') {
		      divider = false;
		      continue;
		    }
		    else if (coeffs[j] == ' ')
		      continue;

		    if (divider)
		      base += coeffs[j];
		    else
		      exponent += coeffs[j];
		  }
		  nasaLOW[i/16] = stod(base) * pow(10,stod(exponent));
		  cout << nasaLOW[i/16] << "   ";
		  //cout << base << exponent << "   ";
		}
*/
