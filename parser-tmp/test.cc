#include <nasa9thermodata.hh>

#include <cstdlib>

int main(int argc, char * argv[]) {
  char const * filename = argv[1];
  FILE * fp = fopen(filename, "r");
  if(fp == nullptr) {
    fprintf(stderr, "Unable to open file: %s", filename);
    return 1;
  }

  NASA9ThermoData * data;
  int speciesCount;
  char * message;
  if(readNASA9ThermoDataFromFile(fp, &data, &speciesCount, &message)) {
    fprintf(stderr, "Error: %s\n", message);
    free(message);
    fclose(fp);
    return 2;
  }
  fprintf(stderr, "#species: %d\n", speciesCount);

  fclose(fp);

  return 0;
}
