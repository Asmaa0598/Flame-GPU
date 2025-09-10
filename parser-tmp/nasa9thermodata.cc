#include <nasa9thermodata.hh>

#include <string>
#include <cstdlib>
#include <cstring>
#include <cmath>

int trimcopy(char const * begin, char const * end, char * buffer) {
  if(begin < end) {
    while(std::isspace(*begin)) { ++begin; }
    do { --end; } while(std::isspace(*end));
  }
  int cnt = 0;
  while(begin <= end) {
    buffer[cnt++] = *begin++;
  }
  buffer[cnt] = '\0';
  return cnt;
}

int toDouble(char const * begin, char const * end, double * value) {
  char * p;
  *value = std::strtod(begin, &p);
  return errno == ERANGE;
}

int toInt(char const * begin, char const * end, int * value) {
  char * p;
  *value = std::strtol(begin, &p, 10);
  return errno == ERANGE;
}

int readNASA9ThermoDataLine1(
  char const * line, NASA9ThermoData * data, char ** message
) {
  trimcopy(&line[0], &line[16], data->name);
  trimcopy(&line[18], &line[80], data->description);
  return 0;
}

int readNASA9ThermoDataLine2(
  char const * line, NASA9ThermoData * data, char ** message
) {
  if(toInt(&line[0], &line[2], &data->numIntervals)) {
    char const * tmp = "Could not read number of temperature intervals";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
    std::strcpy(*message, tmp);
    return 1;
  }

  trimcopy(&line[3], &line[9], data->referenceDateCode);

  data->chemicalFormula.numElements = 0;
  for(int i = 0; i < 5; ++i) {
    int elementNameSize = trimcopy(
      &line[10+i*8], &line[12+i*8], data->chemicalFormula.element[i]
    );
    if(toDouble(&line[12+i*8], &line[18+i*8], &data->chemicalFormula.valency[i])) {
      char const * tmp = "Could not read element valency";
      *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
      std::strcpy(*message, tmp);
      return 2;
    }
    if(elementNameSize > 0) {
      data->chemicalFormula.numElements++;
    }
  }

  if(toInt(&line[50], &line[52], &data->isGas)) {
    char const * tmp = "Could not read phase";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
    std::strcpy(*message, tmp);
    return 3;
  }

  if(toDouble(&line[52], &line[65], &data->molecularWeight)) {
    char const * tmp = "Could not read molecular weight";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
    std::strcpy(*message, tmp);
    return 4;
  }

  if(toDouble(&line[65], &line[80], &data->heatOfFormation)) {
    char const * tmp = "Could not read heat of formation";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
    std::strcpy(*message, tmp);
    return 5;
  }

  return 0;
}

int readNASA9ThermoDataLine3(
  char const * line, int const i, NASA9ThermoData * data, char ** message
) {
  if(toDouble(&line[0], &line[11], &data->temperatureRange[i][0])) {
    char const * tmp = "Could not read start of temperature range";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 1;
  }

  if(toDouble(&line[11], &line[22], &data->temperatureRange[i][1])) {
    char const * tmp = "Could not read end of temperature range";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 2;
  }

  if(toDouble(&line[65], &line[80], &data->H298)) {
    char const * tmp = "Could not read H298";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 3;
  }

  return 0;
}

int readNASA9ThermoDataLine4(
  char const * line, int const i, NASA9ThermoData * data, char ** message
) {
  char tmpline[80];
  std::strncpy(tmpline, line, 80);
  for(int k = 0; k < 80; ++k) {
    if(tmpline[k] == 'D' || tmpline[k] == 'd') {
      tmpline[k] = 'E';
    }
  }

  if(toDouble(&tmpline[0], &tmpline[16], &data->coefficients[i][0])) {
    char const * tmp = "Could not read coefficient 1";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 1;
  }

  if(toDouble(&tmpline[16], &tmpline[32], &data->coefficients[i][1])) {
    char const * tmp = "Could not read coefficient 2";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 2;
  }

  if(toDouble(&tmpline[32], &tmpline[48], &data->coefficients[i][2])) {
    char const * tmp = "Could not read coefficient 3";
    *message - (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 3;
  }

  if(toDouble(&tmpline[48], &tmpline[64], &data->coefficients[i][3])) {
    char const * tmp = "Could not read coefficient 4";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 4;
  }

  if(toDouble(&tmpline[64], &tmpline[80], &data->coefficients[i][4])) {
    char const * tmp = "Could not read coefficient 5";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 5;
  }

  return 0;
}

int readNASA9ThermoDataLine5(
  char const * line, int const i, NASA9ThermoData * data, char ** message
) {
  char tmpline[80];
  std::strncpy(tmpline, line, 80);
  for(int k = 0; k < 80; ++k) {
    if(tmpline[k] == 'D' || tmpline[k] == 'd') {
      tmpline[k] = 'E';
    }
  }

  if(toDouble(&tmpline[0], &tmpline[16], &data->coefficients[i][5])) {
    char const * tmp = "Could not read coefficient 6";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 1;
  }

  if(toDouble(&tmpline[16], &tmpline[32], &data->coefficients[i][6])) {
    char const * tmp = "Could not read coefficient 7";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 2;
  }

  if(toDouble(&tmpline[48], &tmpline[64], &data->coefficients[i][7])) {
    char const * tmp = "Could not read coefficient 8";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 3;
  }

  if(toDouble(&tmpline[64], &tmpline[18], &data->coefficients[i][8])) {
    char const * tmp = "Could not read coefficient 9";
    *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+20));
    sprintf(*message, "%s: %d", tmp, i);
    return 4;
  }

  return 0;
}

int readNASA9ThermoDataFromFile(
  FILE * fp, NASA9ThermoData ** alldata, int * speciesCount,
  char ** message
) {
  int capacity = 100;
  *alldata = (NASA9ThermoData*)malloc(sizeof(NASA9ThermoData)*capacity);
  
  char line[82];
  size_t lineSize = 82;
  NASA9ThermoData data;

  *speciesCount = 0;

  while(std::fgets(line, lineSize, fp) != nullptr) {
    if(line[0] == '!') continue;

    if(std::strncmp(line, "END PRODUCTS", std::strlen("END PRODUCTS")) == 0) {
      continue;
    }

    if(std::strncmp(line, "END REACTANTS", std::strlen("END REACTANTS")) == 0) {
      continue;
    }

    std::memset(&data, sizeof(NASA9ThermoData), 0);

    if(readNASA9ThermoDataLine1(line, &data, message)) {
      return 1;
    }

    if(std::strcmp(data.name, "thermo") == 0) {
      std::fgets(line, lineSize, fp);
      continue;
    }

    if(std::fgets(line, lineSize, fp) != nullptr) {
      if(readNASA9ThermoDataLine2(line, &data, message)) {
        return 2;
      }
    } else {
      char const * tmp = "Truncated record";
      *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
      std::strcpy(*message, tmp);
      return 2;
    }

    if(data.numIntervals == 0) {
      if(std::fgets(line, lineSize, fp) != nullptr) {
        if(readNASA9ThermoDataLine3(line, 0, &data, message)) {
          return 3;
        }
      } else {
        char const * tmp = "Truncated record";
        *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
        std::strcpy(*message, tmp);
        return 3;
      }
    } else {
      for(int i = 0; i < data.numIntervals; ++i) {
        if(std::fgets(line, lineSize, fp) != nullptr) {
          if(readNASA9ThermoDataLine3(line, i, &data, message)) {
            return 3+i*3;
          }
        } else {
          char const * tmp = "Truncated record";
          *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
          std::strcpy(*message, tmp);
          return 3+i*3;
        }
      
        if(std::fgets(line, lineSize, fp) != nullptr) {
          if(readNASA9ThermoDataLine4(line, i, &data, message)) {
            return 4+i*3;
          }
        } else {
          char const * tmp = "Truncated record";
          *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
          std::strcpy(*message, tmp);
          return 4+i*3;
        }
      
        if(std::fgets(line, lineSize, fp) != nullptr) {
          if(readNASA9ThermoDataLine5(line, i, &data, message)) {
            return 5+i*3;
          }
        } else {
          char const * tmp = "Truncated record";
          *message = (char*)malloc(sizeof(char)*(std::strlen(tmp)+1));
          std::strcpy(*message, tmp);
          return 5+i*3;
        }
      }
    }

    ++(*speciesCount);
    if(capacity < *speciesCount) {
      NASA9ThermoData * tmp = (NASA9ThermoData*)malloc(
        sizeof(NASA9ThermoData)*(capacity*2)
      );
      std::memcpy(tmp, *alldata, sizeof(NASA9ThermoData)*capacity);
      free(*alldata);
      *alldata = tmp;
      capacity *= 2;
    }
    std::memcpy(&(*alldata)[*speciesCount-1], &data, sizeof(NASA9ThermoData));
  }

  return 0;
}

void print(FILE * fp, NASA9ThermoData *data, int nSpecies) {
  for(int s = 0; s < nSpecies; ++s) {
    fprintf(fp, "Species Name: '%s'\n", data[s].name);
    fprintf(fp, "Description: '%s'\n", data[s].description);
    fprintf(fp, "Number of Temperature Intervals: %d\n", data[s].numIntervals);
    fprintf(fp, "Reference Date Code: '%s'\n", data[s].referenceDateCode);
    for(int i = 0; i < data[s].chemicalFormula.numElements; ++i) {
      fprintf(fp, "Element Name: '%s'\n", data[s].chemicalFormula.element[i]);
      fprintf(fp, "Element Valency: %lf\n", data[s].chemicalFormula.valency[i]);
    }
    fprintf(fp, "Is Gas?: %d\n", data[s].isGas);
    fprintf(fp, "Molecular Weight: %lf\n", data[s].molecularWeight);
    fprintf(fp, "Heat of Formation: %lf\n", data[s].heatOfFormation);
    for(int i = 0; i < data[s].numIntervals; ++i) {
      fprintf(fp, "Temperature range: [%lf, %lf]\n",
              data[s].temperatureRange[i][0], data[s].temperatureRange[i][1]);
      fprintf(fp, "H(298.15)-H(0): %lf\n", data[s].H298);
      fprintf(fp, "a1: %16.9le\n", data[s].coefficients[i][0]);
      fprintf(fp, "a2: %16.9le\n", data[s].coefficients[i][1]);
      fprintf(fp, "a3: %16.9le\n", data[s].coefficients[i][2]);
      fprintf(fp, "a4: %16.9le\n", data[s].coefficients[i][3]);
      fprintf(fp, "a5: %16.9le\n", data[s].coefficients[i][4]);
      fprintf(fp, "a6: %16.9le\n", data[s].coefficients[i][5]);
      fprintf(fp, "a7: %16.9le\n", data[s].coefficients[i][6]);
      fprintf(fp, "b1: %16.9le\n", data[s].coefficients[i][7]);
      fprintf(fp, "b2: %16.9le\n", data[s].coefficients[i][8]);
    }
  }
}
