#include "FastaReader.h"
#include <fstream>
#include <iostream>
#include <sstream>

std::vector<FastaRecord> FastaReader::read(const std::string& filename) {
  std::vector<FastaRecord> records;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Runtime Error: Cannot open FASTA file '" << filename << "'." << std::endl;
    return records;
  }

  std::string line;
  FastaRecord current;
  bool inRecord = false;

  while (std::getline(file, line)) {
    if (line.empty()) continue;
    if (line[0] == '>') {
      if (inRecord) {
        records.push_back(current);
      }
      current = FastaRecord();
      current.header = line.substr(1);
      std::istringstream iss(current.header);
      iss >> current.sequenceId;
      inRecord = true;
    } else if (inRecord) {
      for (char c : line) {
        if (c != '\r' && c != '\n' && c != ' ') {
          current.sequence += std::toupper(c);
        }
      }
    }
  }
  if (inRecord) {
    records.push_back(current);
  }

  return records;
}
