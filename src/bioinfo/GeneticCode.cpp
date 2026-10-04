#include "GeneticCode.h"
#include <cctype>
#include <map>

namespace {

const std::map<int, std::string> &tables() {
  static const std::map<int, std::string> values = {
      {1, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {2, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"},
      {3, "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {4, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {5, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"},
      {6, "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {9, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
      {10, "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {11, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {12, "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {13, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"},
      {14, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
      {15, "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {16, "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {21, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"},
      {22, "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {23, "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {24, "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"},
      {25, "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {26, "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {27, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {28, "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {29, "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {30, "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {31, "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {32, "FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"},
      {33, "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"}};
  return values;
}

int baseIndex(char base) {
  switch (std::toupper(static_cast<unsigned char>(base))) {
  case 'T':
  case 'U':
    return 0;
  case 'C':
    return 1;
  case 'A':
    return 2;
  case 'G':
    return 3;
  default:
    return -1;
  }
}

}

bool GeneticCode::supports(int table) { return tables().count(table) > 0; }

std::vector<int> GeneticCode::supportedTables() {
  std::vector<int> result;
  for (const auto &entry : tables())
    result.push_back(entry.first);
  return result;
}

bool GeneticCode::translate(const std::string &sequence, int table,
                            size_t frame, std::string &protein,
                            std::string &error) {
  protein.clear();
  error.clear();
  const auto found = tables().find(table);
  if (found == tables().end()) {
    error = "Unsupported NCBI genetic code table " + std::to_string(table) +
            ".";
    return false;
  }
  if (frame > 2) {
    error = "Translation frame must be 0, 1, or 2.";
    return false;
  }
  if (found->second.size() != 64) {
    error = "Invalid internal genetic code table.";
    return false;
  }
  for (size_t offset = frame; offset + 2 < sequence.size(); offset += 3) {
    const int first = baseIndex(sequence[offset]);
    const int second = baseIndex(sequence[offset + 1]);
    const int third = baseIndex(sequence[offset + 2]);
    if (first < 0 || second < 0 || third < 0) {
      protein += 'X';
      continue;
    }
    protein += found->second[static_cast<size_t>(first * 16 + second * 4 +
                                                 third)];
  }
  return true;
}
