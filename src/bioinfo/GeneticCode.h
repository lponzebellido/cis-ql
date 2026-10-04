#ifndef GENETIC_CODE_H
#define GENETIC_CODE_H

#include <cstddef>
#include <string>
#include <vector>

class GeneticCode {
public:
  static bool supports(int table);
  static std::vector<int> supportedTables();
  static bool translate(const std::string &sequence, int table, size_t frame,
                        std::string &protein, std::string &error);
};

#endif
