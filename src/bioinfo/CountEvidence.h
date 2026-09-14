#ifndef COUNT_EVIDENCE_H
#define COUNT_EVIDENCE_H

#include <cstddef>
#include <string>

struct CountEvidence {
  bool present = false;
  std::string relation;
  std::string countedSet;
  std::string containerSet;
  size_t count = 0;
};

#endif
