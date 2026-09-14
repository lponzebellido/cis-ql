#ifndef MODULE_EVIDENCE_H
#define MODULE_EVIDENCE_H

#include "MotifEvidence.h"
#include <cstddef>
#include <string>

struct ModuleMemberEvidence {
  std::string sourceSet;
  std::string chr;
  size_t start = 0;
  size_t end = 0;
  std::string strand;
  std::string type;
  std::string name;
  MotifEvidence motifEvidence;
};

struct ModuleEvidence {
  bool present = false;
  size_t minimumSpacing = 0;
  size_t maximumSpacing = 0;
  size_t observedSpacing = 0;
  std::string orderPolicy;
  std::string observedOrder;
  std::string orientationPolicy;
  std::string observedOrientation;
  ModuleMemberEvidence first;
  ModuleMemberEvidence second;
};

#endif
