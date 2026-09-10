#ifndef MOTIF_EVIDENCE_H
#define MOTIF_EVIDENCE_H

#include <cstddef>
#include <string>

struct MotifEvidence {
  bool present = false;
  std::string matrixAlias;
  std::string matrixId;
  std::string matrixName;
  std::string matrixSource;
  double rawScore = 0.0;
  double scorePercent = 0.0;

  bool hasSourceRegion = false;
  std::string sourceRegionName;
  std::string sourceRegionType;
  size_t sourceRegionStart = 0;
  size_t sourceRegionEnd = 0;
  size_t relativeStart = 0;
};

#endif
