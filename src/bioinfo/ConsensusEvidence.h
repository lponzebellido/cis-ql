#ifndef CONSENSUS_EVIDENCE_H
#define CONSENSUS_EVIDENCE_H

#include <cstddef>
#include <string>
#include <vector>

struct ConsensusEvidence {
  bool present = false;
  std::string anchorSet;
  size_t minimumSupport = 0;
  size_t observedSupport = 0;
  bool hasMinimumReciprocalOverlap = false;
  double minimumReciprocalOverlapPercent = 0.0;
  bool hasMaximumSummitDistance = false;
  size_t maximumSummitDistanceBp = 0;
  std::vector<std::string> inputSets;
};

#endif
