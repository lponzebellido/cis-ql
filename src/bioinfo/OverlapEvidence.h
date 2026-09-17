#ifndef OVERLAP_EVIDENCE_H
#define OVERLAP_EVIDENCE_H

#include "TrackEvidence.h"
#include "ConsensusEvidence.h"
#include <cstddef>
#include <string>
#include <vector>

struct OverlapEvidence {
  std::string referenceSet;
  std::string referenceChr;
  size_t referenceStart = 0;
  size_t referenceEnd = 0;
  std::string referenceStrand;
  std::string referenceType;
  std::string referenceName;
  TrackEvidence trackEvidence;
  ConsensusEvidence consensusEvidence;
  std::vector<OverlapEvidence> supportingEvidence;
};

#endif
