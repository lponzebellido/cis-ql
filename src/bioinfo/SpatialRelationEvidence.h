#ifndef SPATIAL_RELATION_EVIDENCE_H
#define SPATIAL_RELATION_EVIDENCE_H

#include <cstddef>
#include <string>

struct SpatialRelationEvidence {
  bool present = false;
  std::string relation;
  std::string referenceSet;
  std::string referenceChr;
  size_t referenceStart = 0;
  size_t referenceEnd = 0;
  std::string referenceStrand;
  std::string referenceType;
  std::string referenceName;
  size_t distance = 0;
  size_t maximumDistance = 0;
  bool overlaps = false;
};

#endif
