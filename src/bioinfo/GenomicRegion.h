#ifndef GENOMIC_REGION_H
#define GENOMIC_REGION_H

#include "MotifEvidence.h"
#include "SpatialRelationEvidence.h"
#include "CountEvidence.h"
#include "ConsensusEvidence.h"
#include "ModuleEvidence.h"
#include "OverlapEvidence.h"
#include "TrackEvidence.h"
#include <string>
#include <algorithm>
#include <vector>

struct GenomicRegion {
  std::string chr;
  size_t start;
  size_t end;
  std::string strand;
  std::string type;
  std::string name;
  std::string sequence;
  MotifEvidence motifEvidence;
  SpatialRelationEvidence spatialRelation;
  CountEvidence countEvidence;
  ConsensusEvidence consensusEvidence;
  ModuleEvidence moduleEvidence;
  TrackEvidence trackEvidence;
  std::vector<OverlapEvidence> overlapEvidence;

  size_t length() const { return end - start; }

  bool overlaps(const GenomicRegion& other, bool matchStrand = false) const {
    if (chr != other.chr) return false;
    if (matchStrand && !strand.empty() && !other.strand.empty() && strand != other.strand) return false;
    return start < other.end && end > other.start;
  }

  bool operator<(const GenomicRegion& other) const {
    if (chr != other.chr) return chr < other.chr;
    if (start != other.start) return start < other.start;
    return end < other.end;
  }
};

#endif
