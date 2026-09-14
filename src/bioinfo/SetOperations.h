#ifndef SET_OPERATIONS_H
#define SET_OPERATIONS_H

#include "GenomicRegion.h"
#include <vector>

class SetOperations {
public:
  static std::vector<GenomicRegion> intersect(
      std::vector<GenomicRegion> a, std::vector<GenomicRegion> b);
  static std::vector<GenomicRegion> unite(
      std::vector<GenomicRegion> a, std::vector<GenomicRegion> b);
  static std::vector<GenomicRegion> except(
      std::vector<GenomicRegion> a, std::vector<GenomicRegion> b);
  static std::vector<GenomicRegion> selectOverlapping(
      const std::vector<GenomicRegion> &query,
      const std::vector<GenomicRegion> &reference);
  static std::vector<GenomicRegion> selectNear(
      const std::vector<GenomicRegion> &query,
      const std::vector<GenomicRegion> &reference,
      size_t maximumDistance, const std::string &referenceSet);
};

#endif
