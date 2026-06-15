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
};

#endif
