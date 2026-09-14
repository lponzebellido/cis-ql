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
  static std::vector<GenomicRegion> countOverlaps(
      const std::vector<GenomicRegion> &counted,
      const std::vector<GenomicRegion> &containers,
      const std::string &countedSet, const std::string &containerSet);
  static std::vector<GenomicRegion> defineModules(
      const std::vector<GenomicRegion> &first,
      const std::vector<GenomicRegion> &second, size_t minimumSpacing,
      size_t maximumSpacing, const std::string &orderPolicy,
      const std::string &orientationPolicy, const std::string &firstSet,
      const std::string &secondSet);
};

#endif
