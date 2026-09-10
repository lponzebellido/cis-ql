#ifndef REGULATORY_REGIONS_H
#define REGULATORY_REGIONS_H

#include "GenomicRegion.h"
#include <string>
#include <unordered_map>
#include <vector>

class RegulatoryRegions {
public:
  static bool buildPromoters(
      const std::vector<GenomicRegion> &sources,
      const std::unordered_map<std::string, size_t> &chromosomeLengths,
      size_t upstream, size_t downstream,
      std::vector<GenomicRegion> &promoters, std::string &error);
};

#endif
