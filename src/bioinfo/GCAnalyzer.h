#ifndef GC_ANALYZER_H
#define GC_ANALYZER_H

#include "GenomicRegion.h"
#include <string>
#include <vector>

struct GCWindow {
  size_t position;
  double gcPercent;
  std::string chr;
};

class GCAnalyzer {
public:
  static std::vector<GCWindow> gcContentWindowed(
      const std::string& sequence, size_t windowSize = 100,
      size_t stepSize = 50);

  static std::vector<GenomicRegion> findCpGIslands(
      const std::string& sequence, const std::string& chrId,
      size_t minLength = 200, double minGC = 0.50,
      double minObsExpCpG = 0.6);
};

#endif
