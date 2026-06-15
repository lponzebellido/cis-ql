#ifndef GFF_READER_H
#define GFF_READER_H

#include "GenomicRegion.h"
#include <string>
#include <vector>
#include <unordered_map>

class GFFReader {
public:
  static std::vector<GenomicRegion> read(const std::string& filename);
  static std::vector<GenomicRegion> filterByType(
      const std::vector<GenomicRegion>& regions, const std::string& type);
  static std::string mapEntityToGFFType(const std::string& entity);
};

#endif
