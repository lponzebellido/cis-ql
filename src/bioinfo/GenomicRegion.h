#ifndef GENOMIC_REGION_H
#define GENOMIC_REGION_H

#include <string>
#include <algorithm>

struct GenomicRegion {
  std::string chr;
  size_t start;
  size_t end;
  std::string strand;
  std::string type;
  std::string name;
  std::string sequence;

  size_t length() const { return end - start; }

  bool overlaps(const GenomicRegion& other) const {
    return chr == other.chr && start < other.end && end > other.start;
  }

  bool operator<(const GenomicRegion& other) const {
    if (chr != other.chr) return chr < other.chr;
    if (start != other.start) return start < other.start;
    return end < other.end;
  }
};

#endif
