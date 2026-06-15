#include "SetOperations.h"
#include <algorithm>

std::vector<GenomicRegion> SetOperations::intersect(
    std::vector<GenomicRegion> a, std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  for (const auto& ra : a) {
    for (const auto& rb : b) {
      if (ra.overlaps(rb)) {
        GenomicRegion overlap;
        overlap.chr = ra.chr;
        overlap.start = std::max(ra.start, rb.start);
        overlap.end = std::min(ra.end, rb.end);
        overlap.strand = ra.strand;
        overlap.type = ra.type + " ∩ " + rb.type;
        overlap.name = ra.name + " ∩ " + rb.name;
        result.push_back(overlap);
      }
    }
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::unite(
    std::vector<GenomicRegion> a, std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> all;
  all.insert(all.end(), a.begin(), a.end());
  all.insert(all.end(), b.begin(), b.end());
  std::sort(all.begin(), all.end());

  if (all.empty()) return all;

  std::vector<GenomicRegion> result;
  result.push_back(all[0]);

  for (size_t i = 1; i < all.size(); i++) {
    GenomicRegion& last = result.back();
    if (all[i].chr == last.chr && all[i].start <= last.end) {
      last.end = std::max(last.end, all[i].end);
      last.type = "union";
    } else {
      result.push_back(all[i]);
    }
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::except(
    std::vector<GenomicRegion> a, std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  for (const auto& ra : a) {
    bool excluded = false;
    for (const auto& rb : b) {
      if (ra.overlaps(rb)) {
        excluded = true;
        break;
      }
    }
    if (!excluded) {
      result.push_back(ra);
    }
  }
  return result;
}
