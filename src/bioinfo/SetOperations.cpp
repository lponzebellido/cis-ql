#include "SetOperations.h"
#include <algorithm>

std::vector<GenomicRegion> SetOperations::intersect(
    std::vector<GenomicRegion> a, std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  size_t j = 0;
  for (const auto& ra : a) {
    while (j < b.size() && (b[j].chr < ra.chr || (b[j].chr == ra.chr && b[j].end <= ra.start))) {
      j++;
    }

    bool foundOverlap = false;
    std::string matchedType = "";
    for (size_t k = j; k < b.size(); ++k) {
      if (b[k].chr > ra.chr || b[k].start >= ra.end) {
        break;
      }
      if (ra.overlaps(b[k])) {
        foundOverlap = true;
        matchedType = b[k].type;
        break;
      }
    }

    if (foundOverlap) {
      GenomicRegion overlappingA = ra;
      if (!matchedType.empty()) {
        overlappingA.type = ra.type + " ∩ " + matchedType;
      }
      result.push_back(overlappingA);
    }
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::unite(
    std::vector<GenomicRegion> a, std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> all;
  all.reserve(a.size() + b.size());
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

  size_t j = 0;
  for (const auto& ra : a) {
    while (j < b.size() && (b[j].chr < ra.chr || (b[j].chr == ra.chr && b[j].end <= ra.start))) {
      j++;
    }

    bool excluded = false;
    for (size_t k = j; k < b.size(); ++k) {
      if (b[k].chr > ra.chr || b[k].start >= ra.end) {
        break;
      }
      if (ra.overlaps(b[k])) {
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
