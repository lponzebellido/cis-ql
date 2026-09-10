#include "SetOperations.h"
#include <algorithm>

std::vector<GenomicRegion>
SetOperations::intersect(std::vector<GenomicRegion> a,
                         std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  size_t i = 0;
  size_t j = 0;
  while (i < a.size() && j < b.size()) {
    if (a[i].chr < b[j].chr) {
      ++i;
      continue;
    }
    if (b[j].chr < a[i].chr) {
      ++j;
      continue;
    }

    const size_t start = std::max(a[i].start, b[j].start);
    const size_t end = std::min(a[i].end, b[j].end);
    if (start < end) {
      GenomicRegion overlap = a[i];
      overlap.start = start;
      overlap.end = end;
      overlap.type = a[i].type + " ∩ " + b[j].type;
      overlap.name = a[i].name + "_intersect_" + b[j].name;
      if (!overlap.sequence.empty()) {
        const size_t offset = start - a[i].start;
        overlap.sequence = overlap.sequence.substr(offset, end - start);
      }
      if (start != a[i].start || end != a[i].end)
        overlap.motifEvidence = MotifEvidence();
      result.push_back(std::move(overlap));
    }

    if (a[i].end <= b[j].end)
      ++i;
    else
      ++j;
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::unite(std::vector<GenomicRegion> a,
                                                std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> all;
  all.reserve(a.size() + b.size());
  all.insert(all.end(), a.begin(), a.end());
  all.insert(all.end(), b.begin(), b.end());
  std::sort(all.begin(), all.end());

  if (all.empty())
    return all;

  std::vector<GenomicRegion> result;
  result.push_back(all[0]);

  for (size_t i = 1; i < all.size(); i++) {
    GenomicRegion &last = result.back();
    if (all[i].chr == last.chr && all[i].start <= last.end) {
      if (all[i].end > last.end) {
        last.end = all[i].end;
        last.sequence.clear();
      }
      last.type = "union";
      last.motifEvidence = MotifEvidence();
    } else {
      result.push_back(all[i]);
    }
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::except(std::vector<GenomicRegion> a,
                                                 std::vector<GenomicRegion> b) {
  std::vector<GenomicRegion> result;
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  size_t j = 0;
  for (const auto &region : a) {
    while (j < b.size() &&
           (b[j].chr < region.chr ||
            (b[j].chr == region.chr && b[j].end <= region.start))) {
      ++j;
    }

    size_t cursor = region.start;
    size_t k = j;
    while (k < b.size() && b[k].chr == region.chr && b[k].start < region.end) {
      if (b[k].end <= cursor) {
        ++k;
        continue;
      }
      if (b[k].start > cursor) {
        GenomicRegion fragment = region;
        fragment.start = cursor;
        fragment.end = std::min(b[k].start, region.end);
        fragment.name =
            region.name + "_except_" + std::to_string(fragment.start);
        if (!region.sequence.empty()) {
          fragment.sequence = region.sequence.substr(
              fragment.start - region.start, fragment.end - fragment.start);
        }
        fragment.motifEvidence = MotifEvidence();
        result.push_back(std::move(fragment));
      }
      cursor = std::max(cursor, b[k].end);
      if (cursor >= region.end)
        break;
      ++k;
    }

    if (cursor < region.end) {
      GenomicRegion fragment = region;
      fragment.start = cursor;
      fragment.name = region.name + "_except_" + std::to_string(fragment.start);
      if (!region.sequence.empty()) {
        fragment.sequence = region.sequence.substr(
            fragment.start - region.start, fragment.end - fragment.start);
      }
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.motifEvidence = MotifEvidence();
      result.push_back(std::move(fragment));
    }
  }
  return result;
}
