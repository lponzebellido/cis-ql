#include "SetOperations.h"
#include <algorithm>
#include <unordered_map>
#include <utility>

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
      if (start != a[i].start || end != a[i].end)
        overlap.spatialRelation = SpatialRelationEvidence();
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
      last.spatialRelation = SpatialRelationEvidence();
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
        fragment.spatialRelation = SpatialRelationEvidence();
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
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.spatialRelation = SpatialRelationEvidence();
      result.push_back(std::move(fragment));
    }
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::selectNear(
    const std::vector<GenomicRegion> &query,
    const std::vector<GenomicRegion> &reference,
    size_t maximumDistance, const std::string &referenceSet) {
  struct ChromosomeIndex {
    std::vector<const GenomicRegion *> regions;
    std::vector<size_t> starts;
    std::vector<size_t> maximumEnds;
    std::vector<size_t> maximumEndIndexes;
  };

  std::unordered_map<std::string, ChromosomeIndex> indexes;
  for (const auto &region : reference) {
    if (region.start < region.end)
      indexes[region.chr].regions.push_back(&region);
  }

  const auto referenceOrder = [](const GenomicRegion *left,
                                 const GenomicRegion *right) {
    if (left->start != right->start)
      return left->start < right->start;
    if (left->end != right->end)
      return left->end < right->end;
    if (left->name != right->name)
      return left->name < right->name;
    if (left->type != right->type)
      return left->type < right->type;
    return left->strand < right->strand;
  };

  for (auto &entry : indexes) {
    ChromosomeIndex &index = entry.second;
    std::sort(index.regions.begin(), index.regions.end(), referenceOrder);
    size_t maximumEnd = 0;
    size_t maximumEndIndex = 0;
    for (size_t position = 0; position < index.regions.size(); ++position) {
      const GenomicRegion &region = *index.regions[position];
      index.starts.push_back(region.start);
      if (position == 0 || region.end > maximumEnd) {
        maximumEnd = region.end;
        maximumEndIndex = position;
      }
      index.maximumEnds.push_back(maximumEnd);
      index.maximumEndIndexes.push_back(maximumEndIndex);
    }
  }

  std::vector<GenomicRegion> result;
  for (const auto &region : query) {
    if (region.start >= region.end)
      continue;
    const auto chromosome = indexes.find(region.chr);
    if (chromosome == indexes.end())
      continue;
    const ChromosomeIndex &index = chromosome->second;

    const auto firstAtOrAfterEnd =
        std::lower_bound(index.starts.begin(), index.starts.end(), region.end);
    const size_t beforeEnd = static_cast<size_t>(
        std::distance(index.starts.begin(), firstAtOrAfterEnd));

    const GenomicRegion *nearest = nullptr;
    size_t nearestDistance = 0;
    if (beforeEnd > 0) {
      const auto firstOverlap = std::upper_bound(
          index.maximumEnds.begin(), index.maximumEnds.begin() + beforeEnd,
          region.start);
      if (firstOverlap != index.maximumEnds.begin() + beforeEnd) {
        const size_t overlapIndex = static_cast<size_t>(
            std::distance(index.maximumEnds.begin(), firstOverlap));
        nearest = index.regions[overlapIndex];
      } else {
        const size_t leftIndex = index.maximumEndIndexes[beforeEnd - 1];
        nearest = index.regions[leftIndex];
        nearestDistance = region.start - nearest->end;
      }
    }

    if (beforeEnd < index.regions.size()) {
      const GenomicRegion *right = index.regions[beforeEnd];
      const size_t rightDistance = right->start - region.end;
      if (!nearest || rightDistance < nearestDistance ||
          (rightDistance == nearestDistance &&
           referenceOrder(right, nearest))) {
        nearest = right;
        nearestDistance = rightDistance;
      }
    }

    if (!nearest || nearestDistance > maximumDistance)
      continue;

    GenomicRegion selected = region;
    selected.spatialRelation.present = true;
    selected.spatialRelation.relation = "NEAR";
    selected.spatialRelation.referenceSet = referenceSet;
    selected.spatialRelation.referenceChr = nearest->chr;
    selected.spatialRelation.referenceStart = nearest->start;
    selected.spatialRelation.referenceEnd = nearest->end;
    selected.spatialRelation.referenceStrand = nearest->strand;
    selected.spatialRelation.referenceType = nearest->type;
    selected.spatialRelation.referenceName = nearest->name;
    selected.spatialRelation.distance = nearestDistance;
    selected.spatialRelation.maximumDistance = maximumDistance;
    selected.spatialRelation.overlaps = region.overlaps(*nearest);
    result.push_back(std::move(selected));
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::selectOverlapping(
    const std::vector<GenomicRegion> &query,
    const std::vector<GenomicRegion> &reference) {
  struct ChromosomeIndex {
    std::vector<size_t> starts;
    std::vector<size_t> maximumEnds;
  };

  typedef std::pair<size_t, size_t> Interval;
  std::unordered_map<std::string, std::vector<Interval>> byChromosome;
  for (const auto &region : reference)
    byChromosome[region.chr].push_back(
        std::make_pair(region.start, region.end));

  std::unordered_map<std::string, ChromosomeIndex> indexes;
  for (auto &entry : byChromosome) {
    std::sort(entry.second.begin(), entry.second.end());
    ChromosomeIndex &index = indexes[entry.first];
    size_t maximumEnd = 0;
    for (const auto &interval : entry.second) {
      index.starts.push_back(interval.first);
      maximumEnd = std::max(maximumEnd, interval.second);
      index.maximumEnds.push_back(maximumEnd);
    }
  }

  std::vector<GenomicRegion> result;
  for (const auto &region : query) {
    if (region.start >= region.end)
      continue;
    const auto chromosome = indexes.find(region.chr);
    if (chromosome == indexes.end())
      continue;
    const ChromosomeIndex &index = chromosome->second;
    const auto firstTooLate =
        std::lower_bound(index.starts.begin(), index.starts.end(), region.end);
    const size_t candidateCount = static_cast<size_t>(
        std::distance(index.starts.begin(), firstTooLate));
    if (candidateCount > 0 &&
        index.maximumEnds[candidateCount - 1] > region.start) {
      result.push_back(region);
    }
  }
  return result;
}
