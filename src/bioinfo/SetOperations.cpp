#include "SetOperations.h"
#include <algorithm>
#include <limits>
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
      if (start != a[i].start || end != a[i].end)
        overlap.countEvidence = CountEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.moduleEvidence = ModuleEvidence();
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
      last.countEvidence = CountEvidence();
      last.moduleEvidence = ModuleEvidence();
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
        fragment.countEvidence = CountEvidence();
        fragment.moduleEvidence = ModuleEvidence();
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
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.countEvidence = CountEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.moduleEvidence = ModuleEvidence();
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

std::vector<GenomicRegion> SetOperations::countOverlaps(
    const std::vector<GenomicRegion> &counted,
    const std::vector<GenomicRegion> &containers,
    const std::string &countedSet, const std::string &containerSet) {
  struct ChromosomeIndex {
    std::vector<size_t> starts;
    std::vector<size_t> ends;
  };

  std::unordered_map<std::string, ChromosomeIndex> indexes;
  for (const auto &region : counted) {
    if (region.start >= region.end)
      continue;
    ChromosomeIndex &index = indexes[region.chr];
    index.starts.push_back(region.start);
    index.ends.push_back(region.end);
  }
  for (auto &entry : indexes) {
    std::sort(entry.second.starts.begin(), entry.second.starts.end());
    std::sort(entry.second.ends.begin(), entry.second.ends.end());
  }

  std::vector<GenomicRegion> result;
  result.reserve(containers.size());
  for (const auto &container : containers) {
    if (container.start >= container.end)
      continue;
    size_t overlapCount = 0;
    const auto chromosome = indexes.find(container.chr);
    if (chromosome != indexes.end()) {
      const ChromosomeIndex &index = chromosome->second;
      const size_t startedBeforeEnd = static_cast<size_t>(std::distance(
          index.starts.begin(),
          std::lower_bound(index.starts.begin(), index.starts.end(),
                           container.end)));
      const size_t endedBeforeOrAtStart = static_cast<size_t>(std::distance(
          index.ends.begin(),
          std::upper_bound(index.ends.begin(), index.ends.end(),
                           container.start)));
      overlapCount = startedBeforeEnd - endedBeforeOrAtStart;
    }

    GenomicRegion countedContainer = container;
    countedContainer.countEvidence.present = true;
    countedContainer.countEvidence.relation = "OVERLAPS";
    countedContainer.countEvidence.countedSet = countedSet;
    countedContainer.countEvidence.containerSet = containerSet;
    countedContainer.countEvidence.count = overlapCount;
    result.push_back(std::move(countedContainer));
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::defineModules(
    const std::vector<GenomicRegion> &first,
    const std::vector<GenomicRegion> &second, size_t minimumSpacing,
    size_t maximumSpacing, const std::string &orderPolicy,
    const std::string &orientationPolicy, const std::string &firstSet,
    const std::string &secondSet) {
  struct ChromosomeIndex {
    std::vector<const GenomicRegion *> regions;
    std::vector<size_t> starts;
    size_t maximumLength = 0;
  };

  const auto regionOrder = [](const GenomicRegion &left,
                              const GenomicRegion &right) {
    if (left.start != right.start)
      return left.start < right.start;
    if (left.end != right.end)
      return left.end < right.end;
    if (left.strand != right.strand)
      return left.strand < right.strand;
    if (left.name != right.name)
      return left.name < right.name;
    return left.type < right.type;
  };
  const auto pointerOrder = [&regionOrder](const GenomicRegion *left,
                                           const GenomicRegion *right) {
    return regionOrder(*left, *right);
  };

  std::unordered_map<std::string, ChromosomeIndex> indexes;
  for (const auto &region : second) {
    if (region.start >= region.end)
      continue;
    ChromosomeIndex &index = indexes[region.chr];
    index.regions.push_back(&region);
    index.maximumLength = std::max(index.maximumLength, region.length());
  }
  for (auto &entry : indexes) {
    ChromosomeIndex &index = entry.second;
    std::sort(index.regions.begin(), index.regions.end(), pointerOrder);
    for (const auto *region : index.regions)
      index.starts.push_back(region->start);
  }

  const auto fillMember = [](ModuleMemberEvidence &member,
                             const GenomicRegion &region,
                             const std::string &sourceSet) {
    member.sourceSet = sourceSet;
    member.chr = region.chr;
    member.start = region.start;
    member.end = region.end;
    member.strand = region.strand;
    member.type = region.type;
    member.name = region.name;
    member.motifEvidence = region.motifEvidence;
  };

  const bool sameSet = firstSet == secondSet;
  std::vector<GenomicRegion> result;
  for (const auto &firstRegion : first) {
    if (firstRegion.start >= firstRegion.end)
      continue;
    const auto chromosome = indexes.find(firstRegion.chr);
    if (chromosome == indexes.end())
      continue;
    const ChromosomeIndex &index = chromosome->second;

    size_t reach = maximumSpacing;
    if (reach > std::numeric_limits<size_t>::max() - index.maximumLength)
      reach = std::numeric_limits<size_t>::max();
    else
      reach += index.maximumLength;
    const size_t lowerStart = firstRegion.start > reach
                                  ? firstRegion.start - reach
                                  : 0;
    const size_t upperStart =
        firstRegion.end > std::numeric_limits<size_t>::max() - maximumSpacing
            ? std::numeric_limits<size_t>::max()
            : firstRegion.end + maximumSpacing;
    auto candidate = std::lower_bound(index.starts.begin(), index.starts.end(),
                                      lowerStart);
    const auto candidateEnd = std::upper_bound(
        index.starts.begin(), index.starts.end(), upperStart);

    for (; candidate != candidateEnd; ++candidate) {
      const size_t position = static_cast<size_t>(
          std::distance(index.starts.begin(), candidate));
      const GenomicRegion &secondRegion = *index.regions[position];

      if (sameSet) {
        if (firstRegion.start == secondRegion.start &&
            firstRegion.end == secondRegion.end)
          continue;
        if (!regionOrder(firstRegion, secondRegion))
          continue;
      }

      const bool overlaps = firstRegion.overlaps(secondRegion);
      size_t spacing = 0;
      if (firstRegion.end < secondRegion.start)
        spacing = secondRegion.start - firstRegion.end;
      else if (secondRegion.end < firstRegion.start)
        spacing = firstRegion.start - secondRegion.end;
      if (spacing < minimumSpacing || spacing > maximumSpacing)
        continue;

      const bool firstBeforeSecond = firstRegion.start < secondRegion.start;
      if (orderPolicy == "AS_WRITTEN" && !firstBeforeSecond)
        continue;

      std::string observedOrientation = "UNKNOWN";
      const bool firstStranded = firstRegion.strand == "+" ||
                                 firstRegion.strand == "-";
      const bool secondStranded = secondRegion.strand == "+" ||
                                  secondRegion.strand == "-";
      if (firstStranded && secondStranded) {
        observedOrientation = firstRegion.strand == secondRegion.strand
                                  ? "SAME"
                                  : "OPPOSITE";
      }
      if (orientationPolicy != "ANY" &&
          observedOrientation != orientationPolicy)
        continue;

      GenomicRegion module;
      module.chr = firstRegion.chr;
      module.start = std::min(firstRegion.start, secondRegion.start);
      module.end = std::max(firstRegion.end, secondRegion.end);
      module.strand = ".";
      module.type = "cis_regulatory_module";
      module.name = firstRegion.name + "__with__" + secondRegion.name;
      module.moduleEvidence.present = true;
      module.moduleEvidence.minimumSpacing = minimumSpacing;
      module.moduleEvidence.maximumSpacing = maximumSpacing;
      module.moduleEvidence.observedSpacing = spacing;
      module.moduleEvidence.orderPolicy = orderPolicy;
      module.moduleEvidence.observedOrder = overlaps
          ? "OVERLAPPING"
          : firstBeforeSecond ? "FIRST_BEFORE_SECOND"
                              : "SECOND_BEFORE_FIRST";
      module.moduleEvidence.orientationPolicy = orientationPolicy;
      module.moduleEvidence.observedOrientation = observedOrientation;
      fillMember(module.moduleEvidence.first, firstRegion, firstSet);
      fillMember(module.moduleEvidence.second, secondRegion, secondSet);
      result.push_back(std::move(module));
    }
  }

  std::sort(result.begin(), result.end());
  return result;
}
