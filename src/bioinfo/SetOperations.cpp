#include "SetOperations.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
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
        overlap.annotationEvidence = AnnotationEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.motifEvidence = MotifEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.spatialRelation = SpatialRelationEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.countEvidence = CountEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.moduleEvidence = ModuleEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.trackEvidence = TrackEvidence();
      if (start != a[i].start || end != a[i].end)
        overlap.overlapEvidence.clear();
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
      last.annotationEvidence = AnnotationEvidence();
      last.motifEvidence = MotifEvidence();
      last.spatialRelation = SpatialRelationEvidence();
      last.countEvidence = CountEvidence();
      last.moduleEvidence = ModuleEvidence();
      last.trackEvidence = TrackEvidence();
      last.overlapEvidence.clear();
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
        fragment.annotationEvidence = AnnotationEvidence();
        fragment.spatialRelation = SpatialRelationEvidence();
        fragment.countEvidence = CountEvidence();
        fragment.moduleEvidence = ModuleEvidence();
        fragment.trackEvidence = TrackEvidence();
        fragment.overlapEvidence.clear();
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
        fragment.annotationEvidence = AnnotationEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.motifEvidence = MotifEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.spatialRelation = SpatialRelationEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.countEvidence = CountEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.moduleEvidence = ModuleEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.trackEvidence = TrackEvidence();
      if (fragment.start != region.start || fragment.end != region.end)
        fragment.overlapEvidence.clear();
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
    const std::vector<GenomicRegion> &reference,
    const std::string &referenceSet,
    double minimumReciprocalOverlapPercent,
    bool hasMaximumSummitDistance,
    size_t maximumSummitDistanceBp) {
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

  struct IntervalNode {
    size_t center = 0;
    std::vector<const GenomicRegion *> spanningByStart;
    std::vector<const GenomicRegion *> spanningByEnd;
    std::unique_ptr<IntervalNode> left;
    std::unique_ptr<IntervalNode> right;
  };
  std::function<std::unique_ptr<IntervalNode>(
      std::vector<const GenomicRegion *>)>
      buildIndex;
  buildIndex = [&buildIndex, &referenceOrder](
                   std::vector<const GenomicRegion *> regions) {
    if (regions.empty())
      return std::unique_ptr<IntervalNode>();
    std::sort(regions.begin(), regions.end(), referenceOrder);
    std::unique_ptr<IntervalNode> node(new IntervalNode());
    node->center = regions[regions.size() / 2]->start;
    std::vector<const GenomicRegion *> left;
    std::vector<const GenomicRegion *> right;
    for (const auto *region : regions) {
      if (region->end <= node->center)
        left.push_back(region);
      else if (region->start > node->center)
        right.push_back(region);
      else
        node->spanningByStart.push_back(region);
    }
    std::sort(node->spanningByStart.begin(), node->spanningByStart.end(),
              referenceOrder);
    node->spanningByEnd = node->spanningByStart;
    std::sort(node->spanningByEnd.begin(), node->spanningByEnd.end(),
              [&referenceOrder](const GenomicRegion *leftRegion,
                                const GenomicRegion *rightRegion) {
                if (leftRegion->end != rightRegion->end)
                  return leftRegion->end > rightRegion->end;
                return referenceOrder(leftRegion, rightRegion);
              });
    node->left = buildIndex(std::move(left));
    node->right = buildIndex(std::move(right));
    return node;
  };

  std::unordered_map<std::string, std::vector<const GenomicRegion *>> grouped;
  for (const auto &region : reference) {
    if (region.start < region.end)
      grouped[region.chr].push_back(&region);
  }
  std::unordered_map<std::string, std::unique_ptr<IntervalNode>> indexes;
  for (auto &entry : grouped)
    indexes[entry.first] = buildIndex(std::move(entry.second));

  std::function<void(const IntervalNode *, size_t, size_t,
                     std::vector<const GenomicRegion *> &)>
      collectOverlaps;
  collectOverlaps = [&collectOverlaps](
                        const IntervalNode *node, size_t queryStart,
                        size_t queryEnd,
                        std::vector<const GenomicRegion *> &matches) {
    if (!node)
      return;
    if (queryEnd <= node->center) {
      for (const auto *region : node->spanningByStart) {
        if (region->start >= queryEnd)
          break;
        matches.push_back(region);
      }
      collectOverlaps(node->left.get(), queryStart, queryEnd, matches);
      return;
    }
    if (queryStart >= node->center) {
      for (const auto *region : node->spanningByEnd) {
        if (region->end <= queryStart)
          break;
        matches.push_back(region);
      }
      collectOverlaps(node->right.get(), queryStart, queryEnd, matches);
      return;
    }
    matches.insert(matches.end(), node->spanningByStart.begin(),
                   node->spanningByStart.end());
    collectOverlaps(node->left.get(), queryStart, queryEnd, matches);
    collectOverlaps(node->right.get(), queryStart, queryEnd, matches);
  };

  std::vector<GenomicRegion> result;
  for (const auto &region : query) {
    if (region.start >= region.end)
      continue;
    const auto chromosome = indexes.find(region.chr);
    if (chromosome == indexes.end())
      continue;
    std::vector<const GenomicRegion *> matches;
    collectOverlaps(chromosome->second.get(), region.start, region.end,
                    matches);
    if (minimumReciprocalOverlapPercent > 0.0 ||
        hasMaximumSummitDistance) {
      const double minimumFraction =
          minimumReciprocalOverlapPercent / 100.0;
      matches.erase(
          std::remove_if(
              matches.begin(), matches.end(),
              [&region, minimumFraction, hasMaximumSummitDistance,
               maximumSummitDistanceBp](const GenomicRegion *candidate) {
                const size_t overlapStart =
                    std::max(region.start, candidate->start);
                const size_t overlapEnd = std::min(region.end, candidate->end);
                const double overlapLength =
                    static_cast<double>(overlapEnd - overlapStart);
                const double queryFraction =
                    overlapLength / static_cast<double>(region.length());
                const double referenceFraction =
                    overlapLength / static_cast<double>(candidate->length());
                if (queryFraction < minimumFraction ||
                    referenceFraction < minimumFraction)
                  return true;
                if (!hasMaximumSummitDistance)
                  return false;
                if (!region.trackEvidence.hasPeak ||
                    !candidate->trackEvidence.hasPeak)
                  return true;
                const size_t querySummit =
                    region.trackEvidence.peakPosition;
                const size_t referenceSummit =
                    candidate->trackEvidence.peakPosition;
                const size_t summitDistance =
                    querySummit > referenceSummit
                        ? querySummit - referenceSummit
                        : referenceSummit - querySummit;
                return summitDistance > maximumSummitDistanceBp;
              }),
          matches.end());
    }
    std::sort(matches.begin(), matches.end(), referenceOrder);
    GenomicRegion selected = region;
    for (const auto *candidate : matches) {
      const GenomicRegion &referenceRegion = *candidate;
      OverlapEvidence evidence;
      evidence.referenceSet = referenceSet;
      evidence.referenceChr = referenceRegion.chr;
      evidence.referenceStart = referenceRegion.start;
      evidence.referenceEnd = referenceRegion.end;
      evidence.referenceStrand = referenceRegion.strand;
      evidence.referenceType = referenceRegion.type;
      evidence.referenceName = referenceRegion.name;
      evidence.trackEvidence = referenceRegion.trackEvidence;
      evidence.consensusEvidence = referenceRegion.consensusEvidence;
      evidence.supportingEvidence = referenceRegion.overlapEvidence;
      selected.overlapEvidence.push_back(std::move(evidence));
    }
    if (!matches.empty())
      result.push_back(std::move(selected));
  }
  return result;
}

std::vector<GenomicRegion> SetOperations::consensus(
    const std::vector<GenomicRegion> &anchor, const std::string &anchorSet,
    const std::vector<std::pair<std::string, std::vector<GenomicRegion>>>
        &supportSets,
    const std::vector<std::string> &inputSets, size_t minimumSupport,
    double minimumReciprocalOverlapPercent,
    bool hasMaximumSummitDistance,
    size_t maximumSummitDistanceBp) {
  std::vector<GenomicRegion> candidates = anchor;
  std::vector<size_t> observedSupport(anchor.size(), 1);

  const auto sameAnchor = [](const GenomicRegion &left,
                             const GenomicRegion &right) {
    return left.chr == right.chr && left.start == right.start &&
           left.end == right.end && left.strand == right.strand &&
           left.type == right.type && left.name == right.name;
  };

  for (const auto &supportSet : supportSets) {
    const std::vector<GenomicRegion> supported = selectOverlapping(
        anchor, supportSet.second, supportSet.first,
        minimumReciprocalOverlapPercent, hasMaximumSummitDistance,
        maximumSummitDistanceBp);
    size_t supportedIndex = 0;
    for (size_t anchorIndex = 0;
         anchorIndex < anchor.size() && supportedIndex < supported.size();
         ++anchorIndex) {
      if (!sameAnchor(anchor[anchorIndex], supported[supportedIndex]))
        continue;
      ++observedSupport[anchorIndex];
      const size_t inheritedEvidence =
          anchor[anchorIndex].overlapEvidence.size();
      candidates[anchorIndex].overlapEvidence.insert(
          candidates[anchorIndex].overlapEvidence.end(),
          supported[supportedIndex].overlapEvidence.begin() +
              inheritedEvidence,
          supported[supportedIndex].overlapEvidence.end());
      ++supportedIndex;
    }
  }

  std::vector<GenomicRegion> result;
  for (size_t index = 0; index < candidates.size(); ++index) {
    if (observedSupport[index] < minimumSupport)
      continue;
    candidates[index].consensusEvidence.present = true;
    candidates[index].consensusEvidence.anchorSet = anchorSet;
    candidates[index].consensusEvidence.minimumSupport = minimumSupport;
    candidates[index].consensusEvidence.observedSupport =
        observedSupport[index];
    if (minimumReciprocalOverlapPercent > 0.0) {
      candidates[index].consensusEvidence.hasMinimumReciprocalOverlap = true;
      candidates[index].consensusEvidence.minimumReciprocalOverlapPercent =
          minimumReciprocalOverlapPercent;
    }
    if (hasMaximumSummitDistance) {
      candidates[index].consensusEvidence.hasMaximumSummitDistance = true;
      candidates[index].consensusEvidence.maximumSummitDistanceBp =
          maximumSummitDistanceBp;
    }
    candidates[index].consensusEvidence.inputSets = inputSets;
    result.push_back(std::move(candidates[index]));
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
