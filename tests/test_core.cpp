#include "../src/bioinfo/GCAnalyzer.h"
#include "../src/bioinfo/MotifFinder.h"
#include "../src/bioinfo/PWMScanner.h"
#include "../src/bioinfo/RegulatoryRegions.h"
#include "../src/bioinfo/SetOperations.h"
#include "../src/bioinfo/SmithWaterman.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAILED: " << message << std::endl;
    std::exit(1);
  }
}

GenomicRegion region(size_t start, size_t end,
                     const std::string &chr = "chr1") {
  GenomicRegion value;
  value.chr = chr;
  value.start = start;
  value.end = end;
  value.strand = "+";
  value.type = "test";
  value.name = "r" + std::to_string(start);
  return value;
}

}


int main() {
  const auto exact = MotifFinder::findAll("ACGTACGT", "ACG", "chr1", false);
  require(exact.size() == 2 && exact[0].position == 0 &&
              exact[1].position == 4,
          "KMP exact matching");

  const auto iupac =
      MotifFinder::findAll("TATAAATTATATAT", "TATAWAW", "chr1", false);
  require(iupac.size() == 2, "IUPAC matching");
  const auto overlappingIupac =
      MotifFinder::findAll("AAA", "NN", "chr1", false);
  require(overlappingIupac.size() == 2 &&
              overlappingIupac[0].position == 0 &&
              overlappingIupac[1].position == 1,
          "overlapping IUPAC matches");
  require(!MotifFinder::isValidPattern("("),
          "invalid regular expressions are rejected");

  const double identical =
      SmithWaterman::computeSimilarity("ACGTACGT", "ACGTACGT");
  require(std::abs(identical - 100.0) < 1e-9,
          "Smith-Waterman identical sequences");
  const double partial =
      SmithWaterman::computeSimilarity("ACGTACGT", "ACGTAAAA");
  require(partial >= 0.0 && partial <= 100.0,
          "Smith-Waterman normalized score range");
  require(std::abs(
              partial -
              SmithWaterman::align("ACGTACGT", "ACGTAAAA").similarity) <
              1e-9,
          "score-only Smith-Waterman matches the traceback implementation");
  std::mt19937 sequenceGenerator(42);
  const char bases[] = {'A', 'C', 'G', 'T'};
  std::string longQuery(10000, 'A');
  for (char &base : longQuery)
    base = bases[sequenceGenerator() % 4];
  std::string longTarget = longQuery;
  longTarget.erase(longTarget.begin() + 5000);
  require(SmithWaterman::computeSimilarity(longQuery, longTarget) > 99.9,
          "banded Smith-Waterman handles an internal deletion");

  const auto overlap =
      SetOperations::intersect({region(0, 10)}, {region(5, 15)});
  require(overlap.size() == 1 && overlap[0].start == 5 &&
              overlap[0].end == 10,
          "geometric interval intersection");

  const auto difference =
      SetOperations::except({region(0, 10)}, {region(3, 7)});
  require(difference.size() == 2 && difference[0].start == 0 &&
              difference[0].end == 3 && difference[1].start == 7 &&
              difference[1].end == 10,
          "geometric interval subtraction");
  const auto multipleOverlaps = SetOperations::intersect(
      {region(0, 20)}, {region(2, 5), region(8, 12), region(15, 25)});
  require(multipleOverlaps.size() == 3 &&
              multipleOverlaps[0].start == 2 &&
              multipleOverlaps[1].start == 8 &&
              multipleOverlaps[2].start == 15 &&
              multipleOverlaps[2].end == 20,
          "intersection preserves every overlap segment");
  const auto multipleDifference = SetOperations::except(
      {region(0, 20)}, {region(2, 5), region(8, 12), region(15, 25)});
  require(multipleDifference.size() == 3 &&
              multipleDifference[0].start == 0 &&
              multipleDifference[0].end == 2 &&
              multipleDifference[1].start == 5 &&
              multipleDifference[1].end == 8 &&
              multipleDifference[2].start == 12 &&
              multipleDifference[2].end == 15,
          "subtraction handles multiple overlapping intervals");

  GenomicRegion supportedSite = region(0, 10);
  supportedSite.sequence = "AAAAAAAAAA";
  supportedSite.motifEvidence.present = true;
  supportedSite.motifEvidence.matrixId = "TEST";
  const auto supportedSites = SetOperations::selectOverlapping(
      {supportedSite, region(20, 30), region(5, 10, "chr2")},
      {region(3, 6), region(4, 8), region(10, 20)}, "references");
  require(supportedSites.size() == 1 &&
              supportedSites[0].start == 0 && supportedSites[0].end == 10 &&
              supportedSites[0].sequence == "AAAAAAAAAA" &&
              supportedSites[0].motifEvidence.present &&
              supportedSites[0].motifEvidence.matrixId == "TEST" &&
              supportedSites[0].overlapEvidence.size() == 2 &&
              supportedSites[0].overlapEvidence[0].referenceSet ==
                  "references" &&
              supportedSites[0].overlapEvidence[0].referenceStart == 3 &&
              supportedSites[0].overlapEvidence[1].referenceStart == 4,
          "overlap semi-join preserves one query and every matching reference");

  GenomicRegion replicateOne = region(3, 6);
  replicateOne.name = "replicate_one_peak";
  replicateOne.trackEvidence.present = true;
  replicateOne.trackEvidence.replicate = "R1";
  OverlapEvidence replicateTwoSupport;
  replicateTwoSupport.referenceSet = "replicate_two";
  replicateTwoSupport.referenceName = "replicate_two_peak";
  replicateTwoSupport.trackEvidence.present = true;
  replicateTwoSupport.trackEvidence.replicate = "R2";
  replicateOne.overlapEvidence.push_back(replicateTwoSupport);
  const auto nestedSupport = SetOperations::selectOverlapping(
      {region(0, 10)}, {replicateOne}, "replicate_one");
  require(nestedSupport.size() == 1 &&
              nestedSupport[0].overlapEvidence.size() == 1 &&
              nestedSupport[0].overlapEvidence[0].trackEvidence.replicate ==
                  "R1" &&
              nestedSupport[0].overlapEvidence[0]
                      .supportingEvidence.size() == 1 &&
              nestedSupport[0].overlapEvidence[0]
                      .supportingEvidence[0].trackEvidence.replicate == "R2",
          "overlap semi-join preserves nested reference provenance");

  GenomicRegion anchorPeakOne = region(0, 10);
  anchorPeakOne.name = "anchor_one";
  GenomicRegion anchorPeakTwo = region(20, 30);
  anchorPeakTwo.name = "anchor_two";
  GenomicRegion replicateTwoFirst = region(2, 8);
  replicateTwoFirst.name = "replicate_two_first";
  GenomicRegion replicateTwoSecond = region(22, 28);
  replicateTwoSecond.name = "replicate_two_second";
  GenomicRegion replicateThreeFirst = region(4, 12);
  replicateThreeFirst.name = "replicate_three_first";
  GenomicRegion replicateThreeMarginal = region(29, 35);
  replicateThreeMarginal.name = "replicate_three_marginal";
  const auto permissiveConsensus = SetOperations::consensus(
      {anchorPeakOne, anchorPeakTwo}, "replicate_one",
      {{"replicate_two", {replicateTwoFirst, replicateTwoSecond}},
       {"replicate_three",
        {replicateThreeFirst, replicateThreeMarginal}}},
      {"replicate_one", "replicate_two", "replicate_three"}, 3);
  require(permissiveConsensus.size() == 2 &&
              !permissiveConsensus[1].consensusEvidence
                   .hasMinimumReciprocalOverlap,
          "consensus keeps one-base overlaps when no stricter rule is given");
  const auto strictConsensus = SetOperations::consensus(
      {anchorPeakOne, anchorPeakTwo}, "replicate_one",
      {{"replicate_two", {replicateTwoFirst, replicateTwoSecond}},
       {"replicate_three",
        {replicateThreeFirst, replicateThreeMarginal}}},
      {"replicate_one", "replicate_two", "replicate_three"}, 3, 50.0);
  require(strictConsensus.size() == 1 &&
              strictConsensus[0].name == "anchor_one" &&
              strictConsensus[0].consensusEvidence.present &&
              strictConsensus[0].consensusEvidence.minimumSupport == 3 &&
              strictConsensus[0].consensusEvidence.observedSupport == 3 &&
              strictConsensus[0].consensusEvidence
                      .hasMinimumReciprocalOverlap &&
              strictConsensus[0].consensusEvidence
                      .minimumReciprocalOverlapPercent == 50.0 &&
              strictConsensus[0].overlapEvidence.size() == 2 &&
              strictConsensus[0].overlapEvidence[0].referenceSet ==
                  "replicate_two" &&
              strictConsensus[0].overlapEvidence[1].referenceSet ==
                  "replicate_three",
          "consensus applies reciprocal overlap per distinct supporting set");

  GenomicRegion summitAnchorOne = region(0, 10);
  summitAnchorOne.name = "summit_anchor_one";
  summitAnchorOne.trackEvidence.present = true;
  summitAnchorOne.trackEvidence.hasPeak = true;
  summitAnchorOne.trackEvidence.peakPosition = 5;
  GenomicRegion summitAnchorTwo = region(20, 30);
  summitAnchorTwo.name = "summit_anchor_two";
  summitAnchorTwo.trackEvidence.present = true;
  summitAnchorTwo.trackEvidence.hasPeak = true;
  summitAnchorTwo.trackEvidence.peakPosition = 22;
  GenomicRegion summitSupportOne = region(1, 9);
  summitSupportOne.name = "summit_support_one";
  summitSupportOne.trackEvidence.present = true;
  summitSupportOne.trackEvidence.hasPeak = true;
  summitSupportOne.trackEvidence.peakPosition = 5;
  GenomicRegion summitSupportTwo = region(21, 29);
  summitSupportTwo.name = "summit_support_two";
  summitSupportTwo.trackEvidence.present = true;
  summitSupportTwo.trackEvidence.hasPeak = true;
  summitSupportTwo.trackEvidence.peakPosition = 28;
  const auto summitConsensus = SetOperations::consensus(
      {summitAnchorOne, summitAnchorTwo}, "summit_anchor",
      {{"summit_support", {summitSupportOne, summitSupportTwo}}},
      {"summit_anchor", "summit_support"}, 2, 50.0, true, 2);
  require(summitConsensus.size() == 1 &&
              summitConsensus[0].name == "summit_anchor_one" &&
              summitConsensus[0].consensusEvidence
                  .hasMaximumSummitDistance &&
              summitConsensus[0].consensusEvidence
                      .maximumSummitDistanceBp == 2,
          "consensus rejects support with distant or absent summits");

  const auto missingSummitConsensus = SetOperations::consensus(
      {anchorPeakOne}, "anchor_without_summit",
      {{"support_without_summit", {replicateTwoFirst}}},
      {"anchor_without_summit", "support_without_summit"}, 2, 0.0, true,
      2);
  require(missingSummitConsensus.empty(),
          "summit consensus requires summit-bearing intervals");

  std::vector<GenomicRegion> indexedReferences;
  std::vector<GenomicRegion> indexedQueries;
  for (size_t index = 0; index < 250; ++index) {
    const size_t start = sequenceGenerator() % 2000;
    GenomicRegion referenceRegion =
        region(start, start + 1 + sequenceGenerator() % 80,
               index % 3 == 0 ? "chr2" : "chr1");
    referenceRegion.name = "reference_" + std::to_string(index);
    indexedReferences.push_back(referenceRegion);
  }
  for (size_t index = 0; index < 200; ++index) {
    const size_t start = sequenceGenerator() % 2000;
    GenomicRegion queryRegion =
        region(start, start + 1 + sequenceGenerator() % 100,
               index % 4 == 0 ? "chr2" : "chr1");
    queryRegion.name = "query_" + std::to_string(index);
    indexedQueries.push_back(queryRegion);
  }
  const auto indexedOverlaps = SetOperations::selectOverlapping(
      indexedQueries, indexedReferences, "random_references");
  std::unordered_map<std::string, const GenomicRegion *> observedByName;
  for (const auto &queryRegion : indexedOverlaps)
    observedByName[queryRegion.name] = &queryRegion;
  const auto expectedOrder = [](const GenomicRegion *left,
                                const GenomicRegion *right) {
    if (left->start != right->start)
      return left->start < right->start;
    if (left->end != right->end)
      return left->end < right->end;
    return left->name < right->name;
  };
  for (const auto &queryRegion : indexedQueries) {
    std::vector<const GenomicRegion *> expected;
    for (const auto &referenceRegion : indexedReferences) {
      if (queryRegion.overlaps(referenceRegion))
        expected.push_back(&referenceRegion);
    }
    std::sort(expected.begin(), expected.end(), expectedOrder);
    const auto observed = observedByName.find(queryRegion.name);
    require((expected.empty() && observed == observedByName.end()) ||
                (!expected.empty() && observed != observedByName.end()),
            "interval-tree overlap support matches brute-force membership");
    if (expected.empty())
      continue;
    require(observed->second->overlapEvidence.size() == expected.size(),
            "interval-tree overlap support matches brute-force count");
    for (size_t index = 0; index < expected.size(); ++index) {
      require(observed->second->overlapEvidence[index].referenceName ==
                  expected[index]->name,
              "interval-tree overlap support has deterministic reference order");
    }
  }

  GenomicRegion nearQuery = region(10, 20);
  nearQuery.sequence = "AAAAAAAAAA";
  nearQuery.motifEvidence.present = true;
  nearQuery.motifEvidence.matrixId = "TEST";
  GenomicRegion leftReference = region(0, 5);
  leftReference.name = "left_gene";
  leftReference.type = "gene";
  GenomicRegion rightReference = region(25, 30);
  rightReference.name = "right_gene";
  rightReference.type = "gene";
  const auto nearest = SetOperations::selectNear(
      {nearQuery}, {rightReference, leftReference}, 5, "GENE");
  require(nearest.size() == 1 && nearest[0].start == 10 &&
              nearest[0].end == 20 && nearest[0].sequence == "AAAAAAAAAA" &&
              nearest[0].motifEvidence.present &&
              nearest[0].spatialRelation.present &&
              nearest[0].spatialRelation.relation == "NEAR" &&
              nearest[0].spatialRelation.referenceSet == "GENE" &&
              nearest[0].spatialRelation.referenceName == "left_gene" &&
              nearest[0].spatialRelation.distance == 5 &&
              nearest[0].spatialRelation.maximumDistance == 5 &&
              !nearest[0].spatialRelation.overlaps,
          "nearest selection preserves query evidence and breaks ties deterministically");

  GenomicRegion overlappingReference = region(0, 12);
  overlappingReference.name = "overlapping_gene";
  const auto overlappingNearest = SetOperations::selectNear(
      {nearQuery}, {overlappingReference, rightReference}, 0, "genes");
  require(overlappingNearest.size() == 1 &&
              overlappingNearest[0].spatialRelation.distance == 0 &&
              overlappingNearest[0].spatialRelation.overlaps &&
              overlappingNearest[0].spatialRelation.referenceName ==
                  "overlapping_gene",
          "NEAR records overlap as zero-distance spatial evidence");

  const auto adjacentNearest = SetOperations::selectNear(
      {region(10, 20)}, {region(20, 30)}, 0, "features");
  require(adjacentNearest.size() == 1 &&
              adjacentNearest[0].spatialRelation.distance == 0 &&
              !adjacentNearest[0].spatialRelation.overlaps,
          "NEAR distinguishes adjacent intervals from overlap at distance zero");
  require(SetOperations::selectNear(
              {nearQuery}, {rightReference}, 4, "GENE").empty(),
          "NEAR applies its inclusive maximum-distance boundary");

  GenomicRegion countContainer = region(0, 10);
  countContainer.name = "promoter_a";
  countContainer.sequence = "AAAAAAAAAA";
  const auto overlapCounts = SetOperations::countOverlaps(
      {region(0, 5), region(5, 10), region(10, 15),
       region(3, 7, "chr2")},
      {countContainer, region(10, 20), region(0, 10, "chr2"),
       region(30, 40)},
      "sites", "promoters");
  require(overlapCounts.size() == 4 &&
              overlapCounts[0].name == "promoter_a" &&
              overlapCounts[0].sequence == "AAAAAAAAAA" &&
              overlapCounts[0].countEvidence.present &&
              overlapCounts[0].countEvidence.relation == "OVERLAPS" &&
              overlapCounts[0].countEvidence.countedSet == "sites" &&
              overlapCounts[0].countEvidence.containerSet == "promoters" &&
              overlapCounts[0].countEvidence.count == 2 &&
              overlapCounts[1].countEvidence.count == 1 &&
              overlapCounts[2].countEvidence.count == 1 &&
              overlapCounts[3].countEvidence.count == 0,
          "COUNT reports half-open overlap counts for every container");

  GenomicRegion moduleSiteA = region(10, 14);
  moduleSiteA.name = "site_a";
  moduleSiteA.motifEvidence.present = true;
  moduleSiteA.motifEvidence.matrixId = "MYB";
  GenomicRegion moduleSiteB = region(20, 24);
  moduleSiteB.name = "site_b";
  moduleSiteB.motifEvidence.present = true;
  moduleSiteB.motifEvidence.matrixId = "MYB";
  GenomicRegion moduleSiteC = region(30, 34);
  moduleSiteC.name = "site_c";
  moduleSiteC.strand = "-";
  moduleSiteC.motifEvidence.present = true;
  moduleSiteC.motifEvidence.matrixId = "MYB";
  const auto sameOrientationModules = SetOperations::defineModules(
      {moduleSiteA, moduleSiteB, moduleSiteC},
      {moduleSiteA, moduleSiteB, moduleSiteC}, 6, 6, "ANY", "SAME",
      "myb_sites", "myb_sites");
  require(sameOrientationModules.size() == 1 &&
              sameOrientationModules[0].start == 10 &&
              sameOrientationModules[0].end == 24 &&
              sameOrientationModules[0].moduleEvidence.present &&
              sameOrientationModules[0].moduleEvidence.observedSpacing == 6 &&
              sameOrientationModules[0].moduleEvidence.observedOrder ==
                  "FIRST_BEFORE_SECOND" &&
              sameOrientationModules[0].moduleEvidence.observedOrientation ==
                  "SAME" &&
              sameOrientationModules[0].moduleEvidence.first.motifEvidence
                      .matrixId == "MYB" &&
              sameOrientationModules[0].moduleEvidence.second.name ==
                  "site_b",
          "motif modules preserve both members and apply spacing and orientation");
  const auto oppositeOrientationModules = SetOperations::defineModules(
      {moduleSiteA, moduleSiteB, moduleSiteC},
      {moduleSiteA, moduleSiteB, moduleSiteC}, 6, 6, "ANY", "OPPOSITE",
      "myb_sites", "myb_sites");
  require(oppositeOrientationModules.size() == 1 &&
              oppositeOrientationModules[0].start == 20 &&
              oppositeOrientationModules[0].end == 34,
          "homotypic modules emit one canonical pair and reject self-pairs");
  require(SetOperations::defineModules(
              {moduleSiteC}, {moduleSiteA}, 16, 16, "AS_WRITTEN", "ANY",
              "right_sites", "left_sites").empty(),
          "AS_WRITTEN requires the first input member to occur first");
  const auto reverseOrderModules = SetOperations::defineModules(
      {moduleSiteC}, {moduleSiteA}, 16, 16, "ANY", "ANY", "right_sites",
      "left_sites");
  require(reverseOrderModules.size() == 1 &&
              reverseOrderModules[0].moduleEvidence.observedOrder ==
                  "SECOND_BEFORE_FIRST",
          "ANY order records rather than hides reverse input order");

  GenomicRegion firstUnion = region(0, 10);
  firstUnion.sequence = "AAAAAAAAAA";
  const auto merged = SetOperations::unite({firstUnion}, {region(5, 15)});
  require(merged.size() == 1 && merged[0].start == 0 &&
              merged[0].end == 15 && merged[0].sequence.empty(),
          "merged intervals do not retain stale sequence data");

  GenomicRegion positiveGene = region(100, 200);
  positiveGene.name = "positive";
  GenomicRegion negativeGene = region(300, 400);
  negativeGene.name = "negative";
  negativeGene.strand = "-";
  std::vector<GenomicRegion> promoters;
  std::string promoterError;
  require(RegulatoryRegions::buildPromoters(
              {positiveGene, negativeGene}, {{"chr1", 500}}, 50, 10,
              promoters, promoterError),
          "TSS-relative promoter construction succeeds");
  require(promoters.size() == 2 && promoters[0].start == 50 &&
              promoters[0].end == 110 && promoters[0].strand == "+" &&
              promoters[1].start == 390 && promoters[1].end == 450 &&
              promoters[1].strand == "-",
          "promoter coordinates respect source strand");

  GenomicRegion chromosomeEdge = region(5, 25);
  chromosomeEdge.name = "edge";
  require(RegulatoryRegions::buildPromoters(
              {chromosomeEdge}, {{"chr1", 500}}, 50, 10, promoters,
              promoterError) &&
              promoters.size() == 1 && promoters[0].start == 0 &&
              promoters[0].end == 15,
          "promoters are clamped to chromosome bounds");

  GenomicRegion unstranded = region(100, 200);
  unstranded.strand = ".";
  require(!RegulatoryRegions::buildPromoters(
              {unstranded}, {{"chr1", 500}}, 50, 10, promoters,
              promoterError) &&
              promoterError.find("strand") != std::string::npos,
          "unstranded promoter sources are rejected");
  require(!RegulatoryRegions::buildPromoters(
              {positiveGene}, {{"chr2", 500}}, 50, 10, promoters,
              promoterError) &&
              promoterError.find("active FASTA") != std::string::npos,
          "promoter sources must match the active genome");

  const auto gc = GCAnalyzer::gcContentWindowed("GGCCAAAA", 4, 4);
  require(gc.size() == 2 && std::abs(gc[0].gcPercent - 100.0) < 1e-9 &&
              std::abs(gc[1].gcPercent) < 1e-9,
          "non-overlapping GC windows");
  const auto sparseGc =
      GCAnalyzer::gcContentWindowed("GGCCAAAAGGGG", 4, 6);
  require(sparseGc.size() == 2 &&
              std::abs(sparseGc[0].gcPercent - 100.0) < 1e-9 &&
              std::abs(sparseGc[1].gcPercent - 50.0) < 1e-9,
          "GC windows with steps larger than the window");
  require(GCAnalyzer::gcContentWindowed("ACGT", 2, 0).empty(),
          "zero GC step is rejected");

  const auto islands =
      GCAnalyzer::findCpGIslands(std::string(150, 'C') + std::string(150, 'G'),
                                 "chr1");
  require(islands.empty(), "CpG detector does not confuse GC-rich sequence");
  const auto cpg =
      GCAnalyzer::findCpGIslands(std::string(150, 'C'), "chr1");
  require(cpg.empty(), "CpG detector requires observed CpG dinucleotides");
  std::string repeatedCpG;
  for (int i = 0; i < 150; ++i)
    repeatedCpG += "CG";
  require(!GCAnalyzer::findCpGIslands(repeatedCpG, "chr1").empty(),
          "CpG island detection");

  BackgroundModel estimatedBackground;
  std::string backgroundError;
  BackgroundModelAccumulator backgroundAccumulator;
  backgroundAccumulator.addSequence("AAAACCGTNN");
  require(backgroundAccumulator.build("test-sequences", "forward",
                                      estimatedBackground, backgroundError),
          "zero-order background estimation succeeds");
  require(estimatedBackground.observedBases == 8 &&
              estimatedBackground.strandPolicy == "forward" &&
              std::abs(estimatedBackground.a - 4.025 / 8.1) < 1e-12 &&
              std::abs(estimatedBackground.c - 2.025 / 8.1) < 1e-12 &&
              std::abs(estimatedBackground.g - 1.025 / 8.1) < 1e-12 &&
              std::abs(estimatedBackground.t - 1.025 / 8.1) < 1e-12,
          "background estimation ignores ambiguity and applies a total "
          "pseudocount of 0.1");

  BackgroundModel symmetricBackground;
  require(backgroundAccumulator.build("test-sequences", "symmetric",
                                      symmetricBackground, backgroundError) &&
              std::abs(symmetricBackground.a - symmetricBackground.t) <
                  1e-12 &&
              std::abs(symmetricBackground.c - symmetricBackground.g) <
                  1e-12,
          "two-strand backgrounds average reverse-complement frequencies");
  BackgroundModel reverseBackground;
  require(backgroundAccumulator.build("test-sequences", "reverse_complement",
                                      reverseBackground, backgroundError) &&
              std::abs(reverseBackground.a - estimatedBackground.t) < 1e-12 &&
              std::abs(reverseBackground.c - estimatedBackground.g) < 1e-12 &&
              std::abs(reverseBackground.g - estimatedBackground.c) < 1e-12 &&
              std::abs(reverseBackground.t - estimatedBackground.a) < 1e-12,
          "negative-strand backgrounds use reverse-complement frequencies");
  BackgroundModelAccumulator ambiguousBackground;
  ambiguousBackground.addSequence("NNNN");
  require(!ambiguousBackground.build("ambiguous", "forward",
                                     estimatedBackground, backgroundError) &&
              backgroundError.find("no unambiguous DNA bases") !=
                  std::string::npos,
          "empty effective background sources are rejected");

  PWMatrix matrix;
  matrix.id = "TEST";
  matrix.name = "TEST";
  matrix.source = "test-matrix.pwm";
  matrix.length = 2;
  matrix.counts = {{10, 10}, {0, 0}, {0, 0}, {0, 0}};
  const PSSM pssm = PWMScanner::computePSSM(matrix);
  const double expectedBestScore =
      2.0 * std::log2(((10.0 + 0.1 * 0.25) / 10.1) / 0.25);
  require(std::abs(pssm.maxScore - expectedBestScore) < 1e-12,
          "PSSM uses the MEME/FIMO background-weighted motif pseudocount");
  const auto pwmHits = PWMScanner::scan("AAAA", pssm, 90.0, "chrPWM", true,
                                        false);
  require(pwmHits.size() == 3 && pwmHits[0].chr == "chrPWM",
          "PWM scanning preserves chromosome identifiers");
  require(pwmHits[0].evidence.present &&
              pwmHits[0].evidence.matrixId == "TEST" &&
              pwmHits[0].evidence.matrixName == "TEST" &&
              pwmHits[0].evidence.matrixSource == "test-matrix.pwm" &&
              pwmHits[0].evidence.background.mode == "uniform" &&
              pwmHits[0].evidence.motifPseudocount == 0.1 &&
              std::abs(pwmHits[0].evidence.rawScore - pssm.maxScore) < 1e-9 &&
              std::abs(pwmHits[0].evidence.scorePercent - 100.0) < 1e-9,
          "PWM matches retain matrix provenance and scores");

  require(pssm.scoreRange == 1000 && pssm.scoreScale > 0.0 &&
              pssm.pValueByScaledScore.size() == 2001,
          "PSSM builds a FIMO-range dynamic-programming score table");
  PWMScanResult positiveStatistics = PWMScanner::scanWithStatistics(
      "AAAA", pssm, 90.0, "chrPWM", true, false);
  PWMScanner::applyBenjaminiHochberg(positiveStatistics, pssm);
  require(positiveStatistics.testedPositions == 3 &&
              positiveStatistics.matches.size() == 3 &&
              std::abs(positiveStatistics.matches[0]
                           .evidence.statistics.pValue - 0.0625) < 1e-12 &&
              std::abs(positiveStatistics.matches[0]
                           .evidence.statistics.qValue - 0.0625) < 1e-12,
          "PWM p-values match the exact AA null and tied BH correction");
  PWMScanResult twoStrandStatistics = PWMScanner::scanWithStatistics(
      "AAAA", pssm, 90.0, "chrPWM", true, true);
  PWMScanner::applyBenjaminiHochberg(twoStrandStatistics, pssm);
  require(twoStrandStatistics.testedPositions == 6 &&
              twoStrandStatistics.matches.size() == 3 &&
              std::abs(twoStrandStatistics.matches[0]
                           .evidence.statistics.qValue - 0.125) < 1e-12,
          "BH correction covers every selected position-strand test");
  PWMScanResult ambiguousStatistics = PWMScanner::scanWithStatistics(
      "AANA", pssm, 0.0, "chrPWM", true, false);
  require(ambiguousStatistics.testedPositions == 1,
          "ambiguous PWM windows are excluded from the statistical universe");

  GenomicRegion scoredHit = region(0, 10);
  scoredHit.motifEvidence.present = true;
  scoredHit.motifEvidence.matrixId = "TEST";
  scoredHit.spatialRelation.present = true;
  scoredHit.countEvidence.present = true;
  const auto croppedHit =
      SetOperations::intersect({scoredHit}, {region(5, 15)});
  require(croppedHit.size() == 1 &&
              !croppedHit[0].motifEvidence.present &&
              !croppedHit[0].spatialRelation.present &&
              !croppedHit[0].countEvidence.present,
          "interval operations discard evidence when hit geometry changes");
  PWMatrix invalidMatrix = matrix;
  invalidMatrix.counts[3].pop_back();
  require(PWMScanner::computePSSM(invalidMatrix).length == 0,
          "PWM rows must have equal lengths");

  std::cout << "All core tests passed." << std::endl;
  return 0;
}
