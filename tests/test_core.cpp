#include "../src/bioinfo/GCAnalyzer.h"
#include "../src/bioinfo/MotifFinder.h"
#include "../src/bioinfo/PWMScanner.h"
#include "../src/bioinfo/SetOperations.h"
#include "../src/bioinfo/SmithWaterman.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
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

} // namespace

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

  GenomicRegion firstUnion = region(0, 10);
  firstUnion.sequence = "AAAAAAAAAA";
  const auto merged = SetOperations::unite({firstUnion}, {region(5, 15)});
  require(merged.size() == 1 && merged[0].start == 0 &&
              merged[0].end == 15 && merged[0].sequence.empty(),
          "merged intervals do not retain stale sequence data");

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

  PWMatrix matrix;
  matrix.id = "TEST";
  matrix.name = "TEST";
  matrix.length = 2;
  matrix.counts = {{10, 10}, {0, 0}, {0, 0}, {0, 0}};
  const PSSM pssm = PWMScanner::computePSSM(matrix);
  const auto pwmHits = PWMScanner::scan("AAAA", pssm, 90.0, "chrPWM", true,
                                        false);
  require(pwmHits.size() == 3 && pwmHits[0].chr == "chrPWM",
          "PWM scanning preserves chromosome identifiers");
  PWMatrix invalidMatrix = matrix;
  invalidMatrix.counts[3].pop_back();
  require(PWMScanner::computePSSM(invalidMatrix).length == 0,
          "PWM rows must have equal lengths");

  std::cout << "All core tests passed." << std::endl;
  return 0;
}
