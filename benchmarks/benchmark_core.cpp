#include "../src/bioinfo/GCAnalyzer.h"
#include "../src/bioinfo/MotifFinder.h"
#include "../src/bioinfo/PWMScanner.h"
#include "../src/bioinfo/RegulatoryRegions.h"
#include "../src/bioinfo/SetOperations.h"
#include "../src/bioinfo/SmithWaterman.h"

#include <chrono>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

std::string randomDna(size_t length, uint32_t seed) {
  static const char bases[] = {'A', 'C', 'G', 'T'};
  std::mt19937 generator(seed);
  std::uniform_int_distribution<int> distribution(0, 3);
  std::string result(length, 'A');
  for (char &base : result)
    base = bases[distribution(generator)];
  return result;
}

template <typename Function>
double medianElapsed(Function function, size_t repetitions = 5) {
  function();
  std::vector<double> samples;
  samples.reserve(repetitions);
  for (size_t repetition = 0; repetition < repetitions; ++repetition) {
    const auto start = Clock::now();
    function();
    samples.push_back(
        std::chrono::duration<double>(Clock::now() - start).count());
  }
  std::sort(samples.begin(), samples.end());
  return samples[samples.size() / 2];
}

GenomicRegion makeRegion(size_t start, size_t end) {
  GenomicRegion region;
  region.chr = "chr1";
  region.start = start;
  region.end = end;
  region.strand = "+";
  region.type = "benchmark";
  return region;
}

void printResult(const std::string &operation, size_t size, double seconds,
                 double resultValue, const std::string &resultUnit) {
  std::cout << operation << "," << size << "," << seconds << ","
            << resultValue << "," << resultUnit << "\n";
}

}


int main() {
  std::cout << std::setprecision(9);
  std::cout << "operation,input_size,median_seconds,result_value,result_unit\n";

  for (size_t size : {100000U, 1000000U, 5000000U}) {
    const std::string sequence = randomDna(size, 100 + size);
    size_t count = 0;
    const double seconds = medianElapsed([&] {
      count = MotifFinder::findAll(sequence, "ACGTACGT", "chr1", false).size();
    });
    printResult("kmp_exact", size, seconds, count, "hits");
  }

  for (size_t size : {100000U, 1000000U, 5000000U}) {
    const std::string sequence = randomDna(size, 200 + size);
    size_t count = 0;
    const double seconds = medianElapsed([&] {
      count = GCAnalyzer::findCpGIslands(sequence, "chr1").size();
    });
    printResult("cpg_islands", size, seconds, count, "islands");
  }

  PWMatrix matrix;
  matrix.id = "BENCH";
  matrix.name = "BENCH";
  matrix.length = 8;
  matrix.counts = {
      {10, 1, 1, 1, 10, 1, 1, 1},
      {1, 10, 1, 1, 1, 10, 1, 1},
      {1, 1, 10, 1, 1, 1, 10, 1},
      {1, 1, 1, 10, 1, 1, 1, 10},
  };
  const PSSM pssm = PWMScanner::computePSSM(matrix);
  for (size_t size : {100000U, 1000000U, 5000000U}) {
    const std::string sequence = randomDna(size, 300 + size);
    size_t count = 0;
    const double seconds = medianElapsed([&] {
      count =
          PWMScanner::scan(sequence, pssm, 80.0, "chr1", true, true).size();
    });
    printResult("pwm_scan", size, seconds, count, "hits");
  }

  for (size_t size : {250U, 500U, 1000U, 2000U, 10000U}) {
    const std::string first = randomDna(size, 400 + size);
    std::string second = first;
    for (size_t i = 0; i < second.size(); i += 20)
      second[i] = second[i] == 'A' ? 'C' : 'A';
    double score = 0.0;
    const double seconds = medianElapsed(
        [&] { score = SmithWaterman::computeSimilarity(first, second); });
    printResult("smith_waterman", size, seconds, score, "normalized_percent");
  }

  for (size_t size : {10000U, 100000U, 250000U}) {
    std::vector<GenomicRegion> first;
    std::vector<GenomicRegion> second;
    first.reserve(size);
    second.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      first.push_back(makeRegion(i * 20, i * 20 + 12));
      second.push_back(makeRegion(i * 20 + 8, i * 20 + 18));
    }
    size_t count = 0;
    const double seconds = medianElapsed(
        [&] { count = SetOperations::intersect(first, second).size(); });
    printResult("interval_intersect", size, seconds, count, "regions");

    size_t totalOverlaps = 0;
    const double countSeconds = medianElapsed([&] {
      totalOverlaps = 0;
      const auto counts = SetOperations::countOverlaps(
          first, second, "first", "second");
      for (const auto &container : counts)
        totalOverlaps += container.countEvidence.count;
    });
    printResult("interval_count_overlaps", size, countSeconds,
                totalOverlaps, "overlaps");

    const double overlapSeconds = medianElapsed(
        [&] { count = SetOperations::selectOverlapping(
                         first, second, "reference").size(); });
    printResult("interval_overlaps", size, overlapSeconds, count, "regions");

    const double nearSeconds = medianElapsed(
        [&] { count = SetOperations::selectNear(
                         first, second, 25, "reference").size(); });
    printResult("interval_near", size, nearSeconds, count, "regions");
  }

  for (size_t size : {10000U, 100000U, 250000U}) {
    std::vector<GenomicRegion> genes;
    genes.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      GenomicRegion gene = makeRegion(i * 40 + 100, i * 40 + 130);
      gene.name = "gene_" + std::to_string(i);
      if (i % 2 == 1)
        gene.strand = "-";
      genes.push_back(std::move(gene));
    }
    std::vector<GenomicRegion> promoters;
    std::string error;
    bool success = false;
    const std::unordered_map<std::string, size_t> chromosomeLengths = {
        {"chr1", size * 40 + 500}};
    const double seconds = medianElapsed([&] {
      success = RegulatoryRegions::buildPromoters(
          genes, chromosomeLengths, 60, 20, promoters, error);
    });
    printResult("promoter_windows", size, seconds,
                success ? promoters.size() : 0, "regions");
  }

  for (size_t size : {10000U, 50000U, 100000U}) {
    std::vector<GenomicRegion> first;
    std::vector<GenomicRegion> second;
    first.reserve(size);
    second.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      first.push_back(makeRegion(i * 100, i * 100 + 6));
      second.push_back(makeRegion(i * 100 + 16, i * 100 + 22));
    }
    size_t count = 0;
    const double seconds = medianElapsed([&] {
      count = SetOperations::defineModules(
                  first, second, 10, 10, "AS_WRITTEN", "ANY",
                  "first", "second").size();
    });
    printResult("motif_modules", size, seconds, count, "modules");
  }

  for (size_t size : {10000U, 50000U, 100000U}) {
    std::vector<GenomicRegion> anchor;
    std::vector<GenomicRegion> support;
    anchor.reserve(size);
    support.reserve(size);
    for (size_t i = 0; i < size; ++i) {
      anchor.push_back(makeRegion(i * 50, i * 50 + 20));
      support.push_back(makeRegion(i * 50 + 2, i * 50 + 22));
    }
    size_t count = 0;
    const double seconds = medianElapsed([&] {
      count = SetOperations::consensus(
                  anchor, "replicate_1", {{"replicate_2", support}},
                  {"replicate_1", "replicate_2"}, 2, 50.0).size();
    });
    printResult("track_consensus", size, seconds, count, "regions");
  }
  return 0;
}
