#include "../src/bioinfo/GCAnalyzer.h"
#include "../src/bioinfo/MotifFinder.h"
#include "../src/bioinfo/PWMScanner.h"
#include "../src/bioinfo/SetOperations.h"
#include "../src/bioinfo/SmithWaterman.h"

#include <chrono>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
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
  function(); // Warm up instruction and data paths before measuring.
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

} // namespace

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

  for (size_t size : {10000U, 100000U, 500000U}) {
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
  }
  return 0;
}
