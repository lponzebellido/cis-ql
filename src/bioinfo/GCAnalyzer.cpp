#include "GCAnalyzer.h"
#include <cctype>

std::vector<GCWindow> GCAnalyzer::gcContentWindowed(
    const std::string& sequence, size_t windowSize, size_t stepSize) {
  std::vector<GCWindow> results;
  if (windowSize == 0 || stepSize == 0 || sequence.size() < windowSize)
    return results;

  size_t gcCount = 0;
  size_t previousPos = 0;
  for (size_t pos = 0; pos + windowSize <= sequence.size(); pos += stepSize) {
    if (pos == 0 || stepSize >= windowSize) {
      gcCount = 0;
      for (size_t i = pos; i < pos + windowSize; ++i) {
        const char c = static_cast<char>(
            std::toupper(static_cast<unsigned char>(sequence[i])));
        if (c == 'G' || c == 'C')
          ++gcCount;
      }
    } else {
      for (size_t i = previousPos; i < pos; ++i) {
        const char c = static_cast<char>(
            std::toupper(static_cast<unsigned char>(sequence[i])));
        if (c == 'G' || c == 'C')
          --gcCount;
      }
      for (size_t i = previousPos + windowSize; i < pos + windowSize; ++i) {
        const char c = static_cast<char>(
            std::toupper(static_cast<unsigned char>(sequence[i])));
        if (c == 'G' || c == 'C')
          ++gcCount;
      }
    }
    results.push_back(
        {pos, static_cast<double>(gcCount) / windowSize * 100.0});
    previousPos = pos;
  }
  return results;
}

std::vector<GenomicRegion> GCAnalyzer::findCpGIslands(
    const std::string& sequence, const std::string& chrId,
    size_t minLength, double minGC, double minObsExpCpG) {
  std::vector<GenomicRegion> islands;
  const size_t scanWindow = 200;
  const size_t scanStep = 1;

  if (sequence.size() < scanWindow) return islands;

  size_t cCount = 0;
  size_t gCount = 0;
  size_t cpgCount = 0;
  for (size_t i = 0; i < scanWindow; ++i) {
    const char c = std::toupper(sequence[i]);
    if (c == 'C')
      ++cCount;
    else if (c == 'G')
      ++gCount;
    if (i + 1 < scanWindow && c == 'C' &&
        std::toupper(sequence[i + 1]) == 'G') {
      ++cpgCount;
    }
  }

  for (size_t pos = 0; pos + scanWindow <= sequence.size();
       pos += scanStep) {
    double gc = (double)(cCount + gCount) / scanWindow;
    double expectedCpG = (double)(cCount * gCount) / scanWindow;
    double obsExp = expectedCpG > 0 ? (double)cpgCount / expectedCpG : 0;

    bool passes = gc >= minGC && obsExp >= minObsExpCpG;

    if (passes) {
      const size_t candidateEnd = pos + scanWindow;
      if (islands.empty() || pos > islands.back().end) {
        GenomicRegion r;
        r.chr = chrId;
        r.start = pos;
        r.end = candidateEnd;
        r.strand = "+";
        r.type = "CpG_island";
        r.name = "CpG_" + std::to_string(pos);
        islands.push_back(std::move(r));
      } else {
        islands.back().end = std::max(islands.back().end, candidateEnd);
      }
    }

    if (pos + scanWindow >= sequence.size())
      break;

    const char outgoing = std::toupper(sequence[pos]);
    if (outgoing == 'C')
      --cCount;
    else if (outgoing == 'G')
      --gCount;
    if (outgoing == 'C' && std::toupper(sequence[pos + 1]) == 'G')
      --cpgCount;

    const size_t incomingPos = pos + scanWindow;
    const char incoming = std::toupper(sequence[incomingPos]);
    if (incoming == 'C')
      ++cCount;
    else if (incoming == 'G')
      ++gCount;
    if (std::toupper(sequence[incomingPos - 1]) == 'C' && incoming == 'G')
      ++cpgCount;
  }

  std::vector<GenomicRegion> filtered;
  for (auto &island : islands) {
    if (island.length() < minLength)
      continue;
    island.sequence =
        sequence.substr(island.start, island.end - island.start);
    filtered.push_back(std::move(island));
  }
  return filtered;
}
