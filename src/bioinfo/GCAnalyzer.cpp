#include "GCAnalyzer.h"
#include <cctype>

std::vector<GCWindow> GCAnalyzer::gcContentWindowed(
    const std::string& sequence, size_t windowSize, size_t stepSize) {
  std::vector<GCWindow> results;
  if (sequence.size() < windowSize) return results;

  size_t gcCount = 0;
  for (size_t i = 0; i < windowSize; ++i) {
    char c = std::toupper(sequence[i]);
    if (c == 'G' || c == 'C') gcCount++;
  }
  results.push_back({0, (double)gcCount / windowSize * 100.0});

  for (size_t pos = stepSize; pos + windowSize <= sequence.size(); pos += stepSize) {
    for (size_t i = pos - stepSize; i < pos; ++i) {
      char c = std::toupper(sequence[i]);
      if (c == 'G' || c == 'C') gcCount--;
    }
    for (size_t i = pos + windowSize - stepSize; i < pos + windowSize; ++i) {
      char c = std::toupper(sequence[i]);
      if (c == 'G' || c == 'C') gcCount++;
    }
    results.push_back({pos, (double)gcCount / windowSize * 100.0});
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

  bool inIsland = false;
  size_t islandStart = 0;

  for (size_t pos = 0; pos + scanWindow <= sequence.size(); pos += scanStep) {
    size_t cCount = 0, gCount = 0, cpgCount = 0;
    for (size_t i = pos; i < pos + scanWindow; ++i) {
      char c = std::toupper(sequence[i]);
      if (c == 'C') {
        cCount++;
        if (i + 1 < pos + scanWindow && std::toupper(sequence[i + 1]) == 'G') {
          cpgCount++;
        }
      } else if (c == 'G') {
        gCount++;
      }
    }

    double gc = (double)(cCount + gCount) / scanWindow;
    double expectedCpG = (double)(cCount * gCount) / scanWindow;
    double obsExp = expectedCpG > 0 ? (double)cpgCount / expectedCpG : 0;

    bool passes = gc >= minGC && obsExp >= minObsExpCpG;

    if (passes && !inIsland) {
      inIsland = true;
      islandStart = pos;
    } else if (!passes && inIsland) {
      inIsland = false;
      size_t islandEnd = pos + scanWindow - 1;
      if (islandEnd - islandStart >= minLength) {
        GenomicRegion r;
        r.chr = chrId;
        r.start = islandStart;
        r.end = islandEnd;
        r.strand = "+";
        r.type = "CpG_island";
        r.name = "CpG_" + std::to_string(islandStart);
        if (islandEnd <= sequence.size()) {
          r.sequence = sequence.substr(islandStart, islandEnd - islandStart);
        }
        islands.push_back(r);
      }
    }
  }

  if (inIsland) {
    size_t islandEnd = sequence.size();
    if (islandEnd - islandStart >= minLength) {
      GenomicRegion r;
      r.chr = chrId;
      r.start = islandStart;
      r.end = islandEnd;
      r.strand = "+";
      r.type = "CpG_island";
      r.name = "CpG_" + std::to_string(islandStart);
      if (islandEnd <= sequence.size()) {
        r.sequence = sequence.substr(islandStart, islandEnd - islandStart);
      }
      islands.push_back(r);
    }
  }

  return islands;
}
