#ifndef PATTERN_EVIDENCE_H
#define PATTERN_EVIDENCE_H

#include <cstddef>
#include <string>
#include <vector>

struct PatternCapture {
  size_t index = 0;
  bool matched = false;
  std::string value;
  size_t start = 0;
  size_t end = 0;
};

struct PatternEvidence {
  bool present = false;
  std::string pattern;
  std::vector<PatternCapture> captures;
};

#endif
