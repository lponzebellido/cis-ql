#ifndef MOTIF_FINDER_H
#define MOTIF_FINDER_H

#include <string>
#include <vector>

struct MotifMatch {
  size_t position;
  size_t matchLength;
  std::string strand;
  std::string context;
};

class MotifFinder {
public:
  static std::vector<MotifMatch> findAll(
      const std::string& sequence, const std::string& pattern,
      const std::string& chr = "", bool searchNegativeStrand = true);

  static std::vector<MotifMatch> findInWindow(
      const std::string& sequence, const std::string& pattern,
      size_t windowStart, size_t windowEnd, const std::string& chr = "");

  static std::string reverseComplement(const std::string& seq);
  static bool isRegexPattern(const std::string& pattern);

private:
  static std::vector<int> computeKMPTable(const std::string& pattern);
  static std::vector<size_t> kmpSearch(const std::string& text, const std::string& pattern);
  static std::vector<MotifMatch> regexSearch(const std::string& text, const std::string& pattern, const std::string& strand);
};

#endif
