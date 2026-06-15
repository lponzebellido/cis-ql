#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#include <string>

struct AlignmentResult {
  double score;
  double similarity;
  std::string alignedQuery;
  std::string alignedTarget;
  size_t queryStart;
  size_t targetStart;
};

class SmithWaterman {
public:
  static AlignmentResult align(const std::string& query, const std::string& target,
                               int matchScore = 2, int mismatchPenalty = -1,
                               int gapPenalty = -2);

  static double computeSimilarity(const std::string& seq1, const std::string& seq2);

private:
  static AlignmentResult alignBanded(const std::string& query, const std::string& target,
                                     int bandwidth, int matchScore, int mismatchPenalty,
                                     int gapPenalty);
};

#endif
