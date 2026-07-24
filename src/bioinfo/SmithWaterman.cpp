#include "SmithWaterman.h"
#include <algorithm>
#include <cctype>
#include <vector>

AlignmentResult SmithWaterman::align(const std::string &query,
                                     const std::string &target, int matchScore,
                                     int mismatchPenalty, int gapPenalty) {
  size_t totalCells = (size_t)(query.size() + 1) * (size_t)(target.size() + 1);
  if (totalCells > 100000000) {
    int bandwidth =
        std::max(50, (int)std::min(query.size(), target.size()) / 10);
    return alignBanded(query, target, bandwidth, matchScore, mismatchPenalty,
                       gapPenalty);
  }

  int m = (int)query.size();
  int n = (int)target.size();
  std::vector<std::vector<int>> H(m + 1, std::vector<int>(n + 1, 0));

  int maxScore = 0;
  int maxI = 0, maxJ = 0;

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= n; j++) {
      int match = H[i - 1][j - 1] +
                  (std::toupper(query[i - 1]) == std::toupper(target[j - 1])
                       ? matchScore
                       : mismatchPenalty);
      int deleteGap = H[i - 1][j] + gapPenalty;
      int insertGap = H[i][j - 1] + gapPenalty;
      H[i][j] = std::max({0, match, deleteGap, insertGap});
      if (H[i][j] > maxScore) {
        maxScore = H[i][j];
        maxI = i;
        maxJ = j;
      }
    }
  }

  std::string alignedQuery, alignedTarget;
  int i = maxI, j = maxJ;
  while (i > 0 && j > 0 && H[i][j] > 0) {
    if (H[i][j] == H[i - 1][j - 1] + (std::toupper(query[i - 1]) ==
                                              std::toupper(target[j - 1])
                                          ? matchScore
                                          : mismatchPenalty)) {
      alignedQuery = query[i - 1] + alignedQuery;
      alignedTarget = target[j - 1] + alignedTarget;
      i--;
      j--;
    } else if (H[i][j] == H[i - 1][j] + gapPenalty) {
      alignedQuery = query[i - 1] + alignedQuery;
      alignedTarget = '-' + alignedTarget;
      i--;
    } else {
      alignedQuery = '-' + alignedQuery;
      alignedTarget = target[j - 1] + alignedTarget;
      j--;
    }
  }

  AlignmentResult result;
  result.score = maxScore;
  const size_t normalizer = std::min(query.size(), target.size());
  result.similarity =
      normalizer == 0
          ? 0.0
          : static_cast<double>(maxScore) /
                (static_cast<double>(normalizer) * matchScore) * 100.0;
  result.similarity = std::max(0.0, std::min(result.similarity, 100.0));
  result.alignedQuery = alignedQuery;
  result.alignedTarget = alignedTarget;
  result.queryStart = i;
  result.targetStart = j;
  return result;
}

AlignmentResult SmithWaterman::alignBanded(const std::string &query,
                                           const std::string &target,
                                           int bandwidth, int matchScore,
                                           int mismatchPenalty,
                                           int gapPenalty) {
  int m = (int)query.size();
  int n = (int)target.size();
  int W = bandwidth;

  std::vector<std::vector<int>> H(m + 1, std::vector<int>(2 * W + 1, 0));

  int maxScore = 0;
  for (int i = 1; i <= m; i++) {
    for (int b = 0; b < 2 * W + 1; b++) {
      int j = i - W + b;
      if (j < 1 || j > n)
        continue;

      int diagB = b;
      // H(i-1, j) is one diagonal to the right in the previous row because
      // b encodes j = i - W + b.
      int upB = b + 1;
      int leftB = b - 1;

      int match = 0;
      if (diagB >= 0 && diagB < 2 * W + 1) {
        match = H[i - 1][diagB] +
                (std::toupper(query[i - 1]) == std::toupper(target[j - 1])
                     ? matchScore
                     : mismatchPenalty);
      }
      int deleteGap = 0;
      if (upB >= 0 && upB < 2 * W + 1) {
        deleteGap = H[i - 1][upB] + gapPenalty;
      }
      int insertGap = 0;
      if (leftB >= 0 && leftB < 2 * W + 1) {
        insertGap = H[i][leftB] + gapPenalty;
      }

      H[i][b] = std::max({0, match, deleteGap, insertGap});
      if (H[i][b] > maxScore) {
        maxScore = H[i][b];
      }
    }
  }

  AlignmentResult result;
  result.score = maxScore;
  int alignLen = std::min(m, n);
  result.similarity =
      alignLen > 0 ? (double)maxScore / (alignLen * matchScore) * 100.0 : 0.0;
  result.similarity = std::min(result.similarity, 100.0);
  result.queryStart = 0;
  result.targetStart = 0;
  return result;
}

double SmithWaterman::computeSimilarity(const std::string &seq1,
                                        const std::string &seq2) {
  if (seq1.empty() || seq2.empty())
    return 0.0;

  const size_t totalCells = (seq1.size() + 1) * (seq2.size() + 1);
  if (totalCells > 100000000) {
    const int bandwidth =
        std::max(50, static_cast<int>(std::min(seq1.size(), seq2.size()) / 10));
    return alignBanded(seq1, seq2, bandwidth, 2, -1, -2).similarity;
  }

  // Similarity filtering only needs the maximum local-alignment score, not a
  // traceback. Two rows preserve the exact recurrence while reducing memory
  // from O(m*n) to O(min(m,n)).
  const std::string *rows = &seq1;
  const std::string *columns = &seq2;
  if (columns->size() > rows->size())
    std::swap(rows, columns);

  std::vector<int> previous(columns->size() + 1, 0);
  std::vector<int> current(columns->size() + 1, 0);
  int maxScore = 0;
  for (size_t i = 1; i <= rows->size(); ++i) {
    current[0] = 0;
    for (size_t j = 1; j <= columns->size(); ++j) {
      const int diagonal =
          previous[j - 1] +
          (std::toupper(static_cast<unsigned char>((*rows)[i - 1])) ==
                   std::toupper(static_cast<unsigned char>((*columns)[j - 1]))
               ? 2
               : -1);
      current[j] =
          std::max({0, diagonal, previous[j] - 2, current[j - 1] - 2});
      maxScore = std::max(maxScore, current[j]);
    }
    std::swap(previous, current);
  }

  const double maximum =
      static_cast<double>(std::min(seq1.size(), seq2.size())) * 2.0;
  return std::max(0.0, std::min(100.0, maxScore / maximum * 100.0));
}
