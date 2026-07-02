#ifndef PWM_SCANNER_H
#define PWM_SCANNER_H

#include "MotifFinder.h"
#include <string>
#include <vector>

// Position Weight Matrix: stores raw frequency counts per position
// Rows: A=0, C=1, G=2, T=3
// Columns: one per position in the motif
struct PWMatrix {
  std::string name;
  std::string id;
  int length;                          // motif width (number of columns)
  std::vector<std::vector<double>> counts; // [4][length] raw counts
};

// Position-Specific Scoring Matrix: log-odds transformed
struct PSSM {
  std::string name;
  int length;
  std::vector<std::vector<double>> scores; // [4][length] log-odds scores
  double maxScore;  // theoretical max score (all best nucleotides)
  double minScore;  // theoretical min score (all worst nucleotides)
};

class PWMScanner {
public:
  // Load a PWM from JASPAR format file
  static PWMatrix loadJASPAR(const std::string& filename);

  // Convert frequency matrix to log-odds PSSM
  // bgFreqs: background nucleotide frequencies (default: uniform 0.25 each)
  static PSSM computePSSM(const PWMatrix& pwm,
                           double bgA = 0.25, double bgC = 0.25,
                           double bgG = 0.25, double bgT = 0.25);

  // Scan a sequence with a PSSM, returning matches above threshold (0-100%)
  static std::vector<MotifMatch> scan(const std::string& sequence,
                                       const PSSM& pssm,
                                       double thresholdPercent,
                                       const std::string& chrId,
                                       bool searchPositive = true,
                                       bool searchNegative = true);

  // Convert a raw score to percentage (0-100%) relative to PSSM min/max
  static double scoreToPercent(double score, const PSSM& pssm);

private:
  // Scan one strand, returning matches
  static std::vector<MotifMatch> scanStrand(const std::string& sequence,
                                             const PSSM& pssm,
                                             double minRawScore,
                                             const std::string& strand);
};

#endif
