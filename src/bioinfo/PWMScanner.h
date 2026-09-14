#ifndef PWM_SCANNER_H
#define PWM_SCANNER_H

#include "BackgroundModel.h"
#include "MotifFinder.h"
#include <string>
#include <vector>




struct PWMatrix {
  std::string name;
  std::string id;
  std::string source;
  int length;                          
  std::vector<std::vector<double>> counts; 
};


struct PSSM {
  std::string name;
  std::string id;
  std::string source;
  int length;
  std::vector<std::vector<double>> scores; 
  double maxScore;  
  double minScore;  
  BackgroundModel background;
  double motifPseudocount;
  std::vector<std::vector<int>> scaledScores;
  std::vector<double> pValueByScaledScore;
  int scoreRange;
  double scoreScale;
  double scoreOffset;
};

struct PWMScanResult {
  std::vector<MotifMatch> matches;
  std::vector<size_t> testedScoreCounts;
  size_t testedPositions = 0;
};

class PWMScanner {
public:
  
  static PWMatrix loadJASPAR(const std::string& filename);

  
  
  static PSSM computePSSM(
      const PWMatrix &pwm,
      const BackgroundModel &background = BackgroundModel(),
      double motifPseudocount = 0.1);

  
  static std::vector<MotifMatch> scan(const std::string& sequence,
                                       const PSSM& pssm,
                                       double thresholdPercent,
                                       const std::string& chrId,
                                       bool searchPositive = true,
                                       bool searchNegative = true);

  static PWMScanResult scanWithStatistics(
      const std::string &sequence, const PSSM &pssm,
      double thresholdPercent, const std::string &chrId,
      bool searchPositive = true, bool searchNegative = true);

  static void mergeScanResults(PWMScanResult &destination,
                               PWMScanResult source);
  static void applyBenjaminiHochberg(PWMScanResult &result,
                                     const PSSM &pssm);

  
  static double scoreToPercent(double score, const PSSM& pssm);

private:
  
  static PWMScanResult scanStrand(const std::string& sequence,
                                  const PSSM& pssm,
                                  double minRawScore,
                                  const std::string& strand);
};

#endif
