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

  
  static double scoreToPercent(double score, const PSSM& pssm);

private:
  
  static std::vector<MotifMatch> scanStrand(const std::string& sequence,
                                             const PSSM& pssm,
                                             double minRawScore,
                                             const std::string& strand);
};

#endif
