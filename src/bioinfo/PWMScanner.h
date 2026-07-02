#ifndef PWM_SCANNER_H
#define PWM_SCANNER_H

#include "MotifFinder.h"
#include <string>
#include <vector>




struct PWMatrix {
  std::string name;
  std::string id;
  int length;                          
  std::vector<std::vector<double>> counts; 
};


struct PSSM {
  std::string name;
  int length;
  std::vector<std::vector<double>> scores; 
  double maxScore;  
  double minScore;  
};

class PWMScanner {
public:
  
  static PWMatrix loadJASPAR(const std::string& filename);

  
  
  static PSSM computePSSM(const PWMatrix& pwm,
                           double bgA = 0.25, double bgC = 0.25,
                           double bgG = 0.25, double bgT = 0.25);

  
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
