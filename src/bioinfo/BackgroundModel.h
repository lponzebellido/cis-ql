#ifndef BACKGROUND_MODEL_H
#define BACKGROUND_MODEL_H

#include <cstddef>
#include <string>

struct BackgroundModel {
  std::string mode = "uniform";
  std::string source = "default";
  double a = 0.25;
  double c = 0.25;
  double g = 0.25;
  double t = 0.25;
  double estimationPseudocount = 0.0;
  size_t observedBases = 0;
  std::string strandPolicy = "forward";
};

class BackgroundModelEstimator {
public:
  static BackgroundModel uniform(const std::string &source = "explicit");
};

class BackgroundModelAccumulator {
private:
  double counts[4];
  size_t observedBases;

public:
  BackgroundModelAccumulator();
  void addSequence(const std::string &sequence);
  void addRange(const std::string &sequence, size_t start, size_t end);
  bool build(const std::string &source, const std::string &strandPolicy,
             BackgroundModel &model, std::string &error,
             double totalPseudocount = 0.1) const;
};

#endif
