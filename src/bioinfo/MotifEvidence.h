#ifndef MOTIF_EVIDENCE_H
#define MOTIF_EVIDENCE_H

#include "BackgroundModel.h"
#include "TrackEvidence.h"
#include <cstddef>
#include <string>

struct MotifStatisticalEvidence {
  double pValue = 1.0;
  double qValue = 1.0;
  size_t testedPositions = 0;
  std::string pValueMethod = "zero_order_dynamic_programming";
  std::string multipleTestingMethod = "benjamini_hochberg";
  int scaledScore = 0;
  int scoreRange = 1000;
  double scoreScale = 0.0;
  double scoreOffset = 0.0;
};

struct MotifEvidence {
  bool present = false;
  std::string matrixAlias;
  std::string matrixId;
  std::string matrixName;
  std::string matrixSource;
  double rawScore = 0.0;
  double scorePercent = 0.0;
  BackgroundModel background;
  double motifPseudocount = 0.1;
  MotifStatisticalEvidence statistics;

  bool hasSourceRegion = false;
  std::string sourceRegionName;
  std::string sourceRegionType;
  size_t sourceRegionStart = 0;
  size_t sourceRegionEnd = 0;
  size_t relativeStart = 0;
  TrackEvidence sourceTrackEvidence;
};

#endif
