#ifndef TRACK_EVIDENCE_H
#define TRACK_EVIDENCE_H

#include <cstddef>
#include <string>

struct TrackEvidence {
  bool present = false;
  std::string trackAlias;
  std::string source;
  std::string format;
  bool hasScore = false;
  double score = 0.0;
  bool hasSignalValue = false;
  double signalValue = 0.0;
  bool hasMinusLog10PValue = false;
  double minusLog10PValue = 0.0;
  bool hasMinusLog10QValue = false;
  double minusLog10QValue = 0.0;
  bool hasPeak = false;
  size_t peakOffset = 0;
};

#endif
