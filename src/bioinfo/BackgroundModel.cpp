#include "BackgroundModel.h"

#include <cctype>

BackgroundModel
BackgroundModelEstimator::uniform(const std::string &source) {
  BackgroundModel model;
  model.source = source;
  return model;
}

BackgroundModelAccumulator::BackgroundModelAccumulator()
    : counts{0.0, 0.0, 0.0, 0.0}, observedBases(0) {}

void BackgroundModelAccumulator::addSequence(const std::string &sequence) {
  addRange(sequence, 0, sequence.size());
}

void BackgroundModelAccumulator::addRange(const std::string &sequence,
                                          size_t start, size_t end) {
  if (start > end || end > sequence.size())
    return;
  for (size_t position = start; position < end; ++position) {
    const char nucleotide = sequence[position];
    switch (std::toupper(static_cast<unsigned char>(nucleotide))) {
    case 'A': counts[0] += 1.0; break;
    case 'C': counts[1] += 1.0; break;
    case 'G': counts[2] += 1.0; break;
    case 'T': counts[3] += 1.0; break;
    default: continue;
    }
    ++observedBases;
  }
}

bool BackgroundModelAccumulator::build(
    const std::string &source, const std::string &strandPolicy,
    BackgroundModel &model, std::string &error,
    double totalPseudocount) const {
  if (totalPseudocount < 0.0) {
    error = "Background pseudocount cannot be negative.";
    return false;
  }

  if (observedBases == 0) {
    error = "Background source contains no unambiguous DNA bases.";
    return false;
  }

  if (strandPolicy != "forward" &&
      strandPolicy != "reverse_complement" &&
      strandPolicy != "symmetric") {
    error = "Unknown background strand policy '" + strandPolicy + "'.";
    return false;
  }

  double effectiveCounts[4] = {counts[0], counts[1], counts[2], counts[3]};
  if (strandPolicy == "symmetric") {
    const double at = (effectiveCounts[0] + effectiveCounts[3]) / 2.0;
    const double cg = (effectiveCounts[1] + effectiveCounts[2]) / 2.0;
    effectiveCounts[0] = effectiveCounts[3] = at;
    effectiveCounts[1] = effectiveCounts[2] = cg;
  } else if (strandPolicy == "reverse_complement") {
    const double a = effectiveCounts[0];
    const double c = effectiveCounts[1];
    effectiveCounts[0] = effectiveCounts[3];
    effectiveCounts[1] = effectiveCounts[2];
    effectiveCounts[2] = c;
    effectiveCounts[3] = a;
  }

  const double perBasePseudocount = totalPseudocount / 4.0;
  const double denominator =
      static_cast<double>(observedBases) + totalPseudocount;

  model.mode = "estimated_zero_order";
  model.source = source;
  model.a = (effectiveCounts[0] + perBasePseudocount) / denominator;
  model.c = (effectiveCounts[1] + perBasePseudocount) / denominator;
  model.g = (effectiveCounts[2] + perBasePseudocount) / denominator;
  model.t = (effectiveCounts[3] + perBasePseudocount) / denominator;
  model.estimationPseudocount = totalPseudocount;
  model.observedBases = observedBases;
  model.strandPolicy = strandPolicy;
  error.clear();
  return true;
}
