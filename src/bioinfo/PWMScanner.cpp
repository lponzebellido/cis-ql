#include "PWMScanner.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>

namespace {

const int FIMO_SCORE_RANGE = 1000;

bool buildScoreDistribution(PSSM &pssm) {
  double smallest = 1e300;
  double largest = -1e300;
  for (int nucleotide = 0; nucleotide < 4; ++nucleotide) {
    for (int position = 0; position < pssm.length; ++position) {
      const double score = pssm.scores[nucleotide][position];
      if (!std::isfinite(score))
        return false;
      smallest = std::min(smallest, score);
      largest = std::max(largest, score);
    }
  }

  if (largest == smallest)
    smallest = largest - 1.0;
  pssm.scoreRange = FIMO_SCORE_RANGE;
  pssm.scoreOffset = std::floor(smallest);
  pssm.scoreScale = std::floor(
      static_cast<double>(pssm.scoreRange) /
      (largest - pssm.scoreOffset));
  if (!std::isfinite(pssm.scoreScale) || pssm.scoreScale <= 0.0)
    return false;

  pssm.scaledScores.assign(4, std::vector<int>(pssm.length, 0));
  for (int nucleotide = 0; nucleotide < 4; ++nucleotide) {
    for (int position = 0; position < pssm.length; ++position) {
      const double scaled =
          (pssm.scores[nucleotide][position] - pssm.scoreOffset) *
          pssm.scoreScale;
      const int rounded = static_cast<int>(std::floor(scaled + 0.5));
      if (rounded < 0 || rounded > pssm.scoreRange)
        return false;
      pssm.scaledScores[nucleotide][position] = rounded;
    }
  }

  const size_t distributionSize =
      static_cast<size_t>(pssm.length * pssm.scoreRange + 1);
  std::vector<double> probability(distributionSize, 0.0);
  std::vector<double> updated(distributionSize, 0.0);
  probability[0] = 1.0;
  size_t reachableMaximum = 0;
  const double background[4] = {
      pssm.background.a, pssm.background.c,
      pssm.background.g, pssm.background.t};

  for (int position = 0; position < pssm.length; ++position) {
    std::fill(updated.begin(), updated.end(), 0.0);
    size_t columnMaximum = 0;
    for (int nucleotide = 0; nucleotide < 4; ++nucleotide) {
      columnMaximum = std::max(
          columnMaximum,
          static_cast<size_t>(pssm.scaledScores[nucleotide][position]));
      const size_t shift = static_cast<size_t>(
          pssm.scaledScores[nucleotide][position]);
      for (size_t score = 0; score <= reachableMaximum; ++score) {
        if (probability[score] != 0.0) {
          updated[score + shift] +=
              probability[score] * background[nucleotide];
        }
      }
    }
    reachableMaximum += columnMaximum;
    probability.swap(updated);
  }

  pssm.pValueByScaledScore.assign(distributionSize, 0.0);
  double tailProbability = 0.0;
  for (size_t score = distributionSize; score-- > 0;) {
    tailProbability = std::min(1.0, tailProbability + probability[score]);
    pssm.pValueByScaledScore[score] = tailProbability;
  }
  return true;
}

} // namespace

static int nucToIndex(char c) {
  switch (c) {
  case 'A':
  case 'a':
    return 0;
  case 'C':
  case 'c':
    return 1;
  case 'G':
  case 'g':
    return 2;
  case 'T':
  case 't':
    return 3;
  default:
    return -1;
  }
}

static std::vector<double> parseLine(const std::string &line) {
  std::vector<double> values;
  std::string cleaned;

  for (size_t i = 0; i < line.size(); i++) {
    char c = line[i];
    if (c == '[' || c == ']' || c == '>' || c == '<') {
      cleaned += ' ';
    } else {
      cleaned += c;
    }
  }

  std::istringstream iss(cleaned);
  std::string token;
  while (iss >> token) {

    if (token.size() == 1 &&
        (token[0] == 'A' || token[0] == 'C' || token[0] == 'G' ||
         token[0] == 'T' || token[0] == 'a' || token[0] == 'c' ||
         token[0] == 'g' || token[0] == 't')) {
      continue;
    }
    try {
      double val = std::stod(token);
      values.push_back(val);
    } catch (...) {
    }
  }
  return values;
}

PWMatrix PWMScanner::loadJASPAR(const std::string &filename) {
  PWMatrix pwm;
  pwm.length = 0;
  pwm.source = filename;

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "  Error: Cannot open PWM file: " << filename << std::endl;
    return pwm;
  }

  std::string line;
  int rowIndex = 0;
  pwm.counts.resize(4);

  while (std::getline(file, line)) {

    size_t start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos)
      continue;
    line = line.substr(start);

    if (line[0] == '>') {

      std::istringstream iss(line.substr(1));
      iss >> pwm.id;
      if (iss >> pwm.name) {

      } else {
        pwm.name = pwm.id;
      }
      continue;
    }

    char firstAlpha = 0;
    for (char c : line) {
      if (std::isalpha(c)) {
        firstAlpha = std::toupper(c);
        break;
      }
    }

    int idx = -1;
    if (firstAlpha == 'A')
      idx = 0;
    else if (firstAlpha == 'C')
      idx = 1;
    else if (firstAlpha == 'G')
      idx = 2;
    else if (firstAlpha == 'T')
      idx = 3;
    else
      idx = rowIndex;

    std::vector<double> values = parseLine(line);
    if (!values.empty()) {
      if (idx < 0 || idx >= 4 || !pwm.counts[idx].empty()) {
        std::cerr << "  Error: Invalid or duplicate PWM row in " << filename
                  << std::endl;
        pwm.length = 0;
        return pwm;
      }
      pwm.counts[idx] = values;
      if ((int)values.size() > pwm.length) {
        pwm.length = (int)values.size();
      }
      rowIndex++;
    }
  }

  if (pwm.length == 0 || rowIndex < 4) {
    std::cerr << "  Error: Invalid PWM format in " << filename << std::endl;
    pwm.length = 0;
  }

  if (pwm.name.empty()) {
    pwm.name = filename;
  }

  return pwm;
}

PSSM PWMScanner::computePSSM(const PWMatrix &pwm,
                             const BackgroundModel &background,
                             double motifPseudocount) {
  PSSM pssm;
  pssm.name = pwm.name;
  pssm.id = pwm.id;
  pssm.source = pwm.source;
  pssm.length = pwm.length;
  pssm.background = background;
  pssm.motifPseudocount = motifPseudocount;
  pssm.scoreRange = FIMO_SCORE_RANGE;
  pssm.scoreScale = 0.0;
  pssm.scoreOffset = 0.0;
  const double backgroundSum =
      background.a + background.c + background.g + background.t;
  if (pwm.length <= 0 || pwm.counts.size() != 4 || background.a <= 0.0 ||
      background.c <= 0.0 || background.g <= 0.0 || background.t <= 0.0 ||
      !std::isfinite(backgroundSum) ||
      std::abs(backgroundSum - 1.0) > 1e-6 ||
      !std::isfinite(motifPseudocount) || motifPseudocount < 0.0) {
    pssm.length = 0;
    return pssm;
  }
  pssm.scores.resize(4);
  pssm.maxScore = 0.0;
  pssm.minScore = 0.0;

  double bg[4] = {background.a, background.c, background.g, background.t};

  for (int i = 0; i < 4; i++) {
    pssm.scores[i].resize(pwm.length, 0.0);
  }

  for (int j = 0; j < pwm.length; j++) {
    for (int i = 0; i < 4; ++i) {
      if (pwm.counts[i].size() != static_cast<size_t>(pwm.length)) {
        pssm.length = 0;
        pssm.scores.clear();
        return pssm;
      }
    }

    double colTotal = 0.0;
    for (int i = 0; i < 4; i++) {
      if (pwm.counts[i][j] < 0.0) {
        pssm.length = 0;
        pssm.scores.clear();
        return pssm;
      }
      colTotal += pwm.counts[i][j];
    }
    if (colTotal + motifPseudocount <= 0.0) {
      pssm.length = 0;
      pssm.scores.clear();
      return pssm;
    }

    double colMax = -1e9;
    double colMin = 1e9;
    for (int i = 0; i < 4; i++) {
      const double adjustedCount =
          pwm.counts[i][j] + motifPseudocount * bg[i];
      const double freq = adjustedCount / (colTotal + motifPseudocount);
      if (freq <= 0.0 || !std::isfinite(freq)) {
        pssm.length = 0;
        pssm.scores.clear();
        return pssm;
      }
      pssm.scores[i][j] = std::log2(freq / bg[i]);

      if (pssm.scores[i][j] > colMax)
        colMax = pssm.scores[i][j];
      if (pssm.scores[i][j] < colMin)
        colMin = pssm.scores[i][j];
    }
    pssm.maxScore += colMax;
    pssm.minScore += colMin;
  }

  if (!buildScoreDistribution(pssm)) {
    pssm.length = 0;
    pssm.scores.clear();
    pssm.scaledScores.clear();
    pssm.pValueByScaledScore.clear();
  }

  return pssm;
}

double PWMScanner::scoreToPercent(double score, const PSSM &pssm) {
  double range = pssm.maxScore - pssm.minScore;
  if (range <= 0.0)
    return 0.0;
  double pct = (score - pssm.minScore) / range * 100.0;
  return std::max(0.0, std::min(100.0, pct));
}

PWMScanResult PWMScanner::scanStrand(const std::string &sequence,
                                     const PSSM &pssm,
                                     double minRawScore,
                                     const std::string &strand,
                                     double maximumPValueForRetention) {

  PWMScanResult result;
  result.testedScoreCounts.assign(pssm.pValueByScaledScore.size(), 0);
  int seqLen = (int)sequence.size();
  int motifLen = pssm.length;

  if (seqLen < motifLen)
    return result;

  size_t contextSize = 20;

  for (int pos = 0; pos <= seqLen - motifLen; pos++) {
    double score = 0.0;
    int scaledScore = 0;
    bool valid = true;

    for (int j = 0; j < motifLen; j++) {
      int idx = nucToIndex(sequence[pos + j]);
      if (idx < 0) {
        valid = false;
        break;
      }
      score += pssm.scores[idx][j];
      scaledScore += pssm.scaledScores[idx][j];
    }

    if (!valid)
      continue;

    ++result.testedPositions;
    ++result.testedScoreCounts[static_cast<size_t>(scaledScore)];

    const double pValue =
        pssm.pValueByScaledScore[static_cast<size_t>(scaledScore)];
    if (score >= minRawScore && pValue <= maximumPValueForRetention) {
      MotifMatch m;
      m.position = (size_t)pos;
      m.matchLength = (size_t)motifLen;
      m.strand = strand;
      m.evidence.present = true;
      m.evidence.matrixId = pssm.id;
      m.evidence.matrixName = pssm.name;
      m.evidence.matrixSource = pssm.source;
      m.evidence.rawScore = score;
      m.evidence.scorePercent = scoreToPercent(score, pssm);
      m.evidence.background = pssm.background;
      m.evidence.motifPseudocount = pssm.motifPseudocount;
      m.evidence.statistics.pValue = pValue;
      m.evidence.statistics.scaledScore = scaledScore;
      m.evidence.statistics.scoreRange = pssm.scoreRange;
      m.evidence.statistics.scoreScale = pssm.scoreScale;
      m.evidence.statistics.scoreOffset = pssm.scoreOffset;

      size_t ctxStart =
          (m.position > contextSize) ? m.position - contextSize : 0;
      size_t ctxEnd =
          std::min(m.position + m.matchLength + contextSize, sequence.size());
      m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
      result.matches.push_back(m);
    }
  }

  return result;
}

PWMScanResult PWMScanner::scanWithStatistics(
    const std::string &sequence, const PSSM &pssm,
    double thresholdPercent, const std::string &chrId,
    bool searchPositive, bool searchNegative,
    double maximumPValueForRetention) {

  PWMScanResult result;
  result.testedScoreCounts.assign(pssm.pValueByScaledScore.size(), 0);

  double minRawScore = pssm.minScore;
  if (thresholdPercent >= 100.0) {
    minRawScore = pssm.maxScore;
  } else if (thresholdPercent > 0.0) {
    const double range = pssm.maxScore - pssm.minScore;
    minRawScore =
        pssm.minScore + (thresholdPercent / 100.0) * range;
  }

  if (searchPositive) {
    mergeScanResults(
        result, scanStrand(sequence, pssm, minRawScore, "+",
                           maximumPValueForRetention));
  }

  if (searchNegative) {
    std::string rcSeq = MotifFinder::reverseComplement(sequence);
    PWMScanResult negative =
        scanStrand(rcSeq, pssm, minRawScore, "-",
                   maximumPValueForRetention);

    for (auto &m : negative.matches) {
      m.position = sequence.size() - m.position - m.matchLength;

      size_t ctxStart = (m.position > 20) ? m.position - 20 : 0;
      size_t ctxEnd =
          std::min(m.position + m.matchLength + 20, sequence.size());
      m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
    }
    mergeScanResults(result, std::move(negative));
  }

  std::sort(result.matches.begin(), result.matches.end(),
            [](const MotifMatch &a, const MotifMatch &b) {
              return a.position < b.position;
            });

  for (auto &match : result.matches)
    match.chr = chrId;

  return result;
}

void PWMScanner::mergeScanResults(PWMScanResult &destination,
                                  PWMScanResult source) {
  if (destination.testedScoreCounts.empty()) {
    destination.testedScoreCounts.assign(source.testedScoreCounts.size(), 0);
  }
  if (destination.testedScoreCounts.size() !=
      source.testedScoreCounts.size()) {
    return;
  }
  for (size_t score = 0; score < source.testedScoreCounts.size(); ++score)
    destination.testedScoreCounts[score] += source.testedScoreCounts[score];
  destination.testedPositions += source.testedPositions;
  destination.matches.insert(
      destination.matches.end(),
      std::make_move_iterator(source.matches.begin()),
      std::make_move_iterator(source.matches.end()));
}

void PWMScanner::applyBenjaminiHochberg(PWMScanResult &result,
                                        const PSSM &pssm) {
  if (result.testedPositions == 0 ||
      result.testedScoreCounts.size() != pssm.pValueByScaledScore.size()) {
    return;
  }

  std::vector<double> adjustedByScore(result.testedScoreCounts.size(), 1.0);
  size_t rank = 0;
  for (size_t score = result.testedScoreCounts.size(); score-- > 0;) {
    const size_t tiedTests = result.testedScoreCounts[score];
    if (tiedTests == 0)
      continue;
    rank += tiedTests;
    adjustedByScore[score] = std::min(
        1.0, pssm.pValueByScaledScore[score] *
                 static_cast<double>(result.testedPositions) /
                 static_cast<double>(rank));
  }

  double runningMinimum = 1.0;
  for (size_t score = 0; score < adjustedByScore.size(); ++score) {
    if (result.testedScoreCounts[score] == 0)
      continue;
    runningMinimum = std::min(runningMinimum, adjustedByScore[score]);
    adjustedByScore[score] = runningMinimum;
  }

  for (auto &match : result.matches) {
    MotifStatisticalEvidence &statistics = match.evidence.statistics;
    statistics.qValue =
        adjustedByScore[static_cast<size_t>(statistics.scaledScore)];
    statistics.testedPositions = result.testedPositions;
  }
}

std::vector<MotifMatch>
PWMScanner::scan(const std::string &sequence, const PSSM &pssm,
                 double thresholdPercent, const std::string &chrId,
                 bool searchPositive, bool searchNegative) {
  PWMScanResult result = scanWithStatistics(
      sequence, pssm, thresholdPercent, chrId, searchPositive, searchNegative);
  applyBenjaminiHochberg(result, pssm);
  return std::move(result.matches);
}
