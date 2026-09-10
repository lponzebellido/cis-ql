#include "PWMScanner.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

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

PSSM PWMScanner::computePSSM(const PWMatrix &pwm, double bgA, double bgC,
                             double bgG, double bgT) {
  PSSM pssm;
  pssm.name = pwm.name;
  pssm.id = pwm.id;
  pssm.source = pwm.source;
  pssm.length = pwm.length;
  if (pwm.length <= 0 || pwm.counts.size() != 4 || bgA <= 0.0 || bgC <= 0.0 ||
      bgG <= 0.0 || bgT <= 0.0) {
    pssm.length = 0;
    return pssm;
  }
  pssm.scores.resize(4);
  pssm.maxScore = 0.0;
  pssm.minScore = 0.0;

  double bg[4] = {bgA, bgC, bgG, bgT};
  const double pseudocount = 0.1;

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

    double colMax = -1e9;
    double colMin = 1e9;
    for (int i = 0; i < 4; i++) {
      double freq =
          (pwm.counts[i][j] + pseudocount) / (colTotal + 4.0 * pseudocount);
      pssm.scores[i][j] = std::log2(freq / bg[i]);

      if (pssm.scores[i][j] > colMax)
        colMax = pssm.scores[i][j];
      if (pssm.scores[i][j] < colMin)
        colMin = pssm.scores[i][j];
    }
    pssm.maxScore += colMax;
    pssm.minScore += colMin;
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

std::vector<MotifMatch> PWMScanner::scanStrand(const std::string &sequence,
                                               const PSSM &pssm,
                                               double minRawScore,
                                               const std::string &strand) {

  std::vector<MotifMatch> matches;
  int seqLen = (int)sequence.size();
  int motifLen = pssm.length;

  if (seqLen < motifLen)
    return matches;

  size_t contextSize = 20;

  for (int pos = 0; pos <= seqLen - motifLen; pos++) {
    double score = 0.0;
    bool valid = true;

    for (int j = 0; j < motifLen; j++) {
      int idx = nucToIndex(sequence[pos + j]);
      if (idx < 0) {
        valid = false;
        break;
      }
      score += pssm.scores[idx][j];
    }

    if (valid && score >= minRawScore) {
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

      size_t ctxStart =
          (m.position > contextSize) ? m.position - contextSize : 0;
      size_t ctxEnd =
          std::min(m.position + m.matchLength + contextSize, sequence.size());
      m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
      matches.push_back(m);
    }
  }

  return matches;
}

std::vector<MotifMatch>
PWMScanner::scan(const std::string &sequence, const PSSM &pssm,
                 double thresholdPercent, const std::string &chrId,
                 bool searchPositive, bool searchNegative) {

  std::vector<MotifMatch> allMatches;

  double range = pssm.maxScore - pssm.minScore;
  double minRawScore = pssm.minScore + (thresholdPercent / 100.0) * range;

  if (searchPositive) {
    auto posMatches = scanStrand(sequence, pssm, minRawScore, "+");
    allMatches.insert(allMatches.end(), posMatches.begin(), posMatches.end());
  }

  if (searchNegative) {
    std::string rcSeq = MotifFinder::reverseComplement(sequence);
    auto negMatches = scanStrand(rcSeq, pssm, minRawScore, "-");

    for (auto &m : negMatches) {
      m.position = sequence.size() - m.position - m.matchLength;

      size_t ctxStart = (m.position > 20) ? m.position - 20 : 0;
      size_t ctxEnd =
          std::min(m.position + m.matchLength + 20, sequence.size());
      m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
    }
    allMatches.insert(allMatches.end(), negMatches.begin(), negMatches.end());
  }

  std::sort(allMatches.begin(), allMatches.end(),
            [](const MotifMatch &a, const MotifMatch &b) {
              return a.position < b.position;
            });

  for (auto &match : allMatches)
    match.chr = chrId;

  return allMatches;
}
