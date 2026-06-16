#include "MotifFinder.h"
#include <algorithm>
#include <regex>

std::string MotifFinder::reverseComplement(const std::string& seq) {
  std::string rc(seq.size(), 'N');
  for (size_t i = 0; i < seq.size(); ++i) {
    char c = seq[seq.size() - 1 - i];
    switch (c) {
      case 'A': case 'a': rc[i] = 'T'; break;
      case 'T': case 't': rc[i] = 'A'; break;
      case 'C': case 'c': rc[i] = 'G'; break;
      case 'G': case 'g': rc[i] = 'C'; break;
      default: rc[i] = 'N'; break;
    }
  }
  return rc;
}

bool MotifFinder::isRegexPattern(const std::string& pattern) {
  for (char c : pattern) {
    if (c == '[' || c == ']' || c == '(' || c == ')' ||
        c == '|' || c == '*' || c == '+' || c == '?' ||
        c == '.' || c == '{' || c == '}' || c == '^' || c == '$') {
      return true;
    }
  }
  return false;
}

std::vector<int> MotifFinder::computeKMPTable(const std::string& pattern) {
  int m = (int)pattern.size();
  std::vector<int> table(m, 0);
  int len = 0;
  int i = 1;
  while (i < m) {
    if (std::toupper(pattern[i]) == std::toupper(pattern[len])) {
      len++;
      table[i] = len;
      i++;
    } else {
      if (len != 0) {
        len = table[len - 1];
      } else {
        table[i] = 0;
        i++;
      }
    }
  }
  return table;
}

std::vector<size_t> MotifFinder::kmpSearch(const std::string& text, const std::string& pattern) {
  std::vector<size_t> positions;
  int n = (int)text.size();
  int m = (int)pattern.size();
  if (m == 0 || n == 0 || m > n) return positions;

  std::vector<int> table = computeKMPTable(pattern);
  int i = 0, j = 0;

  while (i < n) {
    if (std::toupper(text[i]) == std::toupper(pattern[j])) {
      i++;
      j++;
    }
    if (j == m) {
      positions.push_back((size_t)(i - j));
      j = table[j - 1];
    } else if (i < n && std::toupper(text[i]) != std::toupper(pattern[j])) {
      if (j != 0) {
        j = table[j - 1];
      } else {
        i++;
      }
    }
  }
  return positions;
}

std::vector<MotifMatch> MotifFinder::regexSearch(
    const std::string& text, const std::string& pattern,
    const std::string& strand) {

  std::vector<MotifMatch> matches;
  size_t contextSize = 20;

  try {
    std::regex re(pattern, std::regex::icase | std::regex::optimize);
    auto it = std::sregex_iterator(text.begin(), text.end(), re);
    auto end = std::sregex_iterator();

    for (; it != end; ++it) {
      MotifMatch m;
      m.position = (size_t)it->position();
      m.matchLength = (size_t)it->length();
      m.strand = strand;
      size_t ctxStart = (m.position > contextSize) ? m.position - contextSize : 0;
      size_t ctxEnd = std::min(m.position + m.matchLength + contextSize, text.size());
      m.context = text.substr(ctxStart, ctxEnd - ctxStart);
      matches.push_back(m);
    }
  } catch (const std::regex_error& e) {
    return matches;
  }

  return matches;
}

std::vector<MotifMatch> MotifFinder::findAll(
    const std::string& sequence, const std::string& pattern,
    const std::string& chr, bool searchNegativeStrand) {

  std::vector<MotifMatch> matches;
  size_t contextSize = 20;

  if (isRegexPattern(pattern)) {
    matches = regexSearch(sequence, pattern, "+");

    if (searchNegativeStrand) {
      std::string rcSeq = reverseComplement(sequence);
      auto negMatches = regexSearch(rcSeq, pattern, "-");
      for (auto& m : negMatches) {
        m.position = sequence.size() - m.position - m.matchLength;
      }
      matches.insert(matches.end(), negMatches.begin(), negMatches.end());
    }
  } else {
    std::vector<size_t> posPositive = kmpSearch(sequence, pattern);
    for (size_t pos : posPositive) {
      MotifMatch m;
      m.position = pos;
      m.matchLength = pattern.size();
      m.strand = "+";
      size_t ctxStart = (pos > contextSize) ? pos - contextSize : 0;
      size_t ctxEnd = std::min(pos + pattern.size() + contextSize, sequence.size());
      m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
      matches.push_back(m);
    }

    if (searchNegativeStrand) {
      std::string rc = reverseComplement(pattern);
      if (rc != pattern) {
        std::vector<size_t> posNegative = kmpSearch(sequence, rc);
        for (size_t pos : posNegative) {
          MotifMatch m;
          m.position = pos;
          m.matchLength = rc.size();
          m.strand = "-";
          size_t ctxStart = (pos > contextSize) ? pos - contextSize : 0;
          size_t ctxEnd = std::min(pos + rc.size() + contextSize, sequence.size());
          m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
          matches.push_back(m);
        }
      }
    }
  }

  std::sort(matches.begin(), matches.end(),
            [](const MotifMatch& a, const MotifMatch& b) {
              return a.position < b.position;
            });

  return matches;
}

std::vector<MotifMatch> MotifFinder::findInWindow(
    const std::string& sequence, const std::string& pattern,
    size_t windowStart, size_t windowEnd, const std::string& chr) {

  if (windowStart >= sequence.size()) return {};
  if (windowEnd > sequence.size()) windowEnd = sequence.size();
  if (windowStart >= windowEnd) return {};

  std::string window = sequence.substr(windowStart, windowEnd - windowStart);
  std::vector<MotifMatch> matches = findAll(window, pattern, chr, true);

  for (auto& m : matches) {
    m.position += windowStart;
    size_t ctxStart = (m.position > 20) ? m.position - 20 : 0;
    size_t ctxEnd = std::min(m.position + m.matchLength + 20, sequence.size());
    m.context = sequence.substr(ctxStart, ctxEnd - ctxStart);
  }

  return matches;
}
