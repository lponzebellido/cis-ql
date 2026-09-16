#include "BEDReader.h"
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace {

std::vector<std::string> splitFields(const std::string &line) {
  std::vector<std::string> fields;
  std::istringstream input(line);
  std::string field;
  while (input >> field)
    fields.push_back(field);
  return fields;
}

bool parseSize(const std::string &text, size_t &value) {
  if (text.empty() || text[0] == '-')
    return false;
  char *end = nullptr;
  errno = 0;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno != 0 || end == text.c_str() || *end != '\0')
    return false;
  value = static_cast<size_t>(parsed);
  return static_cast<unsigned long long>(value) == parsed;
}

bool parseDouble(const std::string &text, double &value) {
  char *end = nullptr;
  errno = 0;
  value = std::strtod(text.c_str(), &end);
  return errno == 0 && end != text.c_str() && *end == '\0' &&
         std::isfinite(value);
}

} // namespace

std::vector<GenomicRegion> BEDReader::read(const std::string &filename,
                                           const std::string &format,
                                           const std::string &trackAlias,
                                           std::string *error) {
  std::vector<GenomicRegion> regions;
  std::ifstream file(filename);
  if (!file.is_open()) {
    if (error)
      *error = "Cannot open track file '" + filename + "'.";
    return regions;
  }

  std::string line;
  size_t lineNumber = 0;
  while (std::getline(file, line)) {
    ++lineNumber;
    if (line.empty() || line[0] == '#' || line.compare(0, 5, "track") == 0 ||
        line.compare(0, 7, "browser") == 0)
      continue;

    const std::vector<std::string> fields = splitFields(line);
    const size_t required = format == "NARROWPEAK" ? 10 : 3;
    if (fields.size() < required) {
      if (error)
        *error = "Invalid " + format + " record at line " +
                 std::to_string(lineNumber) + ": expected at least " +
                 std::to_string(required) + " columns.";
      return {};
    }

    size_t start = 0;
    size_t end = 0;
    if (!parseSize(fields[1], start) || !parseSize(fields[2], end) ||
        start >= end) {
      if (error)
        *error = "Invalid zero-based half-open interval at line " +
                 std::to_string(lineNumber) + ".";
      return {};
    }

    GenomicRegion region;
    region.chr = fields[0];
    region.start = start;
    region.end = end;
    region.name = fields.size() >= 4 && fields[3] != "."
                      ? fields[3]
                      : trackAlias + "_" + std::to_string(regions.size());
    region.strand = fields.size() >= 6 ? fields[5] : ".";
    if (region.strand != "+" && region.strand != "-")
      region.strand = ".";
    region.type = format == "NARROWPEAK" ? "narrow_peak" : "bed_region";
    region.trackEvidence.present = true;
    region.trackEvidence.trackAlias = trackAlias;
    region.trackEvidence.source = filename;
    region.trackEvidence.format = format;

    if (fields.size() >= 5) {
      if (!parseDouble(fields[4], region.trackEvidence.score) ||
          region.trackEvidence.score < 0.0 ||
          region.trackEvidence.score > 1000.0) {
        if (error)
          *error = "Invalid BED score at line " +
                   std::to_string(lineNumber) +
                   ": expected a finite value from 0 to 1000.";
        return {};
      }
      region.trackEvidence.hasScore = true;
    }

    if (format == "NARROWPEAK") {
      double value = 0.0;
      if (!parseDouble(fields[6], value)) {
        if (error)
          *error = "Invalid narrowPeak signalValue at line " +
                   std::to_string(lineNumber) + ".";
        return {};
      }
      if (value >= 0.0) {
        region.trackEvidence.hasSignalValue = true;
        region.trackEvidence.signalValue = value;
      }
      if (!parseDouble(fields[7], value)) {
        if (error)
          *error = "Invalid narrowPeak pValue at line " +
                   std::to_string(lineNumber) + ".";
        return {};
      }
      if (value >= 0.0) {
        region.trackEvidence.hasMinusLog10PValue = true;
        region.trackEvidence.minusLog10PValue = value;
      }
      if (!parseDouble(fields[8], value)) {
        if (error)
          *error = "Invalid narrowPeak qValue at line " +
                   std::to_string(lineNumber) + ".";
        return {};
      }
      if (value >= 0.0) {
        region.trackEvidence.hasMinusLog10QValue = true;
        region.trackEvidence.minusLog10QValue = value;
      }
      if (fields[9] != "-1") {
        size_t peak = 0;
        if (!parseSize(fields[9], peak) || peak >= region.length()) {
          if (error)
            *error = "Invalid narrowPeak summit offset at line " +
                     std::to_string(lineNumber) + ".";
          return {};
        }
        region.trackEvidence.hasPeak = true;
        region.trackEvidence.peakOffset = peak;
      }
    }
    regions.push_back(region);
  }
  return regions;
}
