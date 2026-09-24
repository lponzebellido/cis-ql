#include "GFFReader.h"
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

std::string GFFReader::mapEntityToGFFType(const std::string& entity) {
  if (entity == "GENE") return "gene";
  if (entity == "EXON") return "exon";
  if (entity == "INTRON") return "intron";
  if (entity == "PROMOTER") return "promoter";
  if (entity == "ENHANCER") return "enhancer";
  if (entity == "UTR") return "UTR";
  if (entity == "TSS") return "TSS";
  if (entity == "REGION") return "region";
  if (entity == "CDS") return "CDS";
  return entity;
}

static int hexValue(char value) {
  if (value >= '0' && value <= '9') return value - '0';
  if (value >= 'A' && value <= 'F') return value - 'A' + 10;
  if (value >= 'a' && value <= 'f') return value - 'a' + 10;
  return -1;
}

static std::string decodeAttribute(const std::string& value) {
  std::string decoded;
  decoded.reserve(value.size());
  for (size_t index = 0; index < value.size(); ++index) {
    if (value[index] == '%' && index + 2 < value.size()) {
      const int high = hexValue(value[index + 1]);
      const int low = hexValue(value[index + 2]);
      if (high >= 0 && low >= 0) {
        decoded += static_cast<char>((high << 4) | low);
        index += 2;
        continue;
      }
    }
    decoded += value[index];
  }
  return decoded;
}

static std::map<std::string, std::string>
parseAttributes(const std::string& attributes) {
  std::map<std::string, std::string> parsed;
  size_t start = 0;
  while (start <= attributes.size()) {
    const size_t end = attributes.find(';', start);
    const std::string field = attributes.substr(
        start, end == std::string::npos ? std::string::npos : end - start);
    const size_t separator = field.find('=');
    if (separator != std::string::npos && separator > 0) {
      const std::string key = decodeAttribute(field.substr(0, separator));
      const std::string value = decodeAttribute(field.substr(separator + 1));
      parsed[key] = value;
    }
    if (end == std::string::npos) break;
    start = end + 1;
  }
  return parsed;
}

static std::vector<std::string> parseParents(const std::string& attributes) {
  std::vector<std::string> parents;
  size_t fieldStart = 0;
  while (fieldStart <= attributes.size()) {
    const size_t fieldEnd = attributes.find(';', fieldStart);
    const std::string field = attributes.substr(
        fieldStart, fieldEnd == std::string::npos
                        ? std::string::npos
                        : fieldEnd - fieldStart);
    if (field.compare(0, 7, "Parent=") == 0) {
      const std::string rawParents = field.substr(7);
      size_t parentStart = 0;
      while (parentStart <= rawParents.size()) {
        const size_t parentEnd = rawParents.find(',', parentStart);
        parents.push_back(decodeAttribute(rawParents.substr(
            parentStart, parentEnd == std::string::npos
                             ? std::string::npos
                             : parentEnd - parentStart)));
        if (parentEnd == std::string::npos) break;
        parentStart = parentEnd + 1;
      }
      break;
    }
    if (fieldEnd == std::string::npos) break;
    fieldStart = fieldEnd + 1;
  }
  return parents;
}

std::vector<GenomicRegion> GFFReader::read(const std::string& filename) {
  std::vector<GenomicRegion> regions;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Runtime Error: Cannot open GFF3 file '" << filename << "'." << std::endl;
    return regions;
  }

  std::string line;
  while (std::getline(file, line)) {
    if (line.find("##FASTA") != std::string::npos) break;
    if (line.empty() || line[0] == '#') continue;

    std::istringstream iss(line);
    std::string seqid, source, type;
    size_t start, end;
    std::string score, strandStr, phase, attributes;

    if (!(iss >> seqid >> source >> type >> start >> end >> score >> strandStr >> phase)) {
      continue;
    }
    std::getline(iss, attributes);
    if (!attributes.empty() && attributes[0] == '\t') {
      attributes = attributes.substr(1);
    }

    if (start == 0 || end < start) continue;
    if (!attributes.empty() && attributes.back() == '\r') {
      attributes.pop_back();
    }

    GenomicRegion region;
    region.chr = seqid;
    region.start = start - 1;
    region.end = end;
    region.strand = strandStr;
    region.type = type;
    region.annotationEvidence.present = true;
    region.annotationEvidence.source = source;
    region.annotationEvidence.score = score;
    region.annotationEvidence.phase = phase;
    region.annotationEvidence.attributes = parseAttributes(attributes);
    region.annotationEvidence.parents = parseParents(attributes);
    region.annotationEvidence.rawAttributes = attributes;
    const auto id = region.annotationEvidence.attributes.find("ID");
    if (id != region.annotationEvidence.attributes.end()) {
      region.annotationEvidence.id = id->second;
    }
    const auto name = region.annotationEvidence.attributes.find("Name");
    if (name != region.annotationEvidence.attributes.end()) {
      region.annotationEvidence.name = name->second;
    }
    region.name = region.annotationEvidence.name;
    if (region.name.empty()) {
      region.name = region.annotationEvidence.id;
    }

    regions.push_back(region);

  }

  return regions;
}

std::vector<GenomicRegion> GFFReader::filterByType(
    const std::vector<GenomicRegion>& regions, const std::string& type) {
  std::string gffType = mapEntityToGFFType(type);
  std::vector<GenomicRegion> result;
  for (const auto& r : regions) {
    if (type == "FEATURE") {
      result.push_back(r);
      continue;
    }
    std::string rTypeLower = r.type;
    std::string gffTypeLower = gffType;
    std::transform(rTypeLower.begin(), rTypeLower.end(), rTypeLower.begin(), ::tolower);
    std::transform(gffTypeLower.begin(), gffTypeLower.end(), gffTypeLower.begin(), ::tolower);
    if (rTypeLower == gffTypeLower) {
      result.push_back(r);
    }
  }
  return result;
}
