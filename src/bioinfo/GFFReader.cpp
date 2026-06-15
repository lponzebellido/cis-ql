#include "GFFReader.h"
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

static std::string extractAttribute(const std::string& attrs, const std::string& key) {
  size_t pos = attrs.find(key + "=");
  if (pos == std::string::npos) return "";
  pos += key.length() + 1;
  size_t end = attrs.find(';', pos);
  if (end == std::string::npos) end = attrs.length();
  return attrs.substr(pos, end - pos);
}

std::vector<GenomicRegion> GFFReader::read(const std::string& filename) {
  std::vector<GenomicRegion> regions;
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Runtime Error: Cannot open GFF3 file '" << filename << "'." << std::endl;
    return regions;
  }

  std::string line;
  std::vector<GenomicRegion> genes;

  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (line.find("##FASTA") != std::string::npos) break;

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

    GenomicRegion region;
    region.chr = seqid;
    region.start = start - 1;
    region.end = end;
    region.strand = strandStr;
    region.type = type;
    region.name = extractAttribute(attributes, "Name");
    if (region.name.empty()) {
      region.name = extractAttribute(attributes, "ID");
    }

    regions.push_back(region);

    if (type == "gene") {
      genes.push_back(region);
    }
  }

  for (const auto& gene : genes) {
    GenomicRegion promoter;
    promoter.chr = gene.chr;
    promoter.strand = gene.strand;
    promoter.type = "promoter";
    promoter.name = gene.name + "_promoter";
    if (gene.strand == "+") {
      promoter.start = (gene.start > 250) ? gene.start - 250 : 0;
      promoter.end = gene.start;
    } else {
      promoter.start = gene.end;
      promoter.end = gene.end + 250;
    }
    if (promoter.start < promoter.end) {
      regions.push_back(promoter);
    }

    GenomicRegion tss;
    tss.chr = gene.chr;
    tss.strand = gene.strand;
    tss.type = "TSS";
    tss.name = gene.name + "_TSS";
    if (gene.strand == "+") {
      tss.start = gene.start;
      tss.end = gene.start + 1;
    } else {
      tss.start = gene.end - 1;
      tss.end = gene.end;
    }
    regions.push_back(tss);
  }

  return regions;
}

std::vector<GenomicRegion> GFFReader::filterByType(
    const std::vector<GenomicRegion>& regions, const std::string& type) {
  std::string gffType = mapEntityToGFFType(type);
  std::vector<GenomicRegion> result;
  for (const auto& r : regions) {
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
