#include "RegulatoryRegions.h"
#include <algorithm>

namespace {

size_t addAndClamp(size_t anchor, size_t distance, size_t chromosomeLength) {
  if (anchor >= chromosomeLength || distance >= chromosomeLength - anchor)
    return chromosomeLength;
  return anchor + distance;
}

std::string sourceLabel(const GenomicRegion &source) {
  if (!source.name.empty())
    return source.name;
  return source.chr + "_" + std::to_string(source.start) + "_" +
         std::to_string(source.end);
}

} // namespace

bool RegulatoryRegions::buildPromoters(
    const std::vector<GenomicRegion> &sources,
    const std::unordered_map<std::string, size_t> &chromosomeLengths,
    size_t upstream, size_t downstream,
    std::vector<GenomicRegion> &promoters, std::string &error) {
  promoters.clear();
  error.clear();

  if (upstream == 0 && downstream == 0) {
    error = "A promoter window cannot have both distances set to zero.";
    return false;
  }

  // Validate the complete input before producing output. A partially built
  // promoter set would be scientifically ambiguous.
  for (const auto &source : sources) {
    const auto chromosome = chromosomeLengths.find(source.chr);
    if (chromosome == chromosomeLengths.end()) {
      error = "Source sequence identifier '" + source.chr +
              "' is not present in the active FASTA dataset.";
      return false;
    }
    if (source.strand != "+" && source.strand != "-") {
      error = "Cannot define a TSS-relative promoter for '" +
              sourceLabel(source) + "': strand must be '+' or '-'.";
      return false;
    }
    if (source.start > source.end || source.end > chromosome->second) {
      error = "Source region '" + sourceLabel(source) +
              "' lies outside chromosome '" + source.chr + "'.";
      return false;
    }
  }

  promoters.reserve(sources.size());
  for (const auto &source : sources) {
    const size_t chromosomeLength = chromosomeLengths.at(source.chr);
    const size_t tss = source.strand == "+" ? source.start : source.end;

    GenomicRegion promoter;
    promoter.chr = source.chr;
    promoter.strand = source.strand;
    promoter.type = "promoter";
    promoter.name = sourceLabel(source) + "_promoter";

    if (source.strand == "+") {
      promoter.start = tss > upstream ? tss - upstream : 0;
      promoter.end = addAndClamp(tss, downstream, chromosomeLength);
    } else {
      promoter.start = tss > downstream ? tss - downstream : 0;
      promoter.end = addAndClamp(tss, upstream, chromosomeLength);
    }

    if (promoter.start < promoter.end)
      promoters.push_back(std::move(promoter));
  }

  std::sort(promoters.begin(), promoters.end());
  return true;
}
