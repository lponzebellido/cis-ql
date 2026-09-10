#include "Interpreter.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <thread>
#include <fstream>

namespace {

bool conditionUsesOnly(const std::shared_ptr<IRCondition> &condition,
                       const std::set<std::string> &properties) {
  if (!condition)
    return true;
  if (condition->kind == IRCondition::Kind::SIMPLE)
    return properties.count(condition->property) > 0;
  return conditionUsesOnly(condition->left, properties) &&
         conditionUsesOnly(condition->right, properties);
}

void collectSimilarityReferences(
    const std::shared_ptr<IRCondition> &condition,
    std::set<std::string> &references) {
  if (!condition)
    return;
  if (condition->kind == IRCondition::Kind::SIMPLE) {
    if (condition->property == "SIMILARITY" &&
        !condition->reference.empty()) {
      references.insert(condition->reference);
    }
    return;
  }
  collectSimilarityReferences(condition->left, references);
  collectSimilarityReferences(condition->right, references);
}

std::string cleanTabularField(const std::string &value) {
  std::string cleaned = value;
  for (char &c : cleaned) {
    if (c == '\t' || c == '\n' || c == '\r')
      c = ' ';
  }
  return cleaned;
}

std::string gffAttributeEscape(const std::string &value) {
  static const char HEX[] = "0123456789ABCDEF";
  std::string escaped;
  for (const unsigned char c : value) {
    if (std::isalnum(c) || c == '_' || c == '-' || c == '.' || c == ':') {
      escaped += static_cast<char>(c);
    } else {
      escaped += '%';
      escaped += HEX[(c >> 4) & 0x0f];
      escaped += HEX[c & 0x0f];
    }
  }
  return escaped;
}

bool isSafeRelativeExportPath(const std::string &path) {
  if (path.empty() || path.front() == '/' || path.front() == '\\' ||
      (path.size() > 1 &&
       std::isalpha(static_cast<unsigned char>(path[0])) &&
       path[1] == ':')) {
    return false;
  }
  std::string component;
  for (size_t index = 0; index <= path.size(); ++index) {
    const bool separator =
        index == path.size() || path[index] == '/' || path[index] == '\\';
    if (!separator) {
      component += path[index];
      continue;
    }
    if (component == "..")
      return false;
    component.clear();
  }
  return true;
}

} // namespace

std::string Interpreter::stripQuotes(const std::string &s) const {
  if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

std::string Interpreter::jsonEscape(const std::string &s) const {
  std::string escaped;
  escaped.reserve(s.size());
  for (const unsigned char c : s) {
    switch (c) {
    case '"': escaped += "\\\""; break;
    case '\\': escaped += "\\\\"; break;
    case '\b': escaped += "\\b"; break;
    case '\f': escaped += "\\f"; break;
    case '\n': escaped += "\\n"; break;
    case '\r': escaped += "\\r"; break;
    case '\t': escaped += "\\t"; break;
    default:
      if (c >= 0x20)
        escaped += static_cast<char>(c);
      break;
    }
  }
  return escaped;
}

void Interpreter::reportRuntimeError(const std::string &message) {
  runtimeError = true;
  std::cerr << "Runtime Error: " << message << std::endl;
}

size_t Interpreter::toBasePairs(double value, const std::string &unit) {
  if (unit == "KB")
    return static_cast<size_t>(value * 1000.0);
  if (unit == "MB")
    return static_cast<size_t>(value * 1000000.0);
  return static_cast<size_t>(value);
}

void Interpreter::printMotifMatches(const std::vector<MotifMatch> &matches,
                                    int maxShow) {
  if (matches.empty()) {
    std::cout << "  No matches found." << std::endl;
    return;
  }
  std::cout << "  Found " << matches.size() << " match(es):" << std::endl;
  int shown = 0;
  for (const auto &m : matches) {
    if (shown >= maxShow) {
      std::cout << "  ... and " << (matches.size() - maxShow) << " more."
                << std::endl;
      break;
    }
    std::cout << "  [" << (shown + 1) << "] pos:" << m.position
              << "  strand:" << m.strand;
    if (m.evidence.present) {
      std::cout << "  matrix:" << m.evidence.matrixId
                << "  score:" << m.evidence.rawScore
                << "  relative:" << m.evidence.scorePercent << "%";
    }
    if (m.evidence.hasSourceRegion)
      std::cout << "  source:" << m.evidence.sourceRegionName
                << "  source-offset:" << m.evidence.relativeStart;
    std::cout << "  context: ..." << m.context << "..." << std::endl;
    shown++;
  }
}

void Interpreter::printRegions(const std::vector<GenomicRegion> &regions,
                               int maxShow) {
  if (regions.empty()) {
    std::cout << "  No regions found." << std::endl;
    return;
  }
  std::cout << "  Found " << regions.size() << " region(s):" << std::endl;
  int shown = 0;
  for (const auto &r : regions) {
    if (shown >= maxShow) {
      std::cout << "  ... and " << (regions.size() - maxShow) << " more."
                << std::endl;
      break;
    }
    std::cout << "  [" << (shown + 1) << "] " << r.chr << ":" << r.start << ".."
              << r.end << "  length:" << r.length() << " BP"
              << "  strand:" << r.strand << "  type:" << r.type;
    if (!r.name.empty())
      std::cout << "  name:" << r.name;
    std::cout << std::endl;
    if (!r.sequence.empty()) {
      std::string display = r.sequence;
      if (display.size() > 60) {
        display = display.substr(0, 60) + "...";
      }
      std::cout << "      seq: " << display << std::endl;
    }
    shown++;
  }
}

std::vector<GenomicRegion>
Interpreter::resolveEntity(const std::string &entity) {
  if (namedRegions.count(entity)) {
    return namedRegions[entity];
  }
  const auto annotations = annotationDatasets.find(activeAnnotationAlias);
  if (annotations == annotationDatasets.end())
    return {};
  return GFFReader::filterByType(annotations->second, entity);
}

GenomicRegion
Interpreter::motifMatchToRegion(const MotifMatch &match,
                                const std::string &alias) const {
  GenomicRegion region;
  region.chr = match.chr;
  if (region.chr.empty()) {
    const auto active = sequenceDatasets.find(activeSequenceAlias);
    if (active != sequenceDatasets.end() && !active->second.empty())
      region.chr = active->second.front().sequenceId;
  }
  region.start = match.position;
  region.end = match.position + match.matchLength;
  region.strand = match.strand;
  region.type = match.evidence.present ? "motif_hit" : alias;
  region.name = alias + "_" + region.chr + "_" +
                std::to_string(match.position);
  if (match.evidence.present)
    region.name += match.strand == "-" ? "_minus" : "_plus";
  if (match.evidence.hasSourceRegion &&
      !match.evidence.sourceRegionName.empty()) {
    region.name += "_in_" + match.evidence.sourceRegionName;
  }
  region.motifEvidence = match.evidence;

  const auto chromosomeMap = sequenceChrMaps.find(activeSequenceAlias);
  if (chromosomeMap != sequenceChrMaps.end()) {
    const auto sequence = chromosomeMap->second.find(region.chr);
    if (sequence != chromosomeMap->second.end() &&
        region.end <= sequence->second.sequence.size()) {
      region.sequence = sequence->second.sequence.substr(
          region.start, region.end - region.start);
    }
  }
  return region;
}

void Interpreter::executeLoadSeq(const IRInstruction &instr) {
  std::string filename = stripQuotes(instr.arg1);
  std::string alias = instr.arg2;
  if (debugMode) {
    std::cout << "> LOAD SEQUENCE \"" << filename << "\" AS " << alias
              << std::endl;
  }

  auto records = FastaReader::read(filename);
  if (records.empty()) {
    reportRuntimeError("No sequences found in " + filename + ".");
    return;
  }

  sequenceDatasets[alias] = records;
  for (const auto &rec : records) {
    sequenceChrMaps[alias][rec.sequenceId] = rec;
  }
  activeSequenceAlias = alias;
  if (debugMode) {
    size_t totalBP = 0;
    for (const auto &rec : records) totalBP += rec.sequence.size();
    std::cout << "  Loaded " << records.size() << " sequence(s) / chromosome(s), total "
              << totalBP << " base pairs from " << filename << std::endl;
  }
}

void Interpreter::executeLoadAnnot(const IRInstruction &instr) {
  std::string filename = stripQuotes(instr.arg1);
  std::string alias = instr.arg2;
  if (debugMode) {
    std::cout << "> LOAD ANNOTATION \"" << filename << "\" AS " << alias
              << std::endl;
  }

  auto annotations = GFFReader::read(filename);
  if (annotations.empty()) {
    reportRuntimeError("No annotation features found in " + filename + ".");
    return;
  }
  annotationDatasets[alias] = annotations;
  activeAnnotationAlias = alias;
  if (debugMode) {
    std::cout << "  Loaded " << annotations.size()
              << " annotation features from " << filename << std::endl;

    std::unordered_map<std::string, int> typeCounts;
    for (const auto &a : annotations) {
      typeCounts[a.type]++;
    }
    for (const auto &pair : typeCounts) {
      std::cout << "    " << pair.first << ": " << pair.second << std::endl;
    }
  }
}

void Interpreter::executeUseDataset(const IRInstruction &instr) {
  const std::string &alias = instr.arg1;
  if (instr.opcode == IROpCode::USE_SEQUENCE) {
    if (!sequenceDatasets.count(alias)) {
      reportRuntimeError("Sequence dataset '" + alias + "' is not loaded.");
      return;
    }
    activeSequenceAlias = alias;
    if (debugMode)
      std::cout << "> USE SEQUENCE " << alias << std::endl;
    return;
  }

  if (!annotationDatasets.count(alias)) {
    reportRuntimeError("Annotation dataset '" + alias + "' is not loaded.");
    return;
  }
  activeAnnotationAlias = alias;
  if (debugMode)
    std::cout << "> USE ANNOTATION " << alias << std::endl;
}

void Interpreter::executeExport(const IRInstruction &instr) {
  const std::string &alias = instr.arg1;
  const std::string filename = stripQuotes(instr.arg2);
  const std::string &format = instr.arg3;

  if (!isSafeRelativeExportPath(filename)) {
    reportRuntimeError(
        "Export paths must be relative to the query workspace and cannot "
        "contain '..'.");
    return;
  }

  const auto regionsIt = resultSets.find(alias);
  const auto gcIt = gcResults.find(alias);
  if (regionsIt == resultSets.end() && gcIt == gcResults.end()) {
    reportRuntimeError("Result alias '" + alias + "' is not available.");
    return;
  }
  if (gcIt != gcResults.end() && format != "TSV") {
    reportRuntimeError("GC profiles can currently be exported only as TSV.");
    return;
  }
  if (regionsIt != resultSets.end() && format != "BED" &&
      format != "GFF3" && format != "TSV") {
    reportRuntimeError("Unsupported export format '" + format + "'.");
    return;
  }

  std::ofstream out(filename);
  if (!out.is_open()) {
    reportRuntimeError("Could not open export file '" + filename + "'.");
    return;
  }

  if (gcIt != gcResults.end()) {
    out << "chromosome\tstart\tgc_percent\n";
    for (const auto &window : gcIt->second) {
      out << cleanTabularField(window.chr) << '\t' << window.position << '\t'
          << window.gcPercent << '\n';
    }
  } else if (format == "BED") {
    for (const auto &region : regionsIt->second) {
      const std::string name =
          region.name.empty() ? "." : cleanTabularField(region.name);
      const std::string strand =
          (region.strand == "+" || region.strand == "-") ? region.strand : ".";
      int bedScore = 0;
      if (region.motifEvidence.present) {
        bedScore = static_cast<int>(
            std::round(region.motifEvidence.scorePercent * 10.0));
        bedScore = std::max(0, std::min(1000, bedScore));
      }
      out << cleanTabularField(region.chr) << '\t' << region.start << '\t'
          << region.end << '\t' << name << '\t' << bedScore << '\t'
          << strand << '\n';
    }
  } else if (format == "GFF3") {
    out << "##gff-version 3\n";
    size_t generatedId = 0;
    for (const auto &region : regionsIt->second) {
      const std::string type = region.type.empty() ? "region" : region.type;
      const std::string strand =
          (region.strand == "+" || region.strand == "-") ? region.strand : ".";
      const std::string name =
          region.name.empty() ? alias + "_" + std::to_string(generatedId++)
                              : region.name;
      const std::string score = region.motifEvidence.present
                                    ? std::to_string(
                                          region.motifEvidence.rawScore)
                                    : ".";
      out << cleanTabularField(region.chr) << "\tCis-QL\t"
          << cleanTabularField(type) << '\t' << (region.start + 1) << '\t'
          << region.end << '\t' << score << '\t' << strand << "\t.\tID="
          << gffAttributeEscape(name) << ";Name="
          << gffAttributeEscape(name);
      if (region.motifEvidence.present) {
        out << ";MatrixAlias="
            << gffAttributeEscape(region.motifEvidence.matrixAlias)
            << ";MatrixID="
            << gffAttributeEscape(region.motifEvidence.matrixId)
            << ";MatrixName="
            << gffAttributeEscape(region.motifEvidence.matrixName)
            << ";MatrixSource="
            << gffAttributeEscape(region.motifEvidence.matrixSource)
            << ";ScorePercent=" << region.motifEvidence.scorePercent;
        if (region.motifEvidence.hasSourceRegion) {
          out << ";SourceRegion="
              << gffAttributeEscape(
                     region.motifEvidence.sourceRegionName)
              << ";SourceRegionType="
              << gffAttributeEscape(
                     region.motifEvidence.sourceRegionType)
              << ";RelativeStart="
              << region.motifEvidence.relativeStart;
        }
      }
      out << '\n';
    }
  } else if (format == "TSV") {
    out << "chromosome\tstart\tend\tstrand\ttype\tname\tlength"
           "\tmatrix_alias\tmatrix_id\tmatrix_name\tmatrix_source"
           "\traw_score\tscore_percent\tsource_region\tsource_type"
           "\tsource_start\tsource_end\trelative_start\n";
    for (const auto &region : regionsIt->second) {
      out << cleanTabularField(region.chr) << '\t' << region.start << '\t'
          << region.end << '\t' << cleanTabularField(region.strand) << '\t'
          << cleanTabularField(region.type) << '\t'
          << cleanTabularField(region.name) << '\t' << region.length()
          << '\t';
      if (region.motifEvidence.present) {
        out << cleanTabularField(region.motifEvidence.matrixAlias) << '\t'
            << cleanTabularField(region.motifEvidence.matrixId) << '\t'
            << cleanTabularField(region.motifEvidence.matrixName) << '\t'
            << cleanTabularField(region.motifEvidence.matrixSource) << '\t'
            << region.motifEvidence.rawScore << '\t'
            << region.motifEvidence.scorePercent << '\t';
        if (region.motifEvidence.hasSourceRegion) {
          out << cleanTabularField(
                     region.motifEvidence.sourceRegionName)
              << '\t' << cleanTabularField(
                             region.motifEvidence.sourceRegionType)
              << '\t' << region.motifEvidence.sourceRegionStart << '\t'
              << region.motifEvidence.sourceRegionEnd << '\t'
              << region.motifEvidence.relativeStart;
        } else {
          out << "\t\t\t\t";
        }
      } else {
        for (int emptyColumn = 0; emptyColumn < 10; ++emptyColumn)
          out << '\t';
      }
      out << '\n';
    }
  }

  if (!out.good()) {
    reportRuntimeError("Failed while writing export file '" + filename + "'.");
    return;
  }
  if (debugMode) {
    std::cout << "> EXPORT " << alias << " TO \"" << filename
              << "\" FORMAT " << format << std::endl;
  }
}

void Interpreter::executeDefinePromoters(const IRInstruction &instr) {
  const std::string &source = instr.arg1;
  const std::string &alias = instr.arg4;

  auto parseDistance = [this](const std::string &value) {
    const size_t separator = value.find(' ');
    const double number =
        std::atof(value.substr(0, separator).c_str());
    const std::string unit = separator == std::string::npos
                                 ? "BP"
                                 : value.substr(separator + 1);
    return toBasePairs(number, unit);
  };

  const size_t upstream = parseDistance(instr.arg2);
  const size_t downstream = parseDistance(instr.arg3);
  const auto activeGenome = sequenceChrMaps.find(activeSequenceAlias);
  if (activeGenome == sequenceChrMaps.end()) {
    reportRuntimeError("DEFINE PROMOTERS requires an active FASTA dataset.");
    return;
  }

  std::unordered_map<std::string, size_t> chromosomeLengths;
  for (const auto &entry : activeGenome->second)
    chromosomeLengths[entry.first] = entry.second.sequence.size();

  const std::vector<GenomicRegion> sources = resolveEntity(source);
  std::vector<GenomicRegion> promoters;
  std::string error;
  if (!RegulatoryRegions::buildPromoters(sources, chromosomeLengths,
                                         upstream, downstream, promoters,
                                         error)) {
    reportRuntimeError(error);
    return;
  }

  for (auto &promoter : promoters) {
    const auto sequence = activeGenome->second.find(promoter.chr);
    if (sequence != activeGenome->second.end()) {
      promoter.sequence = sequence->second.sequence.substr(
          promoter.start, promoter.end - promoter.start);
    }
  }

  namedRegions[alias] = promoters;
  resultSets[alias] = std::move(promoters);
  if (debugMode) {
    std::cout << "> DEFINE PROMOTERS OF " << source
              << " FROM TSS UPSTREAM " << instr.arg2
              << " DOWNSTREAM " << instr.arg3 << " AS " << alias
              << std::endl;
    std::cout << "  Defined " << resultSets[alias].size()
              << " explicit promoter interval(s)." << std::endl;
  }
}

void Interpreter::executeFindMotif(const IRInstruction &instr) {
  currentFind = FindContext();
  currentFind.pattern = stripQuotes(instr.arg1);
  currentFind.hasWithin = false;
  if (!MotifFinder::isValidPattern(currentFind.pattern)) {
    reportRuntimeError("Invalid or empty motif pattern '" +
                       currentFind.pattern + "'.");
  }
}

void Interpreter::executeFindOptWithin(const IRInstruction &instr) {
  currentFind.hasWithin = true;
  currentFind.withinDistance = std::atof(instr.arg1.c_str());
  currentFind.withinUnit = instr.arg2;
  currentFind.withinDirection = instr.arg3;
  currentFind.withinEntity = instr.arg4;
  currentFind.withinTarget = stripQuotes(instr.arg5);
}

void Interpreter::executeFindOptStrand(const IRInstruction &instr) {
  currentFind.strandFilter = instr.arg1;
}

void Interpreter::executeFindOptChr(const IRInstruction &instr) {
  currentFind.chrFilter = stripQuotes(instr.arg1);
}

void Interpreter::executeFindExec(const IRInstruction &instr) {
  std::string resultId = instr.arg1;
  if (debugMode) {
    std::cout << "> FIND MOTIF \"" << currentFind.pattern << "\"";
    if (currentFind.hasWithin) {
      std::cout << " WITHIN " << currentFind.withinDistance << " "
                << currentFind.withinUnit << " " << currentFind.withinDirection
                << " FROM " << currentFind.withinEntity;
      if (!currentFind.withinTarget.empty()) {
        std::cout << " \"" << currentFind.withinTarget << "\"";
      }
    }
    if (!currentFind.strandFilter.empty()) {
      std::cout << " STRAND " << currentFind.strandFilter;
    }
    if (!currentFind.chrFilter.empty()) {
      std::cout << " CHR \"" << currentFind.chrFilter << "\"";
    }
    std::cout << std::endl;
  }

  if (sequenceDatasets.empty()) {
    reportRuntimeError("No sequence loaded.");
    return;
  }

  const auto &targetDataset = sequenceDatasets.count(activeSequenceAlias)
                                  ? sequenceDatasets[activeSequenceAlias]
                                  : sequenceDatasets.begin()->second;

  bool searchNeg = (currentFind.strandFilter != "POSITIVE");

  std::vector<MotifMatch> matches;

  if (currentFind.hasWithin) {
    size_t distBP =
        toBasePairs(currentFind.withinDistance, currentFind.withinUnit);

    std::vector<GenomicRegion> referenceRegions =
        resolveEntity(currentFind.withinEntity);

    std::vector<GenomicRegion> targets;
    if (!currentFind.withinTarget.empty()) {
      for (const auto &r : referenceRegions) {
        if (r.name == currentFind.withinTarget) {
          targets.push_back(r);
        }
      }
    } else {
      targets = referenceRegions;
    }

    if (targets.empty() && !currentFind.withinTarget.empty()) {
      std::cout << "  Warning: No " << currentFind.withinEntity << " named \""
                << currentFind.withinTarget << "\" found." << std::endl;
      targets = referenceRegions;
    }

    for (const auto &target : targets) {
      if (!currentFind.chrFilter.empty() && target.chr != currentFind.chrFilter)
        continue;

      std::string seqData = "";
      if (sequenceChrMaps[activeSequenceAlias].count(target.chr)) {
        seqData = sequenceChrMaps[activeSequenceAlias][target.chr].sequence;
      } else {
        reportRuntimeError("Annotation sequence identifier '" + target.chr +
                           "' is not present in the active FASTA dataset.");
        continue;
      }
      if (seqData.empty()) continue;

      size_t windowStart, windowEnd;
      std::string effectiveStrand = target.strand.empty() ? "+" : target.strand;
      if (currentFind.withinDirection == "UPSTREAM") {
        if (effectiveStrand == "+") {
          windowStart = (target.start > distBP) ? target.start - distBP : 0;
          windowEnd = target.start;
        } else {
          windowStart = target.end;
          windowEnd = std::min(target.end + distBP, seqData.size());
        }
      } else {
        if (effectiveStrand == "+") {
          windowStart = target.end;
          windowEnd = std::min(target.end + distBP, seqData.size());
        } else {
          windowStart = (target.start > distBP) ? target.start - distBP : 0;
          windowEnd = target.start;
        }
      }
      auto windowMatches =
          MotifFinder::findInWindow(seqData, currentFind.pattern,
                                    windowStart, windowEnd, target.chr);
      matches.insert(matches.end(), windowMatches.begin(), windowMatches.end());
    }
  } else {
    std::vector<std::future<std::vector<MotifMatch>>> futures;
    for (const auto &seqRec : targetDataset) {
      if (!currentFind.chrFilter.empty() && seqRec.sequenceId != currentFind.chrFilter)
        continue;

      const FastaRecord *record = &seqRec;
      std::string chrId = seqRec.sequenceId;
      std::string pat = currentFind.pattern;

      futures.push_back(std::async(std::launch::async, [record, pat, chrId, searchNeg]() {
        return MotifFinder::findAll(record->sequence, pat, chrId, searchNeg);
      }));
    }

    for (auto &f : futures) {
      auto res = f.get();
      matches.insert(matches.end(), res.begin(), res.end());
    }
  }

  if (currentFind.strandFilter == "NEGATIVE") {
    std::vector<MotifMatch> filtered;
    for (const auto &m : matches) {
      if (m.strand == "-") filtered.push_back(m);
    }
    matches = filtered;
  } else if (currentFind.strandFilter == "POSITIVE") {
    std::vector<MotifMatch> filtered;
    for (const auto &m : matches) {
      if (m.strand == "+") filtered.push_back(m);
    }
    matches = filtered;
  }

  std::sort(matches.begin(), matches.end(),
            [](const MotifMatch &left, const MotifMatch &right) {
              if (left.chr != right.chr)
                return left.chr < right.chr;
              if (left.position != right.position)
                return left.position < right.position;
              if (left.matchLength != right.matchLength)
                return left.matchLength < right.matchLength;
              return left.strand < right.strand;
            });
  matches.erase(
      std::unique(matches.begin(), matches.end(),
                  [](const MotifMatch &left, const MotifMatch &right) {
                    return left.chr == right.chr &&
                           left.position == right.position &&
                           left.matchLength == right.matchLength &&
                           left.strand == right.strand;
                  }),
      matches.end());

  motifResults[resultId] = matches;
}

void Interpreter::executeFindAlias(const IRInstruction &instr) {
  std::string resultId = instr.arg1;
  std::string alias = instr.arg2;

  std::vector<GenomicRegion> regions;
  if (motifResults.count(resultId)) {
    for (const auto &match : motifResults[resultId])
      regions.push_back(motifMatchToRegion(match, alias));
  }

  namedRegions[alias] = regions;
  resultSets[alias] = regions;
  if (debugMode) {
    std::cout << "  Stored " << regions.size() << " regions as \"" << alias
              << "\"" << std::endl;
  }
}

void Interpreter::executeExtract(const IRInstruction &instr) {
  std::string entityType = instr.arg1;
  std::string resultId = instr.arg2;
  if (debugMode) {
    std::cout << "> EXTRACT " << entityType << std::endl;
  }

  std::vector<GenomicRegion> regions = resolveEntity(entityType);

  if (!sequenceDatasets.empty()) {
    const auto &chrMap = sequenceChrMaps[activeSequenceAlias];
    for (auto &r : regions) {
      if (r.sequence.empty() && chrMap.count(r.chr)) {
        const std::string &seq = chrMap.at(r.chr).sequence;
        if (r.start < seq.size() && r.end <= seq.size()) {
          r.sequence = seq.substr(r.start, r.end - r.start);
        }
      }
    }
  }

  resultSets[resultId] = regions;
}

void Interpreter::executeFilterLength(const IRInstruction &instr) {
  std::string op = instr.arg1;
  std::string valueStr = instr.arg2;
  std::string resultId = instr.arg3;

  double threshold = 0;
  std::string unit = "BP";
  size_t spacePos = valueStr.find(' ');
  if (spacePos != std::string::npos) {
    threshold = std::atof(valueStr.substr(0, spacePos).c_str());
    unit = valueStr.substr(spacePos + 1);
  } else {
    threshold = std::atof(valueStr.c_str());
  }
  size_t thresholdBP = toBasePairs(threshold, unit);

  if (resultSets.count(resultId)) {
    auto &regions = resultSets[resultId];
    std::vector<GenomicRegion> filtered;
    for (const auto &r : regions) {
      bool pass = false;
      if (op == ">")
        pass = r.length() > thresholdBP;
      else if (op == "<")
        pass = r.length() < thresholdBP;
      else if (op == ">=")
        pass = r.length() >= thresholdBP;
      else if (op == "<=")
        pass = r.length() <= thresholdBP;
      else if (op == "=")
        pass = r.length() == thresholdBP;
      if (pass)
        filtered.push_back(r);
    }
    if (debugMode) {
      std::cout << "  WHERE LENGTH " << op << " " << valueStr << ": "
                << filtered.size() << " of " << regions.size() << " passed."
                << std::endl;
    }
    regions = filtered;
  } else if (motifResults.count(resultId)) {
    auto &matches = motifResults[resultId];
    std::vector<MotifMatch> filtered;
    for (const auto &match : matches) {
      const double length = static_cast<double>(match.matchLength);
      if (compareValues(length, op, std::to_string(thresholdBP)))
        filtered.push_back(match);
    }
    matches = std::move(filtered);
  }
}

void Interpreter::executeFilterSimilarity(const IRInstruction &instr) {
  std::string op = instr.arg1;
  std::string valueStr = instr.arg2;
  std::string resultId = instr.arg3;

  double threshold = 0;
  size_t spacePos = valueStr.find(' ');
  if (spacePos != std::string::npos) {
    threshold = std::atof(valueStr.substr(0, spacePos).c_str());
  } else {
    threshold = std::atof(valueStr.c_str());
  }

  auto &regions = resultSets[resultId];
  if (regions.size() < 2) {
    std::cout << "  WHERE SIMILARITY: Need at least 2 regions to compare."
              << std::endl;
    return;
  }

  std::string referenceSeq = regions[0].sequence;
  std::vector<GenomicRegion> filtered;
  filtered.push_back(regions[0]);

  for (size_t i = 1; i < regions.size(); i++) {
    if (regions[i].sequence.empty())
      continue;
    double sim =
        SmithWaterman::computeSimilarity(referenceSeq, regions[i].sequence);
    bool pass = false;
    if (op == ">")
      pass = sim > threshold;
    else if (op == "<")
      pass = sim < threshold;
    else if (op == ">=")
      pass = sim >= threshold;
    else if (op == "<=")
      pass = sim <= threshold;
    else if (op == "=")
      pass = (sim >= threshold - 0.5 && sim <= threshold + 0.5);
    if (pass)
      filtered.push_back(regions[i]);
  }
  if (debugMode) {
    std::cout << "  WHERE SIMILARITY " << op << " " << valueStr << ": "
              << filtered.size() << " of " << regions.size() << " passed."
              << std::endl;
  }
  regions = filtered;
}

bool Interpreter::compareValues(double left, const std::string &op,
                                const std::string &right) const {
  double rightValue = std::atof(right.c_str());
  const size_t spacePos = right.find(' ');
  if (spacePos != std::string::npos) {
    rightValue = std::atof(right.substr(0, spacePos).c_str());
    const std::string unit = right.substr(spacePos + 1);
    if (unit == "KB")
      rightValue *= 1000.0;
    else if (unit == "MB")
      rightValue *= 1000000.0;
  }

  if (op == ">") return left > rightValue;
  if (op == ">=") return left >= rightValue;
  if (op == "<") return left < rightValue;
  if (op == "<=") return left <= rightValue;
  if (op == "=" || op == "==") return std::abs(left - rightValue) < 1e-9;
  if (op == "!=") return std::abs(left - rightValue) >= 1e-9;
  return false;
}

bool Interpreter::evaluateRegionCondition(
    const std::shared_ptr<IRCondition> &condition,
    const GenomicRegion &region,
    const std::string &referenceSequence) const {
  if (!condition)
    return true;
  if (condition->kind == IRCondition::Kind::AND)
    return evaluateRegionCondition(condition->left, region, referenceSequence) &&
           evaluateRegionCondition(condition->right, region, referenceSequence);
  if (condition->kind == IRCondition::Kind::OR)
    return evaluateRegionCondition(condition->left, region, referenceSequence) ||
           evaluateRegionCondition(condition->right, region, referenceSequence);
  if (condition->kind == IRCondition::Kind::NOT)
    return !evaluateRegionCondition(condition->left, region, referenceSequence);

  if (condition->property == "LENGTH")
    return compareValues(static_cast<double>(region.length()), condition->op,
                         condition->value);
  if (condition->property == "SIMILARITY") {
    if (referenceSequence.empty() || region.sequence.empty())
      return false;
    const double similarity =
        SmithWaterman::computeSimilarity(referenceSequence, region.sequence);
    return compareValues(similarity, condition->op, condition->value);
  }
  if (condition->property == "GC_CONTENT") {
    if (region.sequence.empty())
      return false;
    size_t gc = 0;
    for (char c : region.sequence) {
      c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
      if (c == 'G' || c == 'C')
        ++gc;
    }
    const double percent =
        100.0 * static_cast<double>(gc) / region.sequence.size();
    return compareValues(percent, condition->op, condition->value);
  }
  if (condition->property == "ID" || condition->property == "NAME") {
    const std::string expected = stripQuotes(condition->value);
    if (condition->op == "=" || condition->op == "==")
      return region.name == expected;
    if (condition->op == "!=")
      return region.name != expected;
  }
  return false;
}

bool Interpreter::conditionContainsSimilarity(
    const std::shared_ptr<IRCondition> &condition) const {
  if (!condition)
    return false;
  if (condition->kind == IRCondition::Kind::SIMPLE)
    return condition->property == "SIMILARITY";
  return conditionContainsSimilarity(condition->left) ||
         conditionContainsSimilarity(condition->right);
}

bool Interpreter::evaluateReferenceEligibility(
    const std::shared_ptr<IRCondition> &condition,
    const GenomicRegion &region) const {
  if (!condition)
    return true;
  if (condition->kind == IRCondition::Kind::SIMPLE) {
    if (condition->property == "SIMILARITY")
      return true;
    return evaluateRegionCondition(condition, region, "");
  }
  if (condition->kind == IRCondition::Kind::AND) {
    return evaluateReferenceEligibility(condition->left, region) &&
           evaluateReferenceEligibility(condition->right, region);
  }
  // OR and NOT expressions involving similarity do not define a unique
  // pre-filter. They therefore use the first sequence-bearing candidate.
  if (conditionContainsSimilarity(condition))
    return true;
  return evaluateRegionCondition(condition, region, "");
}

bool Interpreter::evaluateMotifCondition(
    const std::shared_ptr<IRCondition> &condition,
    const MotifMatch &match) const {
  if (!condition)
    return true;
  if (condition->kind == IRCondition::Kind::AND)
    return evaluateMotifCondition(condition->left, match) &&
           evaluateMotifCondition(condition->right, match);
  if (condition->kind == IRCondition::Kind::OR)
    return evaluateMotifCondition(condition->left, match) ||
           evaluateMotifCondition(condition->right, match);
  if (condition->kind == IRCondition::Kind::NOT)
    return !evaluateMotifCondition(condition->left, match);
  if (condition->property == "LENGTH")
    return compareValues(static_cast<double>(match.matchLength), condition->op,
                         condition->value);
  if (condition->property == "GC_CONTENT") {
    const auto dataset = sequenceChrMaps.find(activeSequenceAlias);
    if (dataset == sequenceChrMaps.end())
      return false;
    const auto sequence = dataset->second.find(match.chr);
    if (sequence == dataset->second.end() || match.matchLength == 0 ||
        match.position + match.matchLength >
            sequence->second.sequence.size()) {
      return false;
    }
    size_t gc = 0;
    for (size_t i = match.position;
         i < match.position + match.matchLength; ++i) {
      const char c = static_cast<char>(std::toupper(
          static_cast<unsigned char>(sequence->second.sequence[i])));
      if (c == 'G' || c == 'C')
        ++gc;
    }
    return compareValues(100.0 * static_cast<double>(gc) / match.matchLength,
                         condition->op, condition->value);
  }
  return false;
}

bool Interpreter::evaluateGCCondition(
    const std::shared_ptr<IRCondition> &condition,
    const GCWindow &window) const {
  if (!condition)
    return true;
  if (condition->kind == IRCondition::Kind::AND)
    return evaluateGCCondition(condition->left, window) &&
           evaluateGCCondition(condition->right, window);
  if (condition->kind == IRCondition::Kind::OR)
    return evaluateGCCondition(condition->left, window) ||
           evaluateGCCondition(condition->right, window);
  if (condition->kind == IRCondition::Kind::NOT)
    return !evaluateGCCondition(condition->left, window);
  return condition->property == "GC_CONTENT" &&
         compareValues(window.gcPercent, condition->op, condition->value);
}

void Interpreter::executeFilterCondition(const IRInstruction &instr) {
  const std::string &resultId = instr.arg1;
  if (gcResults.count(resultId)) {
    if (!conditionUsesOnly(instr.condition, {"GC_CONTENT"})) {
      reportRuntimeError(
          "GC profiles can only be filtered by GC_CONTENT.");
      return;
    }
    auto &windows = gcResults[resultId];
    std::vector<GCWindow> filtered;
    for (const auto &window : windows) {
      if (evaluateGCCondition(instr.condition, window))
        filtered.push_back(window);
    }
    windows = std::move(filtered);
  } else if (resultSets.count(resultId)) {
    if (!conditionUsesOnly(
            instr.condition,
            {"LENGTH", "SIMILARITY", "GC_CONTENT", "ID", "NAME"})) {
      reportRuntimeError("Unsupported condition for a genomic region set.");
      return;
    }
    auto &regions = resultSets[resultId];
    std::string referenceSequence;
    std::set<std::string> explicitReferences;
    collectSimilarityReferences(instr.condition, explicitReferences);
    if (explicitReferences.size() > 1) {
      reportRuntimeError(
          "A filter cannot use more than one SIMILARITY reference alias.");
      return;
    }
    if (!explicitReferences.empty()) {
      const std::string &referenceAlias = *explicitReferences.begin();
      const auto referenceSet = resultSets.find(referenceAlias);
      if (referenceSet == resultSets.end()) {
        reportRuntimeError("Similarity reference alias '" + referenceAlias +
                           "' is not available.");
        return;
      }
      if (referenceSet->second.size() != 1) {
        reportRuntimeError(
            "Similarity reference alias '" + referenceAlias +
            "' must contain exactly one region; found " +
            std::to_string(referenceSet->second.size()) + ".");
        return;
      }
      GenomicRegion reference = referenceSet->second.front();
      if (reference.sequence.empty()) {
        const auto dataset = sequenceChrMaps.find(activeSequenceAlias);
        if (dataset != sequenceChrMaps.end()) {
          const auto chromosome = dataset->second.find(reference.chr);
          if (chromosome != dataset->second.end() &&
              reference.start < chromosome->second.sequence.size() &&
              reference.end <= chromosome->second.sequence.size()) {
            reference.sequence = chromosome->second.sequence.substr(
                reference.start, reference.end - reference.start);
          }
        }
      }
      if (reference.sequence.empty()) {
        reportRuntimeError("Similarity reference alias '" + referenceAlias +
                           "' has no resolvable sequence.");
        return;
      }
      referenceSequence = reference.sequence;
    } else if (conditionContainsSimilarity(instr.condition)) {
      if (debugMode) {
        std::cout
            << "  Note: implicit SIMILARITY reference uses the first eligible "
               "sequence. Prefer 'SIMILARITY TO alias' for reproducibility."
            << std::endl;
      }
      for (const auto &region : regions) {
        if (!region.sequence.empty() &&
            evaluateReferenceEligibility(instr.condition, region)) {
          referenceSequence = region.sequence;
          break;
        }
      }
    }
    if (referenceSequence.empty() && explicitReferences.empty()) {
      for (const auto &region : regions) {
        if (!region.sequence.empty()) {
          referenceSequence = region.sequence;
          break;
        }
      }
    }
    std::vector<GenomicRegion> filtered;
    std::vector<unsigned char> keep(regions.size(), 0);
    const bool parallel =
        conditionContainsSimilarity(instr.condition) && regions.size() > 1;
    if (parallel) {
      const size_t workerCount = std::min(
          regions.size(),
          std::max<size_t>(1, std::thread::hardware_concurrency()));
      std::atomic<size_t> next(0);
      std::vector<std::thread> workers;
      workers.reserve(workerCount);
      for (size_t worker = 0; worker < workerCount; ++worker) {
        workers.emplace_back([&, referenceSequence] {
          while (true) {
            const size_t index = next.fetch_add(1);
            if (index >= regions.size())
              break;
            keep[index] = evaluateRegionCondition(
                              instr.condition, regions[index],
                              referenceSequence)
                              ? 1
                              : 0;
          }
        });
      }
      for (auto &worker : workers)
        worker.join();
    } else {
      for (size_t index = 0; index < regions.size(); ++index) {
        keep[index] =
            evaluateRegionCondition(instr.condition, regions[index],
                                    referenceSequence)
                ? 1
                : 0;
      }
    }
    for (size_t index = 0; index < regions.size(); ++index) {
      if (keep[index])
        filtered.push_back(regions[index]);
    }
    regions = std::move(filtered);
  } else if (motifResults.count(resultId)) {
    if (!conditionUsesOnly(instr.condition, {"LENGTH", "GC_CONTENT"})) {
      reportRuntimeError(
          "Motif matches can only be filtered by LENGTH or GC_CONTENT.");
      return;
    }
    auto &matches = motifResults[resultId];
    std::vector<MotifMatch> filtered;
    for (const auto &match : matches) {
      if (evaluateMotifCondition(instr.condition, match))
        filtered.push_back(match);
    }
    matches = std::move(filtered);
  }
}

void Interpreter::executeSetOp(const IRInstruction &instr) {
  std::string entity1 = instr.arg1;
  std::string entity2 = instr.arg2;
  std::string resultId = instr.arg3;

  auto regions1 = resolveEntity(entity1);
  auto regions2 = resolveEntity(entity2);

  if (!sequenceDatasets.empty()) {
    const auto &chrMap = sequenceChrMaps[activeSequenceAlias];
    for (auto &r : regions1) {
      if (r.sequence.empty() && chrMap.count(r.chr)) {
        const std::string &seq = chrMap.at(r.chr).sequence;
        if (r.start < seq.size() && r.end <= seq.size())
          r.sequence = seq.substr(r.start, r.end - r.start);
      }
    }
    for (auto &r : regions2) {
      if (r.sequence.empty() && chrMap.count(r.chr)) {
        const std::string &seq = chrMap.at(r.chr).sequence;
        if (r.start < seq.size() && r.end <= seq.size())
          r.sequence = seq.substr(r.start, r.end - r.start);
      }
    }
  }

  std::vector<GenomicRegion> result;
  if (instr.opcode == IROpCode::SET_INTERSECT) {
    if (debugMode) {
      std::cout << "> INTERSECT " << entity1 << " AND " << entity2 << std::endl;
    }
    result = SetOperations::intersect(regions1, regions2);
  } else if (instr.opcode == IROpCode::SET_UNION) {
    if (debugMode) {
      std::cout << "> UNION " << entity1 << " AND " << entity2 << std::endl;
    }
    result = SetOperations::unite(regions1, regions2);
  } else {
    if (debugMode) {
      std::cout << "> EXCEPT " << entity1 << " FROM " << entity2 << std::endl;
    }
    result = SetOperations::except(regions1, regions2);
  }

  resultSets[resultId] = result;
}

void Interpreter::executePrint(const IRInstruction &instr) {
  if (!debugMode && currentPrintIndex < lastPrintIndex) {
    currentPrintIndex++;
    return;
  }
  currentPrintIndex++;

  std::string resultId = instr.arg1;
  std::string type = instr.arg2;

  if (type == "FIND") {
    if (resultSets.count(resultId) && !resultSets[resultId].empty()) {
      printRegions(resultSets[resultId]);
    } else if (motifResults.count(resultId)) {
      printMotifMatches(motifResults[resultId]);
    }
  } else {
    if (resultSets.count(resultId)) {
      printRegions(resultSets[resultId]);
    }
  }
  std::cout << std::endl;
}

void Interpreter::executeLoadMatrix(const IRInstruction &instr) {
  std::string filename = stripQuotes(instr.arg1);
  std::string alias = instr.arg2;
  if (debugMode) {
    std::cout << "> LOAD MATRIX \"" << filename << "\" AS " << alias
              << std::endl;
  }

  PWMatrix pwm = PWMScanner::loadJASPAR(filename);
  if (pwm.length == 0) {
    reportRuntimeError("Failed to load PWM from " + filename + ".");
    return;
  }

  PSSM pssm = PWMScanner::computePSSM(pwm);
  if (pssm.length == 0) {
    reportRuntimeError("Invalid matrix values or unequal row lengths in " +
                       filename + ".");
    return;
  }
  loadedMatrices[alias] = pwm;
  loadedPSSMs[alias] = pssm;

  if (debugMode) {
    std::cout << "  Loaded PWM \"" << pssm.name << "\" (" << pssm.length
              << " positions, score range: " << pssm.minScore << " to "
              << pssm.maxScore << ")" << std::endl;
  }
}

void Interpreter::executeScanOptStrand(const IRInstruction &instr) {
  currentScan.strandFilter = instr.arg1;
}

void Interpreter::executeScanOptThreshold(const IRInstruction &instr) {
  std::string valueStr = instr.arg1;
  
  size_t spacePos = valueStr.find(' ');
  if (spacePos != std::string::npos) {
    currentScan.threshold = std::atof(valueStr.substr(0, spacePos).c_str());
  } else {
    currentScan.threshold = std::atof(valueStr.c_str());
  }
}

void Interpreter::executeScanExec(const IRInstruction &instr) {
  std::string matrixAlias = instr.arg1;
  std::string resultId = instr.arg2;
  std::string target = instr.arg3;

  if (currentScan.threshold < 0.0) {
    currentScan.threshold = 75.0;
  }

  if (debugMode) {
    std::cout << "> SCAN " << matrixAlias;
    if (!target.empty()) {
      std::cout << " IN " << target;
    }
    if (!currentScan.strandFilter.empty()) {
      std::cout << " STRAND " << currentScan.strandFilter;
    }
    std::cout << " THRESHOLD " << currentScan.threshold << "%";
    std::cout << std::endl;
  }

  if (sequenceDatasets.empty()) {
    reportRuntimeError("No sequence loaded.");
    return;
  }

  if (!loadedPSSMs.count(matrixAlias)) {
    reportRuntimeError("Matrix '" + matrixAlias + "' not loaded.");
    return;
  }

  const PSSM &pssm = loadedPSSMs[matrixAlias];

  bool searchPos = (currentScan.strandFilter != "NEGATIVE");
  bool searchNeg = (currentScan.strandFilter != "POSITIVE");

  std::vector<MotifMatch> matches;
  size_t scannedUnits = 0;
  if (target.empty()) {
    const auto &targetDataset = sequenceDatasets.count(activeSequenceAlias)
                                    ? sequenceDatasets[activeSequenceAlias]
                                    : sequenceDatasets.begin()->second;
    std::vector<std::future<std::vector<MotifMatch>>> futures;
    for (const auto &seqRec : targetDataset) {
      const FastaRecord *record = &seqRec;
      const std::string chrId = seqRec.sequenceId;
      const double threshold = currentScan.threshold;
      const PSSM *matrix = &pssm;

      futures.push_back(std::async(
          std::launch::async,
          [record, matrix, threshold, chrId, searchPos, searchNeg]() {
            return PWMScanner::scan(record->sequence, *matrix, threshold,
                                    chrId, searchPos, searchNeg);
          }));
    }
    for (auto &future : futures) {
      auto chromosomeMatches = future.get();
      matches.insert(matches.end(), chromosomeMatches.begin(),
                     chromosomeMatches.end());
    }
    scannedUnits = targetDataset.size();
  } else {
    const auto chromosomeMap = sequenceChrMaps.find(activeSequenceAlias);
    if (chromosomeMap == sequenceChrMaps.end()) {
      reportRuntimeError("SCAN IN requires an active FASTA dataset.");
      currentScan = ScanContext();
      return;
    }

    const std::vector<GenomicRegion> targetRegions = resolveEntity(target);
    scannedUnits = targetRegions.size();
    for (const auto &source : targetRegions) {
      const auto chromosome = chromosomeMap->second.find(source.chr);
      if (chromosome == chromosomeMap->second.end()) {
        reportRuntimeError("SCAN target sequence identifier '" + source.chr +
                           "' is not present in the active FASTA dataset.");
        currentScan = ScanContext();
        return;
      }
      if (source.start > source.end ||
          source.end > chromosome->second.sequence.size()) {
        reportRuntimeError("SCAN target region '" + source.name +
                           "' lies outside chromosome '" + source.chr + "'.");
        currentScan = ScanContext();
        return;
      }

      const std::string sequence = chromosome->second.sequence.substr(
          source.start, source.end - source.start);
      auto regionMatches = PWMScanner::scan(
          sequence, pssm, currentScan.threshold, source.chr, searchPos,
          searchNeg);
      for (auto &match : regionMatches) {
        const size_t localPosition = match.position;
        match.position = source.start + localPosition;
        match.evidence.hasSourceRegion = true;
        match.evidence.sourceRegionName = source.name;
        match.evidence.sourceRegionType = source.type;
        match.evidence.sourceRegionStart = source.start;
        match.evidence.sourceRegionEnd = source.end;
        match.evidence.relativeStart =
            source.strand == "-"
                ? source.end - (match.position + match.matchLength)
                : localPosition;
      }
      matches.insert(matches.end(), regionMatches.begin(),
                     regionMatches.end());
    }
  }

  for (auto &match : matches)
    match.evidence.matrixAlias = matrixAlias;

  std::sort(matches.begin(), matches.end(),
            [](const MotifMatch &left, const MotifMatch &right) {
              if (left.chr != right.chr)
                return left.chr < right.chr;
              if (left.position != right.position)
                return left.position < right.position;
              if (left.strand != right.strand)
                return left.strand < right.strand;
              return left.evidence.sourceRegionName <
                     right.evidence.sourceRegionName;
            });

  motifResults[resultId] = matches;

  if (debugMode) {
    std::cout << "  PWM scan found " << matches.size() << " site(s) above "
              << currentScan.threshold << "% threshold across "
              << scannedUnits
              << (target.empty() ? " chromosome(s)." : " source region(s).")
              << std::endl;
  }

  currentScan = ScanContext();
}

void Interpreter::executeScanAlias(const IRInstruction &instr) {
  std::string resultId = instr.arg1;
  std::string alias = instr.arg2;

  std::vector<GenomicRegion> regions;
  if (motifResults.count(resultId)) {
    for (const auto &match : motifResults[resultId])
      regions.push_back(motifMatchToRegion(match, alias));
  }

  namedRegions[alias] = regions;
  resultSets[alias] = regions;
  if (debugMode) {
    std::cout << "  Stored " << regions.size() << " regions as \"" << alias
              << "\"" << std::endl;
  }
}

void Interpreter::executeResultAlias(const IRInstruction &instr) {
  const std::string &resultId = instr.arg1;
  const std::string &alias = instr.arg2;

  if (gcResults.count(resultId)) {
    if (alias != resultId) {
      gcResults[alias] = std::move(gcResults[resultId]);
      gcResults.erase(resultId);
    }
    return;
  }
  if (resultSets.count(resultId)) {
    if (alias != resultId) {
      resultSets[alias] = std::move(resultSets[resultId]);
      resultSets.erase(resultId);
    }
    namedRegions[alias] = resultSets[alias];
    return;
  }

  if (!motifResults.count(resultId))
    return;

  std::vector<GenomicRegion> regions;
  for (const auto &match : motifResults[resultId])
    regions.push_back(motifMatchToRegion(match, alias));
  resultSets[alias] = regions;
  namedRegions[alias] = std::move(regions);
}

void Interpreter::executeAnalyzeGC(const IRInstruction &instr) {
  std::string windowSizeStr = instr.arg1;
  std::string resultId = instr.arg2;
  
  size_t windowSize = 100;
  if (!windowSizeStr.empty()) {
    size_t spacePos = windowSizeStr.find(' ');
    if (spacePos != std::string::npos) {
      windowSize = toBasePairs(
          std::stod(windowSizeStr.substr(0, spacePos)),
          windowSizeStr.substr(spacePos + 1));
    } else {
      windowSize = toBasePairs(std::stod(windowSizeStr), "BP");
    }
  }

  const auto &targetDataset = sequenceDatasets.count(activeSequenceAlias)
                                  ? sequenceDatasets[activeSequenceAlias]
                                  : (sequenceDatasets.empty() ? std::vector<FastaRecord>{} : sequenceDatasets.begin()->second);
      
  std::vector<GCWindow> allWindows;
  for (const auto &seqRec : targetDataset) {
    auto windows =
        GCAnalyzer::gcContentWindowed(seqRec.sequence, windowSize, windowSize);
    for (auto &window : windows)
      window.chr = seqRec.sequenceId;
    allWindows.insert(allWindows.end(), windows.begin(), windows.end());
  }

  gcResults[resultId] = allWindows;
  
  if (debugMode) {
    std::cout << "  Computed GC profile with " << allWindows.size() << " windows for \"" << resultId << "\"" << std::endl;
  }
}

void Interpreter::executeAnalyzeCpG(const IRInstruction &instr) {
  std::string resultId = instr.arg2;
  
  const auto &targetDataset = sequenceDatasets.count(activeSequenceAlias)
                                  ? sequenceDatasets[activeSequenceAlias]
                                  : (sequenceDatasets.empty() ? std::vector<FastaRecord>{} : sequenceDatasets.begin()->second);
      
  std::vector<GenomicRegion> allIslands;
  for (const auto &seqRec : targetDataset) {
    auto islands = GCAnalyzer::findCpGIslands(seqRec.sequence, seqRec.sequenceId);
    allIslands.insert(allIslands.end(), islands.begin(), islands.end());
  }

  resultSets[resultId] = allIslands;
  
  if (debugMode) {
    std::cout << "  Found " << allIslands.size() << " CpG islands for \"" << resultId << "\"" << std::endl;
  }
}


void Interpreter::dumpResultsJSON() const {
  std::ofstream out(".cisql_results.json");
  if (!out.is_open()) return;

  out << "{\n"
      << "  \"metadata\": {\n"
      << "    \"coordinateSystem\": \"zero-based-half-open\",\n"
      << "    \"sequenceDataset\": \"" << jsonEscape(activeSequenceAlias)
      << "\",\n"
      << "    \"annotationDataset\": \"" << jsonEscape(activeAnnotationAlias)
      << "\"\n"
      << "  },\n"
      << "  \"resultSets\": {\n";
  bool firstSet = true;
  for (const auto &pair : resultSets) {
    if (!firstSet) out << ",\n";
    firstSet = false;
    out << "    \"" << jsonEscape(pair.first) << "\": [\n";
    
    bool firstRegion = true;
    for (const auto &r : pair.second) {
      if (!firstRegion) out << ",\n";
      firstRegion = false;
      out << "      {\n"
          << "        \"chr\": \"" << jsonEscape(r.chr) << "\",\n"
          << "        \"start\": " << r.start << ",\n"
          << "        \"end\": " << r.end << ",\n"
          << "        \"strand\": \"" << jsonEscape(r.strand) << "\",\n"
          << "        \"type\": \"" << jsonEscape(r.type) << "\",\n"
          << "        \"name\": \"" << jsonEscape(r.name) << "\",\n"
          << "        \"sequence\": \"" << jsonEscape(r.sequence) << "\"";
      if (r.motifEvidence.present) {
        out << ",\n        \"motifEvidence\": {\n"
            << "          \"matrixAlias\": \""
            << jsonEscape(r.motifEvidence.matrixAlias) << "\",\n"
            << "          \"matrixId\": \""
            << jsonEscape(r.motifEvidence.matrixId) << "\",\n"
            << "          \"matrixName\": \""
            << jsonEscape(r.motifEvidence.matrixName) << "\",\n"
            << "          \"matrixSource\": \""
            << jsonEscape(r.motifEvidence.matrixSource) << "\",\n"
            << "          \"rawScore\": "
            << r.motifEvidence.rawScore << ",\n"
            << "          \"scorePercent\": "
            << r.motifEvidence.scorePercent;
        if (r.motifEvidence.hasSourceRegion) {
          out << ",\n          \"sourceRegion\": {\n"
              << "            \"name\": \""
              << jsonEscape(r.motifEvidence.sourceRegionName) << "\",\n"
              << "            \"type\": \""
              << jsonEscape(r.motifEvidence.sourceRegionType) << "\",\n"
              << "            \"start\": "
              << r.motifEvidence.sourceRegionStart << ",\n"
              << "            \"end\": "
              << r.motifEvidence.sourceRegionEnd << ",\n"
              << "            \"relativeStart\": "
              << r.motifEvidence.relativeStart << "\n"
              << "          }";
        }
        out << "\n        }";
      }
      out << "\n      }";
    }
    out << "\n    ]";
  }
  out << "\n  }";

  if (!gcResults.empty()) {
    out << ",\n  \"gcProfiles\": {\n";
    bool firstProfile = true;
    for (const auto &pair : gcResults) {
      if (!firstProfile) out << ",\n";
      firstProfile = false;
      out << "    \"" << jsonEscape(pair.first) << "\": [\n";
      
      bool firstWindow = true;
      for (const auto &w : pair.second) {
        if (!firstWindow) out << ",\n";
        firstWindow = false;
        out << "      {\"chr\": \"" << jsonEscape(w.chr)
            << "\", \"pos\": " << w.position << ", \"gc\": "
            << w.gcPercent << "}";
      }
      out << "\n    ]";
    }
    out << "\n  }";
  }
  
  out << "\n}\n";
}

bool Interpreter::evaluateGlobalCondition(
    const std::shared_ptr<IRCondition> &condition) {
  if (!condition)
    return false;
  if (condition->kind == IRCondition::Kind::AND)
    return evaluateGlobalCondition(condition->left) &&
           evaluateGlobalCondition(condition->right);
  if (condition->kind == IRCondition::Kind::OR)
    return evaluateGlobalCondition(condition->left) ||
           evaluateGlobalCondition(condition->right);
  if (condition->kind == IRCondition::Kind::NOT)
    return !evaluateGlobalCondition(condition->left);

  const std::string &prop = condition->property;
  double leftVal = 0.0;
  if (prop == "GC_CONTENT") {
    size_t totalBases = 0;
    size_t gcBases = 0;
    const auto active = sequenceDatasets.find(activeSequenceAlias);
    if (active != sequenceDatasets.end()) {
      for (const auto &rec : active->second) {
        totalBases += rec.sequence.length();
        for (char c : rec.sequence) {
          if (c == 'G' || c == 'C' || c == 'g' || c == 'c') gcBases++;
        }
      }
    }
    leftVal = totalBases > 0 ? (100.0 * gcBases / totalBases) : 0.0;
  } else if (resultSets.count(prop)) {
    leftVal = static_cast<double>(resultSets[prop].size());
  } else if (motifResults.count(prop)) {
    leftVal = static_cast<double>(motifResults[prop].size());
  } else {
    leftVal = std::atof(prop.c_str());
  }

  return compareValues(leftVal, condition->op, condition->value);
}

void Interpreter::execute(const std::vector<IRInstruction> &program,
                          bool debug) {
  debugMode = debug;
  runtimeError = false;
  currentPrintIndex = 0;
  currentFind = FindContext();
  currentScan = ScanContext();

  lastPrintIndex = -1;
  int printIdx = 0;
  for (const auto &instr : program) {
    if (instr.opcode == IROpCode::PRINT_RESULTS) {
      lastPrintIndex = printIdx;
      printIdx++;
    }
  }

  for (size_t i = 0; i < program.size(); ++i) {
    const auto &instr = program[i];

    if (instr.opcode == IROpCode::IF_BEGIN) {
      if (!conditionUsesOnly(instr.condition, {"GC_CONTENT"})) {
        reportRuntimeError(
            "IF conditions currently support GC_CONTENT only.");
      }
      bool cond = evaluateGlobalCondition(instr.condition);
      if (debugMode) {
        std::cout << "> IF condition -> " << (cond ? "TRUE" : "FALSE")
                  << std::endl;
      }
      if (!cond) {
        int depth = 1;
        while (i + 1 < program.size()) {
          i++;
          if (program[i].opcode == IROpCode::IF_BEGIN) depth++;
          else if (program[i].opcode == IROpCode::IF_END) {
            depth--;
            if (depth == 0) break;
          } else if (program[i].opcode == IROpCode::IF_ELSE && depth == 1) {
            break;
          }
        }
      }
      continue;
    }

    if (instr.opcode == IROpCode::IF_ELSE) {
      int depth = 1;
      while (i + 1 < program.size()) {
        i++;
        if (program[i].opcode == IROpCode::IF_BEGIN) depth++;
        else if (program[i].opcode == IROpCode::IF_END) {
          depth--;
          if (depth == 0) break;
        }
      }
      continue;
    }

    if (instr.opcode == IROpCode::IF_END) {
      continue;
    }

    switch (instr.opcode) {
    case IROpCode::LOAD_SEQ:
      executeLoadSeq(instr);
      break;
    case IROpCode::LOAD_ANNOT:
      executeLoadAnnot(instr);
      break;
    case IROpCode::USE_SEQUENCE:
    case IROpCode::USE_ANNOTATION:
      executeUseDataset(instr);
      break;
    case IROpCode::EXPORT_RESULTS:
      executeExport(instr);
      break;
    case IROpCode::DEFINE_PROMOTERS:
      executeDefinePromoters(instr);
      break;
    case IROpCode::FIND_MOTIF:
      executeFindMotif(instr);
      break;
    case IROpCode::FIND_OPT_WITHIN:
      executeFindOptWithin(instr);
      break;
    case IROpCode::FIND_OPT_STRAND:
      executeFindOptStrand(instr);
      break;
    case IROpCode::FIND_OPT_CHR:
      executeFindOptChr(instr);
      break;
    case IROpCode::FIND_EXEC:
      executeFindExec(instr);
      break;
    case IROpCode::FIND_ALIAS:
      executeFindAlias(instr);
      break;
    case IROpCode::EXTRACT:
      executeExtract(instr);
      break;
    case IROpCode::FILTER_LENGTH:
      executeFilterLength(instr);
      break;
    case IROpCode::FILTER_SIMILARITY:
      executeFilterSimilarity(instr);
      break;
    case IROpCode::FILTER_CONDITION:
      executeFilterCondition(instr);
      break;
    case IROpCode::SET_INTERSECT:
    case IROpCode::SET_UNION:
    case IROpCode::SET_EXCEPT:
      executeSetOp(instr);
      break;
    case IROpCode::PRINT_RESULTS:
      executePrint(instr);
      break;
    case IROpCode::LOAD_MATRIX:
      executeLoadMatrix(instr);
      break;
    case IROpCode::SCAN_OPT_STRAND:
      executeScanOptStrand(instr);
      break;
    case IROpCode::SCAN_OPT_THRESHOLD:
      executeScanOptThreshold(instr);
      break;
    case IROpCode::SCAN_EXEC:
      executeScanExec(instr);
      break;
    case IROpCode::SCAN_ALIAS:
      executeScanAlias(instr);
      break;
    case IROpCode::RESULT_ALIAS:
      executeResultAlias(instr);
      break;
    case IROpCode::ANALYZE_GC:
      executeAnalyzeGC(instr);
      break;
    case IROpCode::ANALYZE_CPG:
      executeAnalyzeCpG(instr);
      break;
    default:
      break;
    }
  }

  
  if (!runtimeError)
    dumpResultsJSON();
}
