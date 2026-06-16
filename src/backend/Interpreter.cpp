#include "Interpreter.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <set>

static const std::set<std::string> BUILTIN_ENTITIES = {
    "GENE", "PROMOTER", "ENHANCER", "EXON", "INTRON",
    "UTR", "TSS", "CDS", "REGION"};

static bool isBuiltinEntity(const std::string &name) {
  return BUILTIN_ENTITIES.count(name) > 0;
}

std::string Interpreter::stripQuotes(const std::string &s) {
  if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
    return s.substr(1, s.size() - 2);
  }
  return s;
}

size_t Interpreter::toBasePairs(size_t value, const std::string &unit) {
  if (unit == "KB")
    return value * 1000;
  if (unit == "MB")
    return value * 1000000;
  return value;
}

void Interpreter::printMotifMatches(const std::vector<MotifMatch> &matches,
                                    const std::string &pattern, int maxShow) {
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
              << "  strand:" << m.strand << "  context: ..." << m.context
              << "..." << std::endl;
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

std::vector<GenomicRegion> Interpreter::resolveEntity(const std::string &entity) {
  if (namedRegions.count(entity)) {
    return namedRegions[entity];
  }
  return GFFReader::filterByType(annotations, entity);
}

void Interpreter::executeLoadSeq(const IRInstruction &instr) {
  std::string filename = stripQuotes(instr.arg1);
  std::string alias = instr.arg2;
  std::cout << "> LOAD SEQUENCE \"" << filename << "\" AS " << alias
            << std::endl;

  auto records = FastaReader::read(filename);
  if (records.empty()) {
    std::cerr << "  Error: No sequences found in " << filename << std::endl;
    return;
  }

  sequences[alias] = records[0];
  activeSequenceAlias = alias;
  std::cout << "  Loaded " << records[0].sequence.size() << " base pairs from "
            << filename << " (" << records[0].sequenceId << ")" << std::endl;
  if (records.size() > 1) {
    std::cout << "  (File contains " << records.size()
              << " sequences, using first: " << records[0].sequenceId << ")"
              << std::endl;
  }
}

void Interpreter::executeLoadAnnot(const IRInstruction &instr) {
  std::string filename = stripQuotes(instr.arg1);
  std::string alias = instr.arg2;
  std::cout << "> LOAD ANNOTATION \"" << filename << "\" AS " << alias
            << std::endl;

  annotations = GFFReader::read(filename);
  std::cout << "  Loaded " << annotations.size() << " annotation features from "
            << filename << std::endl;

  std::unordered_map<std::string, int> typeCounts;
  for (const auto &a : annotations) {
    typeCounts[a.type]++;
  }
  for (const auto &pair : typeCounts) {
    std::cout << "    " << pair.first << ": " << pair.second << std::endl;
  }
}

void Interpreter::executeFindMotif(const IRInstruction &instr) {
  currentFind = FindContext();
  currentFind.pattern = stripQuotes(instr.arg1);
  currentFind.hasWithin = false;
}

void Interpreter::executeFindOptWithin(const IRInstruction &instr) {
  currentFind.hasWithin = true;
  currentFind.withinDistance = std::atol(instr.arg1.c_str());
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

  if (sequences.empty()) {
    std::cerr << "  Error: No sequence loaded." << std::endl;
    return;
  }
  const FastaRecord &seq = sequences.begin()->second;

  bool searchNeg = true;
  if (currentFind.strandFilter == "POSITIVE")
    searchNeg = false;

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
      size_t windowStart, windowEnd;
      std::string effectiveStrand =
          target.strand.empty() ? "+" : target.strand;
      if (currentFind.withinDirection == "UPSTREAM") {
        if (effectiveStrand == "+") {
          windowStart = (target.start > distBP) ? target.start - distBP : 0;
          windowEnd = target.start;
        } else {
          windowStart = target.end;
          windowEnd = std::min(target.end + distBP, seq.sequence.size());
        }
      } else {
        if (effectiveStrand == "+") {
          windowStart = target.end;
          windowEnd = std::min(target.end + distBP, seq.sequence.size());
        } else {
          windowStart = (target.start > distBP) ? target.start - distBP : 0;
          windowEnd = target.start;
        }
      }
      auto windowMatches = MotifFinder::findInWindow(
          seq.sequence, currentFind.pattern, windowStart, windowEnd,
          seq.sequenceId);
      matches.insert(matches.end(), windowMatches.begin(),
                     windowMatches.end());
    }
  } else {
    matches = MotifFinder::findAll(seq.sequence, currentFind.pattern,
                                   seq.sequenceId, searchNeg);
  }

  if (currentFind.strandFilter == "NEGATIVE") {
    std::vector<MotifMatch> filtered;
    for (const auto &m : matches) {
      if (m.strand == "-")
        filtered.push_back(m);
    }
    matches = filtered;
  } else if (currentFind.strandFilter == "POSITIVE") {
    std::vector<MotifMatch> filtered;
    for (const auto &m : matches) {
      if (m.strand == "+")
        filtered.push_back(m);
    }
    matches = filtered;
  }

  motifResults[resultId] = matches;
}

void Interpreter::executeFindAlias(const IRInstruction &instr) {
  std::string resultId = instr.arg1;
  std::string alias = instr.arg2;

  std::string chrId = "";
  if (!sequences.empty()) {
    chrId = sequences.begin()->second.sequenceId;
  }
  const std::string &seqData =
      sequences.empty() ? "" : sequences.begin()->second.sequence;

  std::vector<GenomicRegion> regions;
  if (motifResults.count(resultId)) {
    for (const auto &m : motifResults[resultId]) {
      GenomicRegion r;
      r.chr = chrId;
      r.start = m.position;
      r.end = m.position + currentFind.pattern.size();
      r.strand = m.strand;
      r.type = alias;
      r.name = alias + "_" + std::to_string(m.position);
      if (!seqData.empty() && r.end <= seqData.size()) {
        r.sequence = seqData.substr(r.start, r.end - r.start);
      }
      regions.push_back(r);
    }
  }

  namedRegions[alias] = regions;
  resultSets[resultId] = regions;
  std::cout << "  Stored " << regions.size() << " regions as \"" << alias
            << "\"" << std::endl;
}

void Interpreter::executeExtract(const IRInstruction &instr) {
  std::string entityType = instr.arg1;
  std::string resultId = instr.arg2;
  std::cout << "> EXTRACT " << entityType << std::endl;

  std::vector<GenomicRegion> regions = resolveEntity(entityType);

  if (!sequences.empty()) {
    const std::string &seq = sequences.begin()->second.sequence;
    for (auto &r : regions) {
      if (r.sequence.empty() && r.start < seq.size() && r.end <= seq.size()) {
        r.sequence = seq.substr(r.start, r.end - r.start);
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
  size_t thresholdBP = toBasePairs((size_t)threshold, unit);

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
    std::cout << "  WHERE LENGTH " << op << " " << valueStr << ": "
              << filtered.size() << " of " << regions.size() << " passed."
              << std::endl;
    regions = filtered;
  } else if (motifResults.count(resultId)) {
    std::cout << "  WHERE LENGTH " << op << " " << valueStr
              << ": filtering motif matches by pattern length." << std::endl;
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
  std::cout << "  WHERE SIMILARITY " << op << " " << valueStr << ": "
            << filtered.size() << " of " << regions.size() << " passed."
            << std::endl;
  regions = filtered;
}

void Interpreter::executeSetOp(const IRInstruction &instr) {
  std::string entity1 = instr.arg1;
  std::string entity2 = instr.arg2;
  std::string resultId = instr.arg3;

  auto regions1 = resolveEntity(entity1);
  auto regions2 = resolveEntity(entity2);

  if (!sequences.empty()) {
    const std::string &seq = sequences.begin()->second.sequence;
    for (auto &r : regions1) {
      if (r.sequence.empty() && r.start < seq.size() && r.end <= seq.size())
        r.sequence = seq.substr(r.start, r.end - r.start);
    }
    for (auto &r : regions2) {
      if (r.sequence.empty() && r.start < seq.size() && r.end <= seq.size())
        r.sequence = seq.substr(r.start, r.end - r.start);
    }
  }

  std::vector<GenomicRegion> result;
  if (instr.opcode == IROpCode::SET_INTERSECT) {
    std::cout << "> INTERSECT " << entity1 << " AND " << entity2 << std::endl;
    result = SetOperations::intersect(regions1, regions2);
  } else if (instr.opcode == IROpCode::SET_UNION) {
    std::cout << "> UNION " << entity1 << " AND " << entity2 << std::endl;
    result = SetOperations::unite(regions1, regions2);
  } else {
    std::cout << "> EXCEPT " << entity1 << " FROM " << entity2 << std::endl;
    result = SetOperations::except(regions1, regions2);
  }

  resultSets[resultId] = result;
}

void Interpreter::executePrint(const IRInstruction &instr) {
  std::string resultId = instr.arg1;
  std::string type = instr.arg2;

  if (type == "FIND") {
    if (resultSets.count(resultId) && !resultSets[resultId].empty()) {
      printRegions(resultSets[resultId]);
    } else if (motifResults.count(resultId)) {
      printMotifMatches(motifResults[resultId], currentFind.pattern);
    }
  } else {
    if (resultSets.count(resultId)) {
      printRegions(resultSets[resultId]);
    }
  }
  std::cout << std::endl;
}

void Interpreter::execute(const std::vector<IRInstruction> &program) {
  for (const auto &instr : program) {
    switch (instr.opcode) {
    case IROpCode::LOAD_SEQ:
      executeLoadSeq(instr);
      break;
    case IROpCode::LOAD_ANNOT:
      executeLoadAnnot(instr);
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
    case IROpCode::SET_INTERSECT:
    case IROpCode::SET_UNION:
    case IROpCode::SET_EXCEPT:
      executeSetOp(instr);
      break;
    case IROpCode::PRINT_RESULTS:
      executePrint(instr);
      break;
    }
  }
}
