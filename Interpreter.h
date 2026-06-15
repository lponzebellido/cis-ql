#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "IRGenerator.h"
#include "FastaReader.h"
#include "GFFReader.h"
#include "GenomicRegion.h"
#include "MotifFinder.h"
#include "SmithWaterman.h"
#include "SetOperations.h"
#include <string>
#include <vector>
#include <unordered_map>

struct FindContext {
  std::string pattern;
  std::string strandFilter;
  std::string chrFilter;
  bool hasWithin;
  size_t withinDistance;
  std::string withinUnit;
  std::string withinDirection;
  std::string withinEntity;
  std::string withinTarget;
};

class Interpreter {
private:
  std::unordered_map<std::string, FastaRecord> sequences;
  std::vector<GenomicRegion> annotations;
  std::unordered_map<std::string, std::vector<GenomicRegion>> resultSets;
  std::unordered_map<std::string, std::vector<MotifMatch>> motifResults;

  FindContext currentFind;
  std::string activeSequenceAlias;

  std::string stripQuotes(const std::string& s);
  size_t toBasePairs(size_t value, const std::string& unit);
  void printRegions(const std::vector<GenomicRegion>& regions, int maxShow = 20);
  void printMotifMatches(const std::vector<MotifMatch>& matches, const std::string& pattern, int maxShow = 20);

  void executeLoadSeq(const IRInstruction& instr);
  void executeLoadAnnot(const IRInstruction& instr);
  void executeFindMotif(const IRInstruction& instr);
  void executeFindOptWithin(const IRInstruction& instr);
  void executeFindOptStrand(const IRInstruction& instr);
  void executeFindOptChr(const IRInstruction& instr);
  void executeFindExec(const IRInstruction& instr);
  void executeExtract(const IRInstruction& instr);
  void executeFilterLength(const IRInstruction& instr);
  void executeFilterSimilarity(const IRInstruction& instr);
  void executeSetOp(const IRInstruction& instr);
  void executePrint(const IRInstruction& instr);

public:
  void execute(const std::vector<IRInstruction>& program);
};

#endif
