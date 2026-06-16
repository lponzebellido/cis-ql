#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "../bioinfo/FastaReader.h"
#include "../bioinfo/GFFReader.h"
#include "../bioinfo/GenomicRegion.h"
#include "../bioinfo/MotifFinder.h"
#include "../bioinfo/SetOperations.h"
#include "../bioinfo/SmithWaterman.h"
#include "IRGenerator.h"
#include <string>
#include <unordered_map>
#include <vector>

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
  std::unordered_map<std::string, std::vector<GenomicRegion>> namedRegions;

  FindContext currentFind;
  std::string activeSequenceAlias;

  std::string stripQuotes(const std::string &s);
  size_t toBasePairs(size_t value, const std::string &unit);
  std::vector<GenomicRegion> resolveEntity(const std::string &entity);
  void printRegions(const std::vector<GenomicRegion> &regions,
                    int maxShow = 20);
  void printMotifMatches(const std::vector<MotifMatch> &matches,
                         const std::string &pattern, int maxShow = 20);

  void executeLoadSeq(const IRInstruction &instr);
  void executeLoadAnnot(const IRInstruction &instr);
  void executeFindMotif(const IRInstruction &instr);
  void executeFindOptWithin(const IRInstruction &instr);
  void executeFindOptStrand(const IRInstruction &instr);
  void executeFindOptChr(const IRInstruction &instr);
  void executeFindExec(const IRInstruction &instr);
  void executeFindAlias(const IRInstruction &instr);
  void executeExtract(const IRInstruction &instr);
  void executeFilterLength(const IRInstruction &instr);
  void executeFilterSimilarity(const IRInstruction &instr);
  void executeSetOp(const IRInstruction &instr);
  void executePrint(const IRInstruction &instr);

public:
  void execute(const std::vector<IRInstruction> &program);
};

#endif
