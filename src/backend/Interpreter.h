#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "../bioinfo/FastaReader.h"
#include "../bioinfo/GFFReader.h"
#include "../bioinfo/GenomicRegion.h"
#include "../bioinfo/MotifFinder.h"
#include "../bioinfo/SetOperations.h"
#include "../bioinfo/SmithWaterman.h"
#include "../bioinfo/PWMScanner.h"
#include "../bioinfo/RegulatoryRegions.h"
#include "../bioinfo/GCAnalyzer.h"
#include "IRGenerator.h"
#include <future>
#include <string>
#include <unordered_map>
#include <vector>

struct FindContext {
  std::string pattern;
  std::string strandFilter;
  std::string chrFilter;
  bool hasWithin = false;
  double withinDistance = 0.0;
  std::string withinUnit;
  std::string withinDirection;
  std::string withinEntity;
  std::string withinTarget;
};

struct ScanContext {
  std::string strandFilter;
  double threshold = -1.0;
};

class Interpreter {
private:
  std::unordered_map<std::string, std::vector<FastaRecord>> sequenceDatasets;
  std::unordered_map<std::string, std::unordered_map<std::string, FastaRecord>> sequenceChrMaps;
  std::unordered_map<std::string, std::vector<GenomicRegion>>
      annotationDatasets;
  std::unordered_map<std::string, std::vector<GenomicRegion>> resultSets;
  std::unordered_map<std::string, std::vector<MotifMatch>> motifResults;
  std::unordered_map<std::string, std::vector<GenomicRegion>> namedRegions;
  std::unordered_map<std::string, PWMatrix> loadedMatrices;
  std::unordered_map<std::string, PSSM> loadedPSSMs;
  std::unordered_map<std::string, std::vector<GCWindow>> gcResults;

  FindContext currentFind;
  ScanContext currentScan;
  std::string activeSequenceAlias;
  std::string activeAnnotationAlias;
  bool debugMode;
  bool runtimeError;
  int lastPrintIndex;
  int currentPrintIndex;

  std::string stripQuotes(const std::string &s) const;
  std::string jsonEscape(const std::string &s) const;
  void reportRuntimeError(const std::string &message);
  size_t toBasePairs(double value, const std::string &unit);
  std::vector<GenomicRegion> resolveEntity(const std::string &entity);
  void printRegions(const std::vector<GenomicRegion> &regions,
                    int maxShow = 20);
  void printMotifMatches(const std::vector<MotifMatch> &matches,
                         int maxShow = 20);

  void executeLoadSeq(const IRInstruction &instr);
  void executeLoadAnnot(const IRInstruction &instr);
  void executeUseDataset(const IRInstruction &instr);
  void executeExport(const IRInstruction &instr);
  void executeDefinePromoters(const IRInstruction &instr);
  void executeFindMotif(const IRInstruction &instr);
  void executeFindOptWithin(const IRInstruction &instr);
  void executeFindOptStrand(const IRInstruction &instr);
  void executeFindOptChr(const IRInstruction &instr);
  void executeFindExec(const IRInstruction &instr);
  void executeFindAlias(const IRInstruction &instr);
  void executeExtract(const IRInstruction &instr);
  void executeFilterLength(const IRInstruction &instr);
  void executeFilterSimilarity(const IRInstruction &instr);
  void executeFilterCondition(const IRInstruction &instr);
  void executeSetOp(const IRInstruction &instr);
  void executePrint(const IRInstruction &instr);
  void executeLoadMatrix(const IRInstruction &instr);
  void executeScanOptStrand(const IRInstruction &instr);
  void executeScanOptThreshold(const IRInstruction &instr);
  void executeScanExec(const IRInstruction &instr);
  void executeScanAlias(const IRInstruction &instr);
  void executeResultAlias(const IRInstruction &instr);
  void executeAnalyzeGC(const IRInstruction &instr);
  void executeAnalyzeCpG(const IRInstruction &instr);
  bool compareValues(double left, const std::string &op,
                     const std::string &right) const;
  bool evaluateRegionCondition(const std::shared_ptr<IRCondition> &condition,
                               const GenomicRegion &region,
                               const std::string &referenceSequence) const;
  bool conditionContainsSimilarity(
      const std::shared_ptr<IRCondition> &condition) const;
  bool evaluateReferenceEligibility(
      const std::shared_ptr<IRCondition> &condition,
      const GenomicRegion &region) const;
  bool evaluateMotifCondition(const std::shared_ptr<IRCondition> &condition,
                              const MotifMatch &match) const;
  bool evaluateGCCondition(const std::shared_ptr<IRCondition> &condition,
                           const GCWindow &window) const;
  bool evaluateGlobalCondition(const std::shared_ptr<IRCondition> &condition);
  
  void dumpResultsJSON() const;

public:
  void execute(const std::vector<IRInstruction> &program, bool debug = false);
  bool hadError() const { return runtimeError; }
};

#endif
