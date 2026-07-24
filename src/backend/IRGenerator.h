#ifndef IR_GENERATOR_H
#define IR_GENERATOR_H

#include "../frontend/AST.h"
#include <memory>
#include <string>
#include <vector>

struct IRCondition {
  enum class Kind { SIMPLE, AND, OR, NOT };

  Kind kind = Kind::SIMPLE;
  std::string property;
  std::string reference;
  std::string op;
  std::string value;
  std::shared_ptr<IRCondition> left;
  std::shared_ptr<IRCondition> right;
};

enum class IROpCode {
  LOAD_SEQ,
  LOAD_ANNOT,
  USE_SEQUENCE,
  USE_ANNOTATION,
  EXPORT_RESULTS,
  FIND_MOTIF,
  FIND_OPT_WITHIN,
  FIND_OPT_STRAND,
  FIND_OPT_CHR,
  FIND_EXEC,
  FIND_ALIAS,
  EXTRACT,
  FILTER_LENGTH,
  FILTER_SIMILARITY,
  FILTER_CONDITION,
  SET_INTERSECT,
  SET_UNION,
  SET_EXCEPT,
  PRINT_RESULTS,
  LOAD_MATRIX,
  SCAN_EXEC,
  SCAN_OPT_STRAND,
  SCAN_OPT_THRESHOLD,
  SCAN_ALIAS,
  RESULT_ALIAS,
  ANALYZE_GC,
  ANALYZE_CPG,
  IF_BEGIN,
  IF_ELSE,
  IF_END
};

struct IRInstruction {
  IROpCode opcode;
  std::string arg1;
  std::string arg2;
  std::string arg3;
  std::string arg4;
  std::string arg5;
  std::shared_ptr<IRCondition> condition;
};

std::string irOpcodeToString(IROpCode op);

class IRGenerator : public ASTVisitor {
private:
  std::vector<IRInstruction> instructions;
  int tempCounter;
  std::string currentTemp;

  std::string newTemp(const std::string& prefix = "t");
  std::shared_ptr<IRCondition> lowerCondition(const ConditionNode *node) const;
  void emitFilter(const ConditionNode *node, const std::string &resultId);

public:
  IRGenerator();

  const std::vector<IRInstruction> &generate(ProgramNode *node);
  const std::vector<IRInstruction> &getInstructions() const {
    return instructions;
  }

  void visit(SimpleConditionNode *node) override;
  void visit(BinaryConditionNode *node) override;
  void visit(NotConditionNode *node) override;
  void visit(LoadStmtNode *node) override;
  void visit(UseStmtNode *node) override;
  void visit(ExportStmtNode *node) override;
  void visit(FindOptNode *node) override;
  void visit(FindStmtNode *node) override;
  void visit(ExtractStmtNode *node) override;
  void visit(SetOpStmtNode *node) override;
  void visit(ScanStmtNode *node) override;
  void visit(AnalyzeStmtNode *node) override;
  void visit(IfStmtNode *node) override;
  void visit(ForeachStmtNode *node) override;
  void visit(ProgramNode *node) override;
};

#endif
