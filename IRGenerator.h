#ifndef IR_GENERATOR_H
#define IR_GENERATOR_H

#include "AST.h"
#include <string>
#include <vector>

enum class IROpCode {
  LOAD_SEQ,
  LOAD_ANNOT,
  FIND_MOTIF,
  FIND_OPT_WITHIN,
  FIND_OPT_STRAND,
  FIND_OPT_CHR,
  FIND_EXEC,
  EXTRACT,
  FILTER_LENGTH,
  FILTER_SIMILARITY,
  SET_INTERSECT,
  SET_UNION,
  SET_EXCEPT,
  PRINT_RESULTS
};

struct IRInstruction {
  IROpCode opcode;
  std::string arg1;
  std::string arg2;
  std::string arg3;
  std::string arg4;
  std::string arg5;
};

std::string irOpcodeToString(IROpCode op);

class IRGenerator : public ASTVisitor {
private:
  std::vector<IRInstruction> instructions;
  int tempCounter;
  std::string currentTemp;

  std::string newTemp();

public:
  IRGenerator();

  const std::vector<IRInstruction>& generate(ProgramNode* node);
  const std::vector<IRInstruction>& getInstructions() const { return instructions; }

  void visit(SimpleConditionNode* node) override;
  void visit(BinaryConditionNode* node) override;
  void visit(NotConditionNode* node) override;
  void visit(LoadStmtNode* node) override;
  void visit(FindOptNode* node) override;
  void visit(FindStmtNode* node) override;
  void visit(ExtractStmtNode* node) override;
  void visit(SetOpStmtNode* node) override;
  void visit(ProgramNode* node) override;
};

#endif
