#ifndef SEMANTIC_ANALYZER_H
#define SEMANTIC_ANALYZER_H

#include "../frontend/AST.h"
#include "SymbolTable.h"
#include <string>

class SemanticAnalyzer : public ASTVisitor {
private:
  SymbolTable &symbolTable;
  bool hasError;
  bool annotationLoaded;

  void reportError(const std::string &message);
  double parseValue(const std::string &val);

public:
  SemanticAnalyzer(SymbolTable &symTab);

  void analyze(ProgramNode *node);
  bool hadError() const { return hasError; }

  void visit(SimpleConditionNode *node) override;
  void visit(BinaryConditionNode *node) override;
  void visit(NotConditionNode *node) override;
  void visit(LoadStmtNode *node) override;
  void visit(FindOptNode *node) override;
  void visit(FindStmtNode *node) override;
  void visit(ExtractStmtNode *node) override;
  void visit(SetOpStmtNode *node) override;
  void visit(ProgramNode *node) override;
};

#endif
