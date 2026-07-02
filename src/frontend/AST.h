#ifndef AST_H
#define AST_H

#include "Token.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

class SimpleConditionNode;
class BinaryConditionNode;
class NotConditionNode;
class LoadStmtNode;
class FindOptNode;
class FindStmtNode;
class ExtractStmtNode;
class SetOpStmtNode;
class ScanStmtNode;
class ProgramNode;

class ASTVisitor {
public:
  virtual ~ASTVisitor() = default;
  virtual void visit(SimpleConditionNode* node) = 0;
  virtual void visit(BinaryConditionNode* node) = 0;
  virtual void visit(NotConditionNode* node) = 0;
  virtual void visit(LoadStmtNode* node) = 0;
  virtual void visit(FindOptNode* node) = 0;
  virtual void visit(FindStmtNode* node) = 0;
  virtual void visit(ExtractStmtNode* node) = 0;
  virtual void visit(SetOpStmtNode* node) = 0;
  virtual void visit(ScanStmtNode* node) = 0;
  virtual void visit(ProgramNode* node) = 0;
};

class ASTNode {
public:
  virtual ~ASTNode() = default;
  virtual void print(std::string prefix = "", bool isLast = true) const = 0;
  virtual void accept(ASTVisitor& visitor) = 0;
};

class ConditionNode : public ASTNode {};

class SimpleConditionNode : public ConditionNode {
public:
  std::string property;
  std::string op;
  std::string value;
  SimpleConditionNode(std::string prop, std::string o, std::string val)
      : property(prop), op(o), value(val) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class BinaryConditionNode : public ConditionNode {
public:
  std::string op;
  std::unique_ptr<ConditionNode> left;
  std::unique_ptr<ConditionNode> right;
  BinaryConditionNode(std::string o, std::unique_ptr<ConditionNode> l,
                      std::unique_ptr<ConditionNode> r)
      : op(o), left(std::move(l)), right(std::move(r)) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class NotConditionNode : public ConditionNode {
public:
  std::unique_ptr<ConditionNode> condition;
  NotConditionNode(std::unique_ptr<ConditionNode> cond)
      : condition(std::move(cond)) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class StatementNode : public ASTNode {};

class LoadStmtNode : public StatementNode {
public:
  std::string loadType;
  std::string filename;
  std::string alias;
  LoadStmtNode(std::string lt, std::string f, std::string a)
      : loadType(lt), filename(f), alias(a) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class FindOptNode : public ASTNode {
public:
  std::string type;
  std::string value1;
  std::string value2;
  std::string value3;
  std::string value4;
  std::string value5;
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class FindStmtNode : public StatementNode {
public:
  std::string motif;
  std::vector<std::unique_ptr<FindOptNode>> opts;
  std::string alias;
  std::unique_ptr<ConditionNode> whereClause;
  FindStmtNode(std::string m, std::string a = "")
      : motif(m), alias(a) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class ExtractStmtNode : public StatementNode {
public:
  std::string entity;
  std::unique_ptr<ConditionNode> whereClause;
  ExtractStmtNode(std::string e, std::unique_ptr<ConditionNode> w)
      : entity(e), whereClause(std::move(w)) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class SetOpStmtNode : public StatementNode {
public:
  std::string op;
  std::string entity1;
  std::string entity2;
  std::unique_ptr<ConditionNode> whereClause;
  SetOpStmtNode(std::string o, std::string e1, std::string e2,
                std::unique_ptr<ConditionNode> w)
      : op(o), entity1(e1), entity2(e2), whereClause(std::move(w)) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class ScanStmtNode : public StatementNode {
public:
  std::string matrixAlias;    // alias of the loaded PWM
  std::string strandFilter;   // POSITIVE, NEGATIVE, or "" (both)
  std::string threshold;      // e.g. "80.0 %"
  std::string alias;          // AS alias
  std::unique_ptr<ConditionNode> whereClause;
  ScanStmtNode(std::string ma, std::string sf, std::string th,
               std::string a, std::unique_ptr<ConditionNode> w)
      : matrixAlias(ma), strandFilter(sf), threshold(th),
        alias(a), whereClause(std::move(w)) {}
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

class ProgramNode : public ASTNode {
public:
  std::vector<std::unique_ptr<StatementNode>> statements;
  void print(std::string prefix = "", bool isLast = true) const override;
  void accept(ASTVisitor& visitor) override;
};

#endif
