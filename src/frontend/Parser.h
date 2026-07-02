#ifndef PARSER_H
#define PARSER_H

#include "AST.h"
#include "Token.h"
#include <memory>
#include <string>
#include <vector>

class Parser {
private:
  std::vector<Token> tokens;
  size_t current;
  bool hasError;

  const Token &peek() const;
  const Token &previous() const;
  bool isAtEnd() const;
  Token advance();
  bool check(TokenType type) const;
  bool match(TokenType type);
  void consume(TokenType type, const std::string &message);

  void synchronize();
  void reportError(const Token &token, const std::string &message);

  std::unique_ptr<StatementNode> parseStatement();
  std::unique_ptr<LoadStmtNode> parseLoad();
  std::unique_ptr<FindStmtNode> parseFind();
  std::unique_ptr<ExtractStmtNode> parseExtract();
  std::unique_ptr<SetOpStmtNode> parseSetOperation();
  std::unique_ptr<ScanStmtNode> parseScan();

  std::unique_ptr<ConditionNode> parseWhereClause();
  std::unique_ptr<ConditionNode> parseCondition();
  std::unique_ptr<ConditionNode>
  parseConditionPrime(std::unique_ptr<ConditionNode> left);
  std::unique_ptr<ConditionNode> parseTerm();
  std::unique_ptr<ConditionNode>
  parseTermPrime(std::unique_ptr<ConditionNode> left);
  std::unique_ptr<ConditionNode> parseFactor();
  std::unique_ptr<SimpleConditionNode> parseSimpleCondition();

public:
  Parser(const std::vector<Token> &tokenList);
  std::unique_ptr<ProgramNode> parse();
  bool hadError() const { return hasError; }
};

#endif
