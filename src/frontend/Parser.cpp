#include "Parser.h"
#include <iostream>
#include <stdexcept>

Parser::Parser(const std::vector<Token> &tokenList)
    : tokens(tokenList), current(0), hasError(false) {}

const Token &Parser::peek() const { return tokens[current]; }

const Token &Parser::previous() const { return tokens[current - 1]; }

bool Parser::isAtEnd() const { return peek().type == TokenType::END_OF_FILE; }

Token Parser::advance() {
  if (!isAtEnd())
    current++;
  return previous();
}

bool Parser::check(TokenType type) const {
  if (isAtEnd())
    return false;
  return peek().type == type;
}

bool Parser::match(TokenType type) {
  if (check(type)) {
    advance();
    return true;
  }
  return false;
}

void Parser::consume(TokenType type, const std::string &message) {
  if (check(type)) {
    advance();
    return;
  }
  reportError(peek(), message);
  throw std::runtime_error("Parse error");
}

void Parser::reportError(const Token &token, const std::string &message) {
  hasError = true;
  std::cerr << "Syntax Error at L" << token.line << ":C" << token.column
            << " at '" << token.lexeme << "': " << message << std::endl;
}

void Parser::synchronize() {
  advance();
  while (!isAtEnd()) {
    if (previous().type == TokenType::SEMICOLON)
      return;
    switch (peek().type) {
    case TokenType::LOAD:
    case TokenType::FIND:
    case TokenType::EXTRACT:
    case TokenType::INTERSECT:
    case TokenType::UNION:
    case TokenType::EXCEPT:
      return;
    default:
      break;
    }
    advance();
  }
}

std::unique_ptr<ProgramNode> Parser::parse() {
  auto program = std::unique_ptr<ProgramNode>(new ProgramNode());
  while (!isAtEnd()) {
    try {
      auto stmt = parseStatement();
      if (stmt) {
        program->statements.push_back(std::move(stmt));
      }
    } catch (std::runtime_error &) {
      synchronize();
    }
  }
  return program;
}

std::unique_ptr<StatementNode> Parser::parseStatement() {
  if (match(TokenType::LOAD))
    return parseLoad();
  if (match(TokenType::FIND))
    return parseFind();
  if (match(TokenType::EXTRACT))
    return parseExtract();
  if (match(TokenType::INTERSECT) || match(TokenType::UNION) ||
      match(TokenType::EXCEPT)) {
    current--;
    return parseSetOperation();
  }
  if (match(TokenType::SCAN))
    return parseScan();

  if (match(TokenType::ERROR_TOKEN)) {
    return nullptr;
  }

  reportError(peek(), "Expected start of a statement (LOAD, FIND, EXTRACT, "
                      "INTERSECT, UNION, EXCEPT, SCAN)");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<LoadStmtNode> Parser::parseLoad() {
  std::string loadType;
  if (match(TokenType::SEQUENCE)) {
    loadType = "SEQUENCE";
  } else if (match(TokenType::ANNOTATION)) {
    loadType = "ANNOTATION";
  } else if (match(TokenType::MATRIX)) {
    loadType = "MATRIX";
  } else {
    reportError(peek(), "Expected 'SEQUENCE', 'ANNOTATION', or 'MATRIX' after LOAD.");
    throw std::runtime_error("Parse error");
  }
  consume(TokenType::STRING, "Expected a file name (string).");
  std::string file = previous().lexeme;
  consume(TokenType::AS, "Expected 'AS' after the file name.");
  consume(TokenType::ID, "Expected an alias identifier.");
  std::string alias = previous().lexeme;
  consume(TokenType::SEMICOLON,
          "Expected ';' at the end of the LOAD statement.");
  return std::unique_ptr<LoadStmtNode>(new LoadStmtNode(loadType, file, alias));
}

std::unique_ptr<FindStmtNode> Parser::parseFind() {
  consume(TokenType::MOTIF, "Expected 'MOTIF' after FIND.");
  consume(TokenType::STRING, "Expected the motif name (string).");
  auto node =
      std::unique_ptr<FindStmtNode>(new FindStmtNode(previous().lexeme));

  while (!check(TokenType::SEMICOLON) && !check(TokenType::AS) &&
         !check(TokenType::WHERE) && !isAtEnd()) {
    auto opt = std::unique_ptr<FindOptNode>(new FindOptNode());
    if (match(TokenType::WITHIN)) {
      opt->type = "WITHIN";
      consume(TokenType::NUM, "Expected a number for WITHIN.");
      opt->value1 = previous().lexeme;

      if (match(TokenType::BP) || match(TokenType::KB) ||
          match(TokenType::MB)) {
        opt->value2 = previous().lexeme;
      }

      if (match(TokenType::UPSTREAM) || match(TokenType::DOWNSTREAM)) {
        opt->value3 = previous().lexeme;
      } else {
        reportError(peek(), "Expected UPSTREAM or DOWNSTREAM.");
        throw std::runtime_error("Parse error");
      }

      consume(TokenType::FROM, "Expected 'FROM' after the direction.");

      if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
          match(TokenType::ENHANCER) || match(TokenType::EXON) ||
          match(TokenType::INTRON) || match(TokenType::UTR) ||
          match(TokenType::TSS) || match(TokenType::CDS) ||
          match(TokenType::REGION) || match(TokenType::ID)) {
        opt->value4 = previous().lexeme;
      } else {
        reportError(peek(), "Expected a biological entity or alias.");
        throw std::runtime_error("Parse error");
      }

      if (check(TokenType::STRING)) {
        consume(TokenType::STRING, "Expected the entity name.");
        opt->value5 = previous().lexeme;
      }

    } else if (match(TokenType::STRAND)) {
      opt->type = "STRAND";
      if (match(TokenType::POSITIVE) || match(TokenType::NEGATIVE)) {
        opt->value1 = previous().lexeme;
      } else {
        reportError(peek(), "Expected POSITIVE or NEGATIVE.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::CHR)) {
      opt->type = "CHR";
      consume(TokenType::STRING, "Expected the chromosome name.");
      opt->value1 = previous().lexeme;
    } else {
      break;
    }
    node->opts.push_back(std::move(opt));
  }

  if (match(TokenType::AS)) {
    consume(TokenType::ID, "Expected an alias identifier after AS.");
    node->alias = previous().lexeme;
  }

  node->whereClause = parseWhereClause();

  consume(TokenType::SEMICOLON, "Expected ';' at the end of FIND.");
  return node;
}

std::unique_ptr<ExtractStmtNode> Parser::parseExtract() {
  if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
      match(TokenType::ENHANCER) || match(TokenType::EXON) ||
      match(TokenType::INTRON) || match(TokenType::UTR) ||
      match(TokenType::TSS) || match(TokenType::CDS) ||
      match(TokenType::REGION) || match(TokenType::ID)) {
    std::string entity = previous().lexeme;
    auto whereClause = parseWhereClause();
    consume(TokenType::SEMICOLON, "Expected ';' at the end of EXTRACT.");
    return std::unique_ptr<ExtractStmtNode>(
        new ExtractStmtNode(entity, std::move(whereClause)));
  }
  reportError(peek(), "Expected an entity or alias for EXTRACT.");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<SetOpStmtNode> Parser::parseSetOperation() {
  std::string op;
  if (match(TokenType::INTERSECT))
    op = "INTERSECT";
  else if (match(TokenType::UNION))
    op = "UNION";
  else if (match(TokenType::EXCEPT))
    op = "EXCEPT";

  if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
      match(TokenType::ENHANCER) || match(TokenType::EXON) ||
      match(TokenType::INTRON) || match(TokenType::UTR) ||
      match(TokenType::TSS) || match(TokenType::CDS) ||
      match(TokenType::REGION) || match(TokenType::ID)) {
    std::string e1 = previous().lexeme;

    std::string sepError =
        "Expected 'AND' (for INTERSECT/UNION) or 'FROM' (for EXCEPT).";
    if (op == "EXCEPT")
      consume(TokenType::FROM, sepError);
    else
      consume(TokenType::AND, sepError);

    if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
        match(TokenType::ENHANCER) || match(TokenType::EXON) ||
        match(TokenType::INTRON) || match(TokenType::UTR) ||
        match(TokenType::TSS) || match(TokenType::CDS) ||
        match(TokenType::REGION) || match(TokenType::ID)) {
      std::string e2 = previous().lexeme;
      auto whereClause = parseWhereClause();
      consume(TokenType::SEMICOLON,
              "Expected ';' at the end of the set operation.");
      return std::unique_ptr<SetOpStmtNode>(
          new SetOpStmtNode(op, e1, e2, std::move(whereClause)));
    }
  }
  reportError(peek(), "Expected an entity or alias in the set operation.");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<ScanStmtNode> Parser::parseScan() {
  
  consume(TokenType::ID, "Expected a matrix alias after SCAN.");
  std::string matrixAlias = previous().lexeme;

  std::string strandFilter;
  std::string threshold;

  
  while (!check(TokenType::SEMICOLON) && !check(TokenType::AS) &&
         !check(TokenType::WHERE) && !isAtEnd()) {
    if (match(TokenType::STRAND)) {
      if (match(TokenType::POSITIVE) || match(TokenType::NEGATIVE)) {
        strandFilter = previous().lexeme;
      } else {
        reportError(peek(), "Expected POSITIVE or NEGATIVE after STRAND.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::THRESHOLD)) {
      if (match(TokenType::NUM) || match(TokenType::FLOAT)) {
        threshold = previous().lexeme;
        if (match(TokenType::PERCENT)) {
          threshold += " %";
        }
      } else {
        reportError(peek(), "Expected a numeric value for THRESHOLD.");
        throw std::runtime_error("Parse error");
      }
    } else {
      break;
    }
  }

  
  std::string alias;
  if (match(TokenType::AS)) {
    consume(TokenType::ID, "Expected an alias identifier after AS.");
    alias = previous().lexeme;
  }

  
  auto whereClause = parseWhereClause();

  consume(TokenType::SEMICOLON, "Expected ';' at the end of SCAN.");

  return std::unique_ptr<ScanStmtNode>(
      new ScanStmtNode(matrixAlias, strandFilter, threshold,
                       alias, std::move(whereClause)));
}

std::unique_ptr<ConditionNode> Parser::parseWhereClause() {
  if (match(TokenType::WHERE)) {
    return parseCondition();
  }
  return nullptr;
}

std::unique_ptr<ConditionNode> Parser::parseCondition() {
  auto term = parseTerm();
  return parseConditionPrime(std::move(term));
}

std::unique_ptr<ConditionNode>
Parser::parseConditionPrime(std::unique_ptr<ConditionNode> left) {
  if (match(TokenType::OR)) {
    auto rightTerm = parseTerm();
    auto newLeft = std::unique_ptr<BinaryConditionNode>(
        new BinaryConditionNode("OR", std::move(left), std::move(rightTerm)));
    return parseConditionPrime(std::move(newLeft));
  }
  return left;
}

std::unique_ptr<ConditionNode> Parser::parseTerm() {
  auto factor = parseFactor();
  return parseTermPrime(std::move(factor));
}

std::unique_ptr<ConditionNode>
Parser::parseTermPrime(std::unique_ptr<ConditionNode> left) {
  if (match(TokenType::AND)) {
    auto rightFactor = parseFactor();
    auto newLeft = std::unique_ptr<BinaryConditionNode>(new BinaryConditionNode(
        "AND", std::move(left), std::move(rightFactor)));
    return parseTermPrime(std::move(newLeft));
  }
  return left;
}

std::unique_ptr<ConditionNode> Parser::parseFactor() {
  if (match(TokenType::NOT)) {
    return std::unique_ptr<NotConditionNode>(
        new NotConditionNode(parseFactor()));
  }
  if (match(TokenType::LPAREN)) {
    auto cond = parseCondition();
    consume(TokenType::RPAREN, "Expected ')' after the condition.");
    return cond;
  }
  return parseSimpleCondition();
}

std::unique_ptr<SimpleConditionNode> Parser::parseSimpleCondition() {
  std::string prop;
  if (match(TokenType::LENGTH) || match(TokenType::SIMILARITY)) {
    prop = previous().lexeme;
  } else {
    reportError(peek(), "Expected 'LENGTH' or 'SIMILARITY'.");
    throw std::runtime_error("Parse error");
  }

  std::string op;
  if (match(TokenType::GREATER) || match(TokenType::LESS) ||
      match(TokenType::GREATER_EQ) || match(TokenType::LESS_EQ) ||
      match(TokenType::ASSIGN)) {
    op = previous().lexeme;
  } else {
    reportError(peek(), "Expected a relational operator (> < >= <= =).");
    throw std::runtime_error("Parse error");
  }

  std::string val;
  if (match(TokenType::NUM) || match(TokenType::FLOAT) ||
      match(TokenType::STRING)) {
    val = previous().lexeme;
    if (match(TokenType::BP) || match(TokenType::KB) || match(TokenType::MB) ||
        match(TokenType::PERCENT)) {
      val += " " + previous().lexeme;
    }
  } else {
    reportError(peek(), "Expected a value (NUM, FLOAT, STRING).");
    throw std::runtime_error("Parse error");
  }

  return std::unique_ptr<SimpleConditionNode>(
      new SimpleConditionNode(prop, op, val));
}
