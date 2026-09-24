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
    case TokenType::USE:
    case TokenType::EXPORT:
    case TokenType::DEFINE:
    case TokenType::FIND:
    case TokenType::EXTRACT:
    case TokenType::INTERSECT:
    case TokenType::UNION:
    case TokenType::EXCEPT:
    case TokenType::OVERLAPS:
    case TokenType::NEAR:
    case TokenType::CONSENSUS:
    case TokenType::COUNT:
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
  if (match(TokenType::USE))
    return parseUse();
  if (match(TokenType::EXPORT))
    return parseExport();
  if (match(TokenType::DEFINE)) {
    if (check(TokenType::PROMOTERS))
      return parseDefinePromoters();
    if (check(TokenType::MODULE))
      return parseDefineModule();
    reportError(peek(), "Expected 'PROMOTERS' or 'MODULE' after DEFINE.");
    throw std::runtime_error("Parse error");
  }
  if (match(TokenType::FIND))
    return parseFind();
  if (match(TokenType::EXTRACT))
    return parseExtract();
  if (match(TokenType::COUNT))
    return parseCount();
  if (match(TokenType::INTERSECT) || match(TokenType::UNION) ||
      match(TokenType::EXCEPT) || match(TokenType::OVERLAPS) ||
      match(TokenType::NEAR) || match(TokenType::CONSENSUS)) {
    current--;
    return parseSetOperation();
  }
  if (match(TokenType::SCAN))
    return parseScan();
  if (match(TokenType::ANALYZE))
    return parseAnalyze();
  if (match(TokenType::IF))
    return parseIf();
  if (match(TokenType::FOREACH))
    return parseForeach();

  if (match(TokenType::ERROR_TOKEN)) {
    return nullptr;
  }

  reportError(peek(), "Expected start of a statement (LOAD, USE, EXPORT, FIND, "
                      "EXTRACT, INTERSECT, UNION, EXCEPT, OVERLAPS, NEAR, "
                      "CONSENSUS, SCAN, COUNT, ANALYZE, IF, FOREACH, DEFINE)");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<DefinePromotersStmtNode> Parser::parseDefinePromoters() {
  consume(TokenType::PROMOTERS, "Expected 'PROMOTERS' after DEFINE.");
  consume(TokenType::OF, "Expected 'OF' after DEFINE PROMOTERS.");

  std::string source;
  if (match(TokenType::GENE) || match(TokenType::TSS) ||
      match(TokenType::ID)) {
    source = previous().lexeme;
  } else {
    reportError(peek(), "Expected GENE, TSS, or a result-set alias after OF.");
    throw std::runtime_error("Parse error");
  }

  consume(TokenType::FROM, "Expected 'FROM TSS' after the promoter source.");
  consume(TokenType::TSS, "Expected 'TSS' after FROM.");
  consume(TokenType::UPSTREAM, "Expected 'UPSTREAM' after FROM TSS.");
  if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
    reportError(peek(), "Expected an upstream distance.");
    throw std::runtime_error("Parse error");
  }
  const std::string upstreamValue = previous().lexeme;
  if (!match(TokenType::BP) && !match(TokenType::KB) &&
      !match(TokenType::MB)) {
    reportError(peek(), "Expected BP, KB, or MB after the upstream distance.");
    throw std::runtime_error("Parse error");
  }
  const std::string upstreamUnit = previous().lexeme;

  consume(TokenType::DOWNSTREAM,
          "Expected 'DOWNSTREAM' after the upstream window.");
  if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
    reportError(peek(), "Expected a downstream distance.");
    throw std::runtime_error("Parse error");
  }
  const std::string downstreamValue = previous().lexeme;
  if (!match(TokenType::BP) && !match(TokenType::KB) &&
      !match(TokenType::MB)) {
    reportError(peek(),
                "Expected BP, KB, or MB after the downstream distance.");
    throw std::runtime_error("Parse error");
  }
  const std::string downstreamUnit = previous().lexeme;

  consume(TokenType::AS, "Expected 'AS' after the promoter window.");
  consume(TokenType::ID, "Expected an alias after AS.");
  const std::string alias = previous().lexeme;
  consume(TokenType::SEMICOLON,
          "Expected ';' at the end of DEFINE PROMOTERS.");

  return std::unique_ptr<DefinePromotersStmtNode>(
      new DefinePromotersStmtNode(source, upstreamValue, upstreamUnit,
                                  downstreamValue, downstreamUnit, alias));
}

std::unique_ptr<DefineModuleStmtNode> Parser::parseDefineModule() {
  consume(TokenType::MODULE, "Expected 'MODULE' after DEFINE.");
  consume(TokenType::FROM, "Expected 'FROM' after DEFINE MODULE.");

  if (!match(TokenType::ID)) {
    reportError(peek(), "Expected a motif-hit alias after FROM.");
    throw std::runtime_error("Parse error");
  }
  const std::string firstSet = previous().lexeme;

  consume(TokenType::WITH, "Expected 'WITH' after the first motif-hit alias.");
  if (!match(TokenType::ID)) {
    reportError(peek(), "Expected a second motif-hit alias after WITH.");
    throw std::runtime_error("Parse error");
  }
  const std::string secondSet = previous().lexeme;

  consume(TokenType::SPACING,
          "Expected 'SPACING' after the module member aliases.");
  if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
    reportError(peek(), "Expected a minimum module spacing.");
    throw std::runtime_error("Parse error");
  }
  const std::string minimumSpacingValue = previous().lexeme;
  if (!match(TokenType::BP) && !match(TokenType::KB) &&
      !match(TokenType::MB)) {
    reportError(peek(),
                "Expected BP, KB, or MB after the minimum module spacing.");
    throw std::runtime_error("Parse error");
  }
  const std::string minimumSpacingUnit = previous().lexeme;

  consume(TokenType::TO,
          "Expected 'TO' between the minimum and maximum module spacing.");
  if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
    reportError(peek(), "Expected a maximum module spacing.");
    throw std::runtime_error("Parse error");
  }
  const std::string maximumSpacingValue = previous().lexeme;
  if (!match(TokenType::BP) && !match(TokenType::KB) &&
      !match(TokenType::MB)) {
    reportError(peek(),
                "Expected BP, KB, or MB after the maximum module spacing.");
    throw std::runtime_error("Parse error");
  }
  const std::string maximumSpacingUnit = previous().lexeme;

  consume(TokenType::ORDER, "Expected 'ORDER' after the spacing range.");
  if (!match(TokenType::ANY) && !match(TokenType::AS_WRITTEN)) {
    reportError(peek(), "Expected ANY or AS_WRITTEN after ORDER.");
    throw std::runtime_error("Parse error");
  }
  const std::string orderPolicy = previous().lexeme;

  consume(TokenType::ORIENTATION,
          "Expected 'ORIENTATION' after the order policy.");
  if (!match(TokenType::ANY) && !match(TokenType::SAME) &&
      !match(TokenType::OPPOSITE)) {
    reportError(peek(),
                "Expected ANY, SAME, or OPPOSITE after ORIENTATION.");
    throw std::runtime_error("Parse error");
  }
  const std::string orientationPolicy = previous().lexeme;

  consume(TokenType::AS, "Expected 'AS' after the orientation policy.");
  consume(TokenType::ID, "Expected a module result alias after AS.");
  const std::string alias = previous().lexeme;
  consume(TokenType::SEMICOLON,
          "Expected ';' at the end of DEFINE MODULE.");

  return std::unique_ptr<DefineModuleStmtNode>(new DefineModuleStmtNode(
      firstSet, secondSet, minimumSpacingValue, minimumSpacingUnit,
      maximumSpacingValue, maximumSpacingUnit, orderPolicy,
      orientationPolicy, alias));
}

std::unique_ptr<UseStmtNode> Parser::parseUse() {
  std::string datasetType;
  if (match(TokenType::SEQUENCE)) {
    datasetType = "SEQUENCE";
  } else if (match(TokenType::ANNOTATION)) {
    datasetType = "ANNOTATION";
  } else {
    reportError(peek(), "Expected 'SEQUENCE' or 'ANNOTATION' after USE.");
    throw std::runtime_error("Parse error");
  }
  consume(TokenType::ID, "Expected a dataset alias after USE.");
  const std::string alias = previous().lexeme;
  consume(TokenType::SEMICOLON, "Expected ';' at the end of USE.");
  return std::unique_ptr<UseStmtNode>(new UseStmtNode(datasetType, alias));
}

std::unique_ptr<ExportStmtNode> Parser::parseExport() {
  consume(TokenType::ID, "Expected a result alias after EXPORT.");
  const std::string alias = previous().lexeme;
  consume(TokenType::TO, "Expected 'TO' after the result alias.");
  consume(TokenType::STRING, "Expected an output file name after TO.");
  const std::string filename = previous().lexeme;
  consume(TokenType::FORMAT, "Expected 'FORMAT' after the output file name.");

  std::string format;
  if (match(TokenType::BED) || match(TokenType::GFF3) ||
      match(TokenType::TSV)) {
    format = previous().lexeme;
  } else {
    reportError(peek(), "Expected BED, GFF3, or TSV after FORMAT.");
    throw std::runtime_error("Parse error");
  }
  consume(TokenType::SEMICOLON, "Expected ';' at the end of EXPORT.");
  return std::unique_ptr<ExportStmtNode>(
      new ExportStmtNode(alias, filename, format));
}

std::unique_ptr<LoadStmtNode> Parser::parseLoad() {
  std::string loadType;
  if (match(TokenType::SEQUENCE)) {
    loadType = "SEQUENCE";
  } else if (match(TokenType::ANNOTATION)) {
    loadType = "ANNOTATION";
  } else if (match(TokenType::MATRIX)) {
    loadType = "MATRIX";
  } else if (match(TokenType::TRACK)) {
    loadType = "TRACK";
  } else {
    reportError(peek(), "Expected 'SEQUENCE', 'ANNOTATION', 'MATRIX', or 'TRACK' after LOAD.");
    throw std::runtime_error("Parse error");
  }
  consume(TokenType::STRING, "Expected a file name (string).");
  std::string file = previous().lexeme;
  std::string format;
  std::string evidenceClass;
  std::string assay;
  std::string sample;
  std::string condition;
  std::string replicate;
  std::string control;
  if (loadType == "TRACK") {
    consume(TokenType::FORMAT, "Expected 'FORMAT' after the track file name.");
    if (match(TokenType::BED) || match(TokenType::NARROWPEAK)) {
      format = previous().lexeme;
    } else {
      reportError(peek(), "Expected BED or NARROWPEAK after FORMAT.");
      throw std::runtime_error("Parse error");
    }
    consume(TokenType::EVIDENCE,
            "Expected 'EVIDENCE' after the track input format.");
    if (match(TokenType::ACCESSIBILITY) || match(TokenType::BINDING) ||
        match(TokenType::OTHER)) {
      evidenceClass = previous().lexeme;
    } else {
      reportError(peek(),
                  "Expected ACCESSIBILITY, BINDING, or OTHER after EVIDENCE.");
      throw std::runtime_error("Parse error");
    }
    bool hasAssay = false;
    bool hasSample = false;
    bool hasCondition = false;
    bool hasReplicate = false;
    bool hasControl = false;
    while (check(TokenType::ASSAY) || check(TokenType::SAMPLE) ||
           check(TokenType::CONDITION) || check(TokenType::REPLICATE) ||
           check(TokenType::CONTROL)) {
      const Token metadata = advance();
      bool *seen = nullptr;
      std::string *destination = nullptr;
      if (metadata.type == TokenType::ASSAY) {
        seen = &hasAssay;
        destination = &assay;
      } else if (metadata.type == TokenType::SAMPLE) {
        seen = &hasSample;
        destination = &sample;
      } else if (metadata.type == TokenType::CONDITION) {
        seen = &hasCondition;
        destination = &condition;
      } else if (metadata.type == TokenType::REPLICATE) {
        seen = &hasReplicate;
        destination = &replicate;
      } else {
        seen = &hasControl;
        destination = &control;
      }
      if (*seen) {
        reportError(metadata, "Duplicate " + metadata.lexeme +
                                  " clause in LOAD TRACK.");
        throw std::runtime_error("Parse error");
      }
      *seen = true;
      consume(TokenType::STRING,
              "Expected a string after '" + metadata.lexeme + "'.");
      *destination = previous().lexeme;
    }
  }
  consume(TokenType::AS, "Expected 'AS' after the file name.");
  consume(TokenType::ID, "Expected an alias identifier.");
  std::string alias = previous().lexeme;
  consume(TokenType::SEMICOLON,
          "Expected ';' at the end of the LOAD statement.");
  return std::unique_ptr<LoadStmtNode>(
      new LoadStmtNode(loadType, file, format, evidenceClass, assay, sample,
                       condition, replicate, control, alias));
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
      if (match(TokenType::NUM) || match(TokenType::FLOAT)) {
        opt->value1 = previous().lexeme;
      } else {
        reportError(peek(), "Expected a number for WITHIN.");
        throw std::runtime_error("Parse error");
      }

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
    std::string alias;
    if (match(TokenType::AS)) {
      consume(TokenType::ID, "Expected an alias identifier after AS.");
      alias = previous().lexeme;
    }
    auto whereClause = parseWhereClause();
    consume(TokenType::SEMICOLON, "Expected ';' at the end of EXTRACT.");
    return std::unique_ptr<ExtractStmtNode>(
        new ExtractStmtNode(entity, alias, std::move(whereClause)));
  }
  reportError(peek(), "Expected an entity or alias for EXTRACT.");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<SetOpStmtNode> Parser::parseSetOperation() {
  std::string op;
  if (match(TokenType::CONSENSUS)) {
    consume(TokenType::FROM, "Expected 'FROM' after CONSENSUS.");
    consume(TokenType::LBRACKET,
            "Expected '[' before the CONSENSUS input aliases.");
    std::vector<std::string> entities;
    consume(TokenType::ID,
            "Expected at least one input alias in CONSENSUS.");
    entities.push_back(previous().lexeme);
    while (match(TokenType::COMMA)) {
      consume(TokenType::ID, "Expected an input alias after ','.");
      entities.push_back(previous().lexeme);
    }
    consume(TokenType::RBRACKET,
            "Expected ']' after the CONSENSUS input aliases.");
    consume(TokenType::ANCHOR,
            "Expected 'ANCHOR' after the CONSENSUS input list.");
    consume(TokenType::ID, "Expected an anchor alias after ANCHOR.");
    const std::string anchor = previous().lexeme;
    consume(TokenType::MIN_SUPPORT,
            "Expected 'MIN_SUPPORT' after the CONSENSUS anchor.");
    consume(TokenType::NUM,
            "Expected a whole-number support threshold after MIN_SUPPORT.");
    const std::string minimumSupport = previous().lexeme;
    std::string minimumReciprocalOverlap;
    std::string maximumSummitDistanceValue;
    std::string maximumSummitDistanceUnit;
    while (check(TokenType::MIN_RECIPROCAL_OVERLAP) ||
           check(TokenType::MAX_SUMMIT_DISTANCE)) {
      if (match(TokenType::MIN_RECIPROCAL_OVERLAP)) {
        if (!minimumReciprocalOverlap.empty()) {
          reportError(previous(),
                      "Duplicate MIN_RECIPROCAL_OVERLAP clause.");
          throw std::runtime_error("Parse error");
        }
        if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
          reportError(peek(), "Expected a percentage after "
                              "MIN_RECIPROCAL_OVERLAP.");
          throw std::runtime_error("Parse error");
        }
        minimumReciprocalOverlap = previous().lexeme;
        consume(TokenType::PERCENT,
                "Expected '%' after MIN_RECIPROCAL_OVERLAP percentage.");
      } else {
        match(TokenType::MAX_SUMMIT_DISTANCE);
        if (!maximumSummitDistanceValue.empty()) {
          reportError(previous(), "Duplicate MAX_SUMMIT_DISTANCE clause.");
          throw std::runtime_error("Parse error");
        }
        if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
          reportError(peek(),
                      "Expected a distance after MAX_SUMMIT_DISTANCE.");
          throw std::runtime_error("Parse error");
        }
        maximumSummitDistanceValue = previous().lexeme;
        if (!match(TokenType::BP) && !match(TokenType::KB) &&
            !match(TokenType::MB)) {
          reportError(peek(), "Expected BP, KB, or MB after "
                              "MAX_SUMMIT_DISTANCE.");
          throw std::runtime_error("Parse error");
        }
        maximumSummitDistanceUnit = previous().lexeme;
      }
    }
    consume(TokenType::AS, "Expected 'AS' after the CONSENSUS threshold.");
    consume(TokenType::ID, "Expected an alias identifier after AS.");
    const std::string alias = previous().lexeme;
    auto whereClause = parseWhereClause();
    consume(TokenType::SEMICOLON,
            "Expected ';' at the end of CONSENSUS.");
    auto node = std::unique_ptr<SetOpStmtNode>(new SetOpStmtNode(
        "CONSENSUS", "", "", "", "", alias, std::move(whereClause)));
    node->entities = std::move(entities);
    node->anchor = anchor;
    node->minimumSupport = minimumSupport;
    node->minimumReciprocalOverlap = minimumReciprocalOverlap;
    node->maximumSummitDistanceValue = maximumSummitDistanceValue;
    node->maximumSummitDistanceUnit = maximumSummitDistanceUnit;
    return node;
  }
  if (match(TokenType::INTERSECT))
    op = "INTERSECT";
  else if (match(TokenType::UNION))
    op = "UNION";
  else if (match(TokenType::EXCEPT))
    op = "EXCEPT";
  else if (match(TokenType::OVERLAPS))
    op = "OVERLAPS";
  else if (match(TokenType::NEAR))
    op = "NEAR";

  if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
      match(TokenType::ENHANCER) || match(TokenType::EXON) ||
      match(TokenType::INTRON) || match(TokenType::UTR) ||
      match(TokenType::TSS) || match(TokenType::CDS) ||
      match(TokenType::REGION) || match(TokenType::ID)) {
    std::string e1 = previous().lexeme;

    std::string sepError = "Expected 'AND' (for INTERSECT/UNION), 'FROM' "
                           "(for EXCEPT), 'WITH' (for OVERLAPS), or 'TO' "
                           "(for NEAR).";
    if (op == "EXCEPT") {
      consume(TokenType::FROM, sepError);
    } else if (op == "OVERLAPS") {
      consume(TokenType::WITH, sepError);
    } else if (op == "NEAR") {
      consume(TokenType::TO, sepError);
    } else {
      consume(TokenType::AND, sepError);
    }

    if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
        match(TokenType::ENHANCER) || match(TokenType::EXON) ||
        match(TokenType::INTRON) || match(TokenType::UTR) ||
        match(TokenType::TSS) || match(TokenType::CDS) ||
        match(TokenType::REGION) || match(TokenType::ID)) {
      std::string e2 = previous().lexeme;
      std::string distanceValue;
      std::string distanceUnit;
      if (op == "NEAR") {
        consume(TokenType::WITHIN,
                "Expected 'WITHIN' after the NEAR reference entity.");
        if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
          reportError(peek(), "Expected a maximum distance after WITHIN.");
          throw std::runtime_error("Parse error");
        }
        distanceValue = previous().lexeme;
        if (!match(TokenType::BP) && !match(TokenType::KB) &&
            !match(TokenType::MB)) {
          reportError(peek(),
                      "Expected BP, KB, or MB after the NEAR distance.");
          throw std::runtime_error("Parse error");
        }
        distanceUnit = previous().lexeme;
      }
      std::string alias;
      if (match(TokenType::AS)) {
        consume(TokenType::ID, "Expected an alias identifier after AS.");
        alias = previous().lexeme;
      }
      auto whereClause = parseWhereClause();
      consume(TokenType::SEMICOLON,
              "Expected ';' at the end of the set operation.");
      return std::unique_ptr<SetOpStmtNode>(
          new SetOpStmtNode(op, e1, e2, distanceValue, distanceUnit, alias,
                            std::move(whereClause)));
    }
  }
  reportError(peek(), "Expected an entity or alias in the set operation.");
  throw std::runtime_error("Parse error");
}

std::unique_ptr<CountStmtNode> Parser::parseCount() {
  std::string countedEntity;
  if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
      match(TokenType::ENHANCER) || match(TokenType::EXON) ||
      match(TokenType::INTRON) || match(TokenType::UTR) ||
      match(TokenType::TSS) || match(TokenType::CDS) ||
      match(TokenType::REGION) || match(TokenType::ID)) {
    countedEntity = previous().lexeme;
  } else {
    reportError(peek(), "Expected an entity or alias after COUNT.");
    throw std::runtime_error("Parse error");
  }

  consume(TokenType::IN, "Expected 'IN' after the counted entity.");
  std::string containerEntity;
  if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
      match(TokenType::ENHANCER) || match(TokenType::EXON) ||
      match(TokenType::INTRON) || match(TokenType::UTR) ||
      match(TokenType::TSS) || match(TokenType::CDS) ||
      match(TokenType::REGION) || match(TokenType::ID)) {
    containerEntity = previous().lexeme;
  } else {
    reportError(peek(), "Expected a container entity or alias after IN.");
    throw std::runtime_error("Parse error");
  }

  consume(TokenType::AS, "Expected 'AS' after the container entity.");
  consume(TokenType::ID, "Expected a result alias after AS.");
  const std::string alias = previous().lexeme;
  auto whereClause = parseWhereClause();
  consume(TokenType::SEMICOLON, "Expected ';' at the end of COUNT.");
  return std::unique_ptr<CountStmtNode>(new CountStmtNode(
      countedEntity, containerEntity, alias, std::move(whereClause)));
}

std::unique_ptr<ScanStmtNode> Parser::parseScan() {
  
  consume(TokenType::ID, "Expected a matrix alias after SCAN.");
  std::string matrixAlias = previous().lexeme;

  std::string target;
  std::string strandFilter;
  std::string threshold;
  std::string significanceMetric;
  std::string significanceOperator;
  std::string significanceThreshold;
  std::string backgroundMode;
  std::string backgroundSource;

  
  while (!check(TokenType::SEMICOLON) && !check(TokenType::AS) &&
         !check(TokenType::WHERE) && !isAtEnd()) {
    if (match(TokenType::IN)) {
      if (!target.empty()) {
        reportError(previous(), "SCAN accepts only one IN target.");
        throw std::runtime_error("Parse error");
      }
      if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
          match(TokenType::ENHANCER) || match(TokenType::EXON) ||
          match(TokenType::INTRON) || match(TokenType::UTR) ||
          match(TokenType::TSS) || match(TokenType::CDS) ||
          match(TokenType::REGION) || match(TokenType::ID)) {
        target = previous().lexeme;
      } else {
        reportError(peek(), "Expected a biological entity or alias after IN.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::STRAND)) {
      if (match(TokenType::POSITIVE) || match(TokenType::NEGATIVE)) {
        strandFilter = previous().lexeme;
      } else {
        reportError(peek(), "Expected POSITIVE or NEGATIVE after STRAND.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::THRESHOLD)) {
      if (!threshold.empty()) {
        reportError(previous(), "SCAN accepts only one THRESHOLD option.");
        throw std::runtime_error("Parse error");
      }
      if (match(TokenType::NUM) || match(TokenType::FLOAT)) {
        threshold = previous().lexeme;
        consume(TokenType::PERCENT,
                "Expected '%' after the relative THRESHOLD value.");
        threshold += " %";
      } else {
        reportError(peek(), "Expected a numeric value for THRESHOLD.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::PVALUE) || match(TokenType::QVALUE)) {
      if (!significanceMetric.empty()) {
        reportError(previous(),
                    "SCAN accepts only one PVALUE or QVALUE option.");
        throw std::runtime_error("Parse error");
      }
      significanceMetric = previous().lexeme;
      if (match(TokenType::LESS) || match(TokenType::LESS_EQ)) {
        significanceOperator = previous().lexeme;
      } else {
        reportError(peek(), "Expected '<' or '<=' after " +
                                significanceMetric + ".");
        throw std::runtime_error("Parse error");
      }
      if (match(TokenType::NUM) || match(TokenType::FLOAT)) {
        significanceThreshold = previous().lexeme;
      } else {
        reportError(peek(), "Expected a numeric significance threshold.");
        throw std::runtime_error("Parse error");
      }
    } else if (match(TokenType::BACKGROUND)) {
      if (!backgroundMode.empty()) {
        reportError(previous(), "SCAN accepts only one BACKGROUND option.");
        throw std::runtime_error("Parse error");
      }
      if (match(TokenType::UNIFORM)) {
        backgroundMode = "UNIFORM";
      } else if (match(TokenType::FROM)) {
        backgroundMode = "FROM";
        if (match(TokenType::GENE) || match(TokenType::PROMOTER) ||
            match(TokenType::ENHANCER) || match(TokenType::EXON) ||
            match(TokenType::INTRON) || match(TokenType::UTR) ||
            match(TokenType::TSS) || match(TokenType::CDS) ||
            match(TokenType::REGION) || match(TokenType::ID)) {
          backgroundSource = previous().lexeme;
        } else {
          reportError(peek(),
                      "Expected a sequence or region alias after "
                      "BACKGROUND FROM.");
          throw std::runtime_error("Parse error");
        }
      } else {
        reportError(peek(), "Expected UNIFORM or FROM after BACKGROUND.");
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
      new ScanStmtNode(matrixAlias, target, strandFilter, threshold,
                       significanceMetric, significanceOperator,
                       significanceThreshold,
                       backgroundMode, backgroundSource, alias,
                       std::move(whereClause)));
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
  if (match(TokenType::LENGTH) || match(TokenType::START) ||
      match(TokenType::END) || match(TokenType::STRAND) ||
      match(TokenType::SIMILARITY) ||
      match(TokenType::GC_CONTENT) || match(TokenType::COUNT) ||
      match(TokenType::TRACK_SCORE) || match(TokenType::SIGNAL_VALUE) ||
      match(TokenType::MINUS_LOG10_PVALUE) ||
      match(TokenType::MINUS_LOG10_QVALUE) ||
      match(TokenType::EVIDENCE_CLASS) || match(TokenType::SUPPORT_COUNT) ||
      match(TokenType::ASSAY) ||
      match(TokenType::SAMPLE) || match(TokenType::CONDITION) ||
      match(TokenType::REPLICATE) || match(TokenType::CONTROL) ||
      match(TokenType::ID)) {
    prop = previous().lexeme;
  } else {
    reportError(peek(), "Expected a supported condition property.");
    throw std::runtime_error("Parse error");
  }

  std::string reference;
  if (prop == "SIMILARITY" && match(TokenType::TO)) {
    consume(TokenType::ID,
            "Expected a result-set alias after 'SIMILARITY TO'.");
    reference = previous().lexeme;
  }

  std::string modifier;
  std::string modifierValue;
  if (match(TokenType::MOD)) {
    modifier = previous().lexeme;
    if (!match(TokenType::NUM) && !match(TokenType::FLOAT)) {
      reportError(peek(), "Expected a numeric divisor after MOD.");
      throw std::runtime_error("Parse error");
    }
    modifierValue = previous().lexeme;
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
      new SimpleConditionNode(prop, op, val, reference, modifier,
                              modifierValue));
}

std::unique_ptr<AnalyzeStmtNode> Parser::parseAnalyze() {
  std::string analysisType;
  if (match(TokenType::GC_CONTENT)) {
    analysisType = "GC_CONTENT";
  } else if (match(TokenType::CPG_ISLANDS)) {
    analysisType = "CPG_ISLANDS";
  } else {
    reportError(peek(), "Expected 'GC_CONTENT' or 'CPG_ISLANDS' after ANALYZE.");
    throw std::runtime_error("Parse error");
  }

  std::string windowSize;
  if (match(TokenType::WINDOW)) {
    if (match(TokenType::NUM) || match(TokenType::FLOAT)) {
      windowSize = previous().lexeme;
      if (match(TokenType::BP) || match(TokenType::KB) || match(TokenType::MB)) {
        windowSize += " " + previous().lexeme;
      }
    } else {
      reportError(peek(), "Expected window size (number) after WINDOW.");
      throw std::runtime_error("Parse error");
    }
  }

  std::string alias;
  if (match(TokenType::AS)) {
    consume(TokenType::ID, "Expected identifier after AS.");
    alias = previous().lexeme;
  }

  auto where = parseWhereClause();

  consume(TokenType::SEMICOLON, "Expected ';' after ANALYZE statement.");

  return std::unique_ptr<AnalyzeStmtNode>(
      new AnalyzeStmtNode(analysisType, windowSize, alias, std::move(where)));
}

std::unique_ptr<IfStmtNode> Parser::parseIf() {
  auto cond = parseCondition();
  consume(TokenType::THEN, "Expected 'THEN' after IF condition.");
  std::vector<std::unique_ptr<StatementNode>> thenStmts;
  while (!check(TokenType::ELSE) && !check(TokenType::ENDIF) && !isAtEnd()) {
    auto s = parseStatement();
    if (s) thenStmts.push_back(std::move(s));
  }
  std::vector<std::unique_ptr<StatementNode>> elseStmts;
  if (match(TokenType::ELSE)) {
    while (!check(TokenType::ENDIF) && !isAtEnd()) {
      auto s = parseStatement();
      if (s) elseStmts.push_back(std::move(s));
    }
  }
  consume(TokenType::ENDIF, "Expected 'ENDIF' after IF statement.");
  match(TokenType::SEMICOLON);
  return std::unique_ptr<IfStmtNode>(new IfStmtNode(std::move(cond), std::move(thenStmts), std::move(elseStmts)));
}

std::unique_ptr<ForeachStmtNode> Parser::parseForeach() {
  consume(TokenType::ID, "Expected loop variable after FOREACH.");
  std::string var = previous().lexeme;
  consume(TokenType::IN, "Expected 'IN' after loop variable.");
  consume(TokenType::LBRACKET, "Expected '[' before FOREACH collection.");
  std::vector<std::string> coll;
  if (!check(TokenType::RBRACKET)) {
    do {
      if (match(TokenType::ID)) {
        coll.push_back(previous().lexeme);
      } else {
        reportError(peek(),
                    "Expected a matrix alias in the FOREACH collection.");
        throw std::runtime_error("Parse error");
      }
    } while (match(TokenType::COMMA));
  }
  consume(TokenType::RBRACKET, "Expected ']' after FOREACH collection.");
  consume(TokenType::DO, "Expected 'DO' after FOREACH collection.");

  std::vector<std::unique_ptr<StatementNode>> bodyStmts;
  while (!check(TokenType::ENDFOR) && !isAtEnd()) {
    auto s = parseStatement();
    if (s) bodyStmts.push_back(std::move(s));
  }
  consume(TokenType::ENDFOR, "Expected 'ENDFOR' after FOREACH body.");
  match(TokenType::SEMICOLON);
  return std::unique_ptr<ForeachStmtNode>(new ForeachStmtNode(var, coll, std::move(bodyStmts)));
}
