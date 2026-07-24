#include "SemanticAnalyzer.h"
#include <fstream>
#include <iostream>
#include <set>

static const std::set<std::string> BUILTIN_ENTITIES = {
    "GENE", "PROMOTER", "ENHANCER", "EXON", "INTRON",
    "UTR", "TSS", "CDS", "REGION"};

static bool isBuiltinEntity(const std::string &name) {
  return BUILTIN_ENTITIES.count(name) > 0;
}

static bool isResultAlias(const SymbolTable &symbolTable,
                          const std::string &name) {
  return symbolTable.typeOf(name) == "RESULT_SET";
}

SemanticAnalyzer::SemanticAnalyzer(SymbolTable &symTab)
    : symbolTable(symTab), hasError(false), annotationLoaded(false),
      sequenceLoaded(false) {}

void SemanticAnalyzer::reportError(const std::string &message) {
  std::cerr << "Semantic Error: " << message << std::endl;
  hasError = true;
}

double SemanticAnalyzer::parseValue(const std::string &val) {
  try {
    return std::stod(val);
  } catch (...) {
    return 0.0;
  }
}

void SemanticAnalyzer::analyze(ProgramNode *node) {
  if (node) {
    node->accept(*this);
  }
}

void SemanticAnalyzer::visit(ProgramNode *node) {
  for (const auto &stmt : node->statements) {
    stmt->accept(*this);
  }
}

void SemanticAnalyzer::visit(LoadStmtNode *node) {
  std::string actualFilename = node->filename;
  if (actualFilename.size() >= 2 && actualFilename.front() == '"' &&
      actualFilename.back() == '"') {
    actualFilename = actualFilename.substr(1, actualFilename.size() - 2);
  }

  std::ifstream file(actualFilename);
  if (!file.good()) {
    reportError("File '" + actualFilename + "' not found or cannot be opened.");
  }

  if (symbolTable.lookup(node->alias)) {
    reportError("Alias '" + node->alias + "' is already defined.");
  } else {
    if (node->loadType == "SEQUENCE") {
      symbolTable.insert(node->alias, "GENOME_DATA");
      sequenceLoaded = true;
    } else if (node->loadType == "ANNOTATION") {
      symbolTable.insert(node->alias, "ANNOTATION_DATA");
      annotationLoaded = true;
    } else if (node->loadType == "MATRIX") {
      symbolTable.insert(node->alias, "PWM_DATA");
    }
  }
}

void SemanticAnalyzer::visit(FindStmtNode *node) {
  if (!sequenceLoaded) {
    reportError("FIND requires sequence data. Use: LOAD SEQUENCE "
                "\"file.fasta\" AS alias;");
  }
  if (!node->alias.empty()) {
    if (symbolTable.lookup(node->alias)) {
      reportError("Alias '" + node->alias + "' is already defined.");
    } else {
      symbolTable.insert(node->alias, "RESULT_SET");
    }
  }
  for (const auto &opt : node->opts) {
    opt->accept(*this);
  }
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(FindOptNode *node) {
  if (node->type == "WITHIN") {
    double distance = parseValue(node->value1);
    if (distance < 0) {
      reportError("Distance in WITHIN option cannot be negative.");
    }
    if (!node->value4.empty() && isBuiltinEntity(node->value4) &&
        !annotationLoaded) {
      reportError("WITHIN ... FROM " + node->value4 +
                  " requires annotation data or a user-defined alias.");
    }
    if (!node->value4.empty() && !isBuiltinEntity(node->value4) &&
        !symbolTable.lookup(node->value4)) {
      reportError("Alias '" + node->value4 + "' is not defined.");
    } else if (!node->value4.empty() && !isBuiltinEntity(node->value4) &&
               !isResultAlias(symbolTable, node->value4)) {
      reportError("WITHIN ... FROM expects a result-set alias, but '" +
                  node->value4 + "' has type " +
                  symbolTable.typeOf(node->value4) + ".");
    }
  }
}

void SemanticAnalyzer::visit(ExtractStmtNode *node) {
  if (isBuiltinEntity(node->entity) && !annotationLoaded) {
    reportError("EXTRACT " + node->entity +
                " requires annotation data. Use: LOAD ANNOTATION \"file.gff3\" "
                "AS alias;");
  }
  if (!isBuiltinEntity(node->entity) && !symbolTable.lookup(node->entity)) {
    reportError("Alias '" + node->entity + "' is not defined.");
  } else if (!isBuiltinEntity(node->entity) &&
             !isResultAlias(symbolTable, node->entity)) {
    reportError("EXTRACT expects a result-set alias, but '" + node->entity +
                "' has type " + symbolTable.typeOf(node->entity) + ".");
  }
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(SetOpStmtNode *node) {
  if (isBuiltinEntity(node->entity1) && !annotationLoaded) {
    reportError("Set operation on " + node->entity1 +
                " requires annotation data.");
  }
  if (!isBuiltinEntity(node->entity1) &&
      !symbolTable.lookup(node->entity1)) {
    reportError("Alias '" + node->entity1 + "' is not defined.");
  } else if (!isBuiltinEntity(node->entity1) &&
             !isResultAlias(symbolTable, node->entity1)) {
    reportError("Set operations require result-set aliases; '" +
                node->entity1 + "' has type " +
                symbolTable.typeOf(node->entity1) + ".");
  }
  if (isBuiltinEntity(node->entity2) && !annotationLoaded) {
    reportError("Set operation on " + node->entity2 +
                " requires annotation data.");
  }
  if (!isBuiltinEntity(node->entity2) &&
      !symbolTable.lookup(node->entity2)) {
    reportError("Alias '" + node->entity2 + "' is not defined.");
  } else if (!isBuiltinEntity(node->entity2) &&
             !isResultAlias(symbolTable, node->entity2)) {
    reportError("Set operations require result-set aliases; '" +
                node->entity2 + "' has type " +
                symbolTable.typeOf(node->entity2) + ".");
  }
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(BinaryConditionNode *node) {
  if (node->left)
    node->left->accept(*this);
  if (node->right)
    node->right->accept(*this);
}

void SemanticAnalyzer::visit(NotConditionNode *node) {
  if (node->condition)
    node->condition->accept(*this);
}

void SemanticAnalyzer::visit(SimpleConditionNode *node) {
  static const std::set<std::string> supportedProperties = {
      "LENGTH", "SIMILARITY", "GC_CONTENT", "ID", "NAME"};
  if (!supportedProperties.count(node->property)) {
    reportError("Unsupported condition property '" + node->property + "'.");
    return;
  }
  if ((node->property == "SIMILARITY" ||
       node->property == "GC_CONTENT") &&
      !sequenceLoaded) {
    reportError(node->property + " requires sequence data.");
  }
  if (node->property == "LENGTH") {
    double length = parseValue(node->value);
    if (length < 0) {
      reportError("LENGTH cannot be negative.");
    }
  } else if (node->property == "SIMILARITY") {
    double similarity = parseValue(node->value);
    if (similarity < 0 || similarity > 100) {
      reportError("SIMILARITY must be between 0 and 100.");
    }
  } else if (node->property == "GC_CONTENT") {
    double gc = parseValue(node->value);
    if (gc < 0 || gc > 100) {
      reportError("GC_CONTENT must be between 0 and 100.");
    }
  } else if ((node->property == "ID" || node->property == "NAME") &&
             (node->value.size() < 2 || node->value.front() != '"' ||
              node->value.back() != '"')) {
    reportError(node->property + " must be compared with a string value.");
  }
}

void SemanticAnalyzer::visit(ScanStmtNode *node) {
  if (!sequenceLoaded) {
    reportError("SCAN requires sequence data.");
  }

  if (!symbolTable.lookup(node->matrixAlias)) {
    reportError("Matrix alias '" + node->matrixAlias +
                "' is not defined. Use: LOAD MATRIX \"file.pwm\" AS alias;");
  } else {
    const std::string type = symbolTable.typeOf(node->matrixAlias);
    if (type != "PWM_DATA" && type != "ITERATOR_VAR") {
      reportError("SCAN expects a matrix alias, but '" + node->matrixAlias +
                  "' has type " + type + ".");
    }
  }

  
  if (!node->threshold.empty()) {
    double threshold = parseValue(node->threshold);
    if (threshold < 0 || threshold > 100) {
      reportError("THRESHOLD must be between 0 and 100.");
    }
  }

  
  if (!node->alias.empty()) {
    if (symbolTable.lookup(node->alias)) {
      reportError("Alias '" + node->alias + "' is already defined.");
    } else {
      symbolTable.insert(node->alias, "RESULT_SET");
    }
  }

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(AnalyzeStmtNode *node) {
  if (!sequenceLoaded) {
    reportError("ANALYZE requires sequence data.");
  }
  if (!node->windowSize.empty() && parseValue(node->windowSize) <= 0.0) {
    reportError("ANALYZE WINDOW must be greater than zero.");
  }
  if (!node->alias.empty()) {
    if (symbolTable.lookup(node->alias)) {
      reportError("Alias '" + node->alias + "' is already defined.");
    } else {
      symbolTable.insert(node->alias, node->analysisType == "GC_CONTENT"
                                          ? "GC_PROFILE"
                                          : "RESULT_SET");
    }
  }

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(IfStmtNode *node) {
  if (node->condition) node->condition->accept(*this);
  for (auto &stmt : node->thenStatements) {
    if (stmt) stmt->accept(*this);
  }
  for (auto &stmt : node->elseStatements) {
    if (stmt) stmt->accept(*this);
  }
}

void SemanticAnalyzer::visit(ForeachStmtNode *node) {
  if (symbolTable.lookup(node->iteratorVar)) {
    reportError("FOREACH iterator '" + node->iteratorVar +
                "' is already defined.");
  }
  for (const auto &item : node->collection) {
    if (!symbolTable.lookup(item)) {
      reportError("FOREACH collection item '" + item +
                  "' is not a defined alias.");
    } else if (symbolTable.typeOf(item) != "PWM_DATA") {
      reportError("FOREACH currently supports matrix aliases; '" + item +
                  "' has type " + symbolTable.typeOf(item) + ".");
    }
  }
  symbolTable.insert(node->iteratorVar, "ITERATOR_VAR");
  for (auto &stmt : node->bodyStatements) {
    if (stmt) stmt->accept(*this);
  }
}
