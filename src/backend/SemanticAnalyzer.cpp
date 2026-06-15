#include "SemanticAnalyzer.h"
#include <fstream>
#include <iostream>

SemanticAnalyzer::SemanticAnalyzer(SymbolTable& symTab)
    : symbolTable(symTab), hasError(false), annotationLoaded(false) {}

void SemanticAnalyzer::reportError(const std::string& message) {
  std::cerr << "Semantic Error: " << message << std::endl;
  hasError = true;
}

double SemanticAnalyzer::parseValue(const std::string& val) {
  try {
    return std::stod(val);
  } catch (...) {
    return 0.0;
  }
}

void SemanticAnalyzer::analyze(ProgramNode* node) {
  if (node) {
    node->accept(*this);
  }
}

void SemanticAnalyzer::visit(ProgramNode* node) {
  for (const auto& stmt : node->statements) {
    stmt->accept(*this);
  }
}

void SemanticAnalyzer::visit(LoadStmtNode* node) {
  std::string actualFilename = node->filename;
  if (actualFilename.size() >= 2 && actualFilename.front() == '"' && actualFilename.back() == '"') {
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
    } else if (node->loadType == "ANNOTATION") {
      symbolTable.insert(node->alias, "ANNOTATION_DATA");
      annotationLoaded = true;
    }
  }
}

void SemanticAnalyzer::visit(FindStmtNode* node) {
  for (const auto& opt : node->opts) {
    opt->accept(*this);
  }
}

void SemanticAnalyzer::visit(FindOptNode* node) {
  if (node->type == "WITHIN") {
    double distance = parseValue(node->value1);
    if (distance < 0) {
      reportError("Distance in WITHIN option cannot be negative.");
    }
    if (!node->value4.empty() && !annotationLoaded) {
      reportError("WITHIN ... FROM " + node->value4 + " requires annotation data. Use: LOAD ANNOTATION \"file.gff3\" AS alias;");
    }
  }
}

void SemanticAnalyzer::visit(ExtractStmtNode* node) {
  if (!annotationLoaded) {
    reportError("EXTRACT " + node->entity + " requires annotation data. Use: LOAD ANNOTATION \"file.gff3\" AS alias;");
  }
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(SetOpStmtNode* node) {
  if (!annotationLoaded) {
    reportError("Set operation " + node->op + " requires annotation data. Use: LOAD ANNOTATION \"file.gff3\" AS alias;");
  }
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(BinaryConditionNode* node) {
  if (node->left) node->left->accept(*this);
  if (node->right) node->right->accept(*this);
}

void SemanticAnalyzer::visit(NotConditionNode* node) {
  if (node->condition) node->condition->accept(*this);
}

void SemanticAnalyzer::visit(SimpleConditionNode* node) {
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
  }
}
