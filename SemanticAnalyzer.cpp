#include "SemanticAnalyzer.h"
#include <iostream>

SemanticAnalyzer::SemanticAnalyzer(SymbolTable& symTab)
    : symbolTable(symTab), hasError(false) {}

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
  if (symbolTable.lookup(node->alias)) {
    reportError("Alias '" + node->alias + "' is already defined.");
  } else {
    symbolTable.insert(node->alias, "GENOME_DATA");
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
  }
}

void SemanticAnalyzer::visit(ExtractStmtNode* node) {
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }
}

void SemanticAnalyzer::visit(SetOpStmtNode* node) {
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
