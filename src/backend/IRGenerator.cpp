#include "IRGenerator.h"
#include <algorithm>

std::string irOpcodeToString(IROpCode op) {
  switch (op) {
    case IROpCode::LOAD_SEQ:         return "LOAD_SEQ";
    case IROpCode::LOAD_ANNOT:       return "LOAD_ANNOT";
    case IROpCode::FIND_MOTIF:       return "FIND_MOTIF";
    case IROpCode::FIND_OPT_WITHIN:  return "FIND_OPT_WITHIN";
    case IROpCode::FIND_OPT_STRAND:  return "FIND_OPT_STRAND";
    case IROpCode::FIND_OPT_CHR:     return "FIND_OPT_CHR";
    case IROpCode::FIND_EXEC:        return "FIND_EXEC";
    case IROpCode::FIND_ALIAS:       return "FIND_ALIAS";
    case IROpCode::EXTRACT:          return "EXTRACT";
    case IROpCode::FILTER_LENGTH:    return "FILTER_LENGTH";
    case IROpCode::FILTER_SIMILARITY:return "FILTER_SIMILARITY";
    case IROpCode::SET_INTERSECT:    return "SET_INTERSECT";
    case IROpCode::SET_UNION:        return "SET_UNION";
    case IROpCode::SET_EXCEPT:       return "SET_EXCEPT";
    case IROpCode::PRINT_RESULTS:    return "PRINT_RESULTS";
    case IROpCode::LOAD_MATRIX:      return "LOAD_MATRIX";
    case IROpCode::SCAN_EXEC:        return "SCAN_EXEC";
    case IROpCode::SCAN_OPT_STRAND:  return "SCAN_OPT_STRAND";
    case IROpCode::SCAN_OPT_THRESHOLD:return "SCAN_OPT_THRESHOLD";
    case IROpCode::SCAN_ALIAS:       return "SCAN_ALIAS";
    case IROpCode::ANALYZE_GC:       return "ANALYZE_GC";
    case IROpCode::ANALYZE_CPG:      return "ANALYZE_CPG";
    case IROpCode::IF_BEGIN:         return "IF_BEGIN";
    case IROpCode::IF_ELSE:          return "IF_ELSE";
    case IROpCode::IF_END:           return "IF_END";
    default: return "UNKNOWN";
  }
}

IRGenerator::IRGenerator() : tempCounter(0) {}

std::string IRGenerator::newTemp(const std::string& prefix) {
  std::string sanitized = prefix;
  sanitized.erase(std::remove(sanitized.begin(), sanitized.end(), '"'), sanitized.end());
  return sanitized + "_" + std::to_string(tempCounter++);
}

const std::vector<IRInstruction>& IRGenerator::generate(ProgramNode* node) {
  instructions.clear();
  tempCounter = 0;
  if (node) node->accept(*this);
  return instructions;
}

void IRGenerator::visit(ProgramNode* node) {
  for (const auto& stmt : node->statements) {
    stmt->accept(*this);
  }
}

void IRGenerator::visit(LoadStmtNode* node) {
  IRInstruction instr;
  if (node->loadType == "SEQUENCE") {
    instr.opcode = IROpCode::LOAD_SEQ;
  } else if (node->loadType == "MATRIX") {
    instr.opcode = IROpCode::LOAD_MATRIX;
  } else {
    instr.opcode = IROpCode::LOAD_ANNOT;
  }
  instr.arg1 = node->filename;
  instr.arg2 = node->alias;
  instructions.push_back(instr);
}

void IRGenerator::visit(FindStmtNode* node) {
  currentTemp = newTemp("Find_" + node->motif);

  IRInstruction findInstr;
  findInstr.opcode = IROpCode::FIND_MOTIF;
  findInstr.arg1 = node->motif;
  findInstr.arg2 = currentTemp;
  instructions.push_back(findInstr);

  for (const auto& opt : node->opts) {
    opt->accept(*this);
  }

  IRInstruction execInstr;
  execInstr.opcode = IROpCode::FIND_EXEC;
  execInstr.arg1 = currentTemp;
  instructions.push_back(execInstr);

  if (!node->alias.empty()) {
    IRInstruction aliasInstr;
    aliasInstr.opcode = IROpCode::FIND_ALIAS;
    aliasInstr.arg1 = currentTemp;
    aliasInstr.arg2 = node->alias;
    instructions.push_back(aliasInstr);
  }

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = node->alias.empty() ? currentTemp : node->alias;
  printInstr.arg2 = "FIND";
  instructions.push_back(printInstr);
}

void IRGenerator::visit(FindOptNode* node) {
  IRInstruction instr;
  if (node->type == "WITHIN") {
    instr.opcode = IROpCode::FIND_OPT_WITHIN;
    instr.arg1 = node->value1;
    instr.arg2 = node->value2;
    instr.arg3 = node->value3;
    instr.arg4 = node->value4;
    instr.arg5 = node->value5;
  } else if (node->type == "STRAND") {
    instr.opcode = IROpCode::FIND_OPT_STRAND;
    instr.arg1 = node->value1;
  } else if (node->type == "CHR") {
    instr.opcode = IROpCode::FIND_OPT_CHR;
    instr.arg1 = node->value1;
  }
  instructions.push_back(instr);
}

void IRGenerator::visit(ExtractStmtNode* node) {
  currentTemp = newTemp("Extract_" + node->entity);

  IRInstruction extractInstr;
  extractInstr.opcode = IROpCode::EXTRACT;
  extractInstr.arg1 = node->entity;
  extractInstr.arg2 = currentTemp;
  instructions.push_back(extractInstr);

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = currentTemp;
  printInstr.arg2 = "EXTRACT";
  instructions.push_back(printInstr);
}

void IRGenerator::visit(SetOpStmtNode* node) {
  currentTemp = newTemp(node->op + "_" + node->entity1 + "_" + node->entity2);

  IRInstruction setInstr;
  if (node->op == "INTERSECT") setInstr.opcode = IROpCode::SET_INTERSECT;
  else if (node->op == "UNION") setInstr.opcode = IROpCode::SET_UNION;
  else setInstr.opcode = IROpCode::SET_EXCEPT;

  setInstr.arg1 = node->entity1;
  setInstr.arg2 = node->entity2;
  setInstr.arg3 = currentTemp;
  instructions.push_back(setInstr);

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = currentTemp;
  printInstr.arg2 = "SET_OP";
  instructions.push_back(printInstr);
}

void IRGenerator::visit(SimpleConditionNode* node) {
  IRInstruction filterInstr;
  if (node->property == "LENGTH") {
    filterInstr.opcode = IROpCode::FILTER_LENGTH;
  } else {
    filterInstr.opcode = IROpCode::FILTER_SIMILARITY;
  }
  filterInstr.arg1 = node->op;
  filterInstr.arg2 = node->value;
  filterInstr.arg3 = currentTemp;
  instructions.push_back(filterInstr);
}

void IRGenerator::visit(BinaryConditionNode* node) {
  if (node->left) node->left->accept(*this);
  if (node->right) node->right->accept(*this);
}

void IRGenerator::visit(NotConditionNode* node) {
  if (node->condition) node->condition->accept(*this);
}

void IRGenerator::visit(ScanStmtNode* node) {
  currentTemp = newTemp("Scan_" + node->matrixAlias);

  
  if (!node->strandFilter.empty()) {
    IRInstruction strandInstr;
    strandInstr.opcode = IROpCode::SCAN_OPT_STRAND;
    strandInstr.arg1 = node->strandFilter;
    instructions.push_back(strandInstr);
  }

  
  if (!node->threshold.empty()) {
    IRInstruction threshInstr;
    threshInstr.opcode = IROpCode::SCAN_OPT_THRESHOLD;
    threshInstr.arg1 = node->threshold;
    instructions.push_back(threshInstr);
  }

  
  IRInstruction scanInstr;
  scanInstr.opcode = IROpCode::SCAN_EXEC;
  scanInstr.arg1 = node->matrixAlias;
  scanInstr.arg2 = currentTemp;
  instructions.push_back(scanInstr);

  
  if (!node->alias.empty()) {
    IRInstruction aliasInstr;
    aliasInstr.opcode = IROpCode::SCAN_ALIAS;
    aliasInstr.arg1 = currentTemp;
    aliasInstr.arg2 = node->alias;
    instructions.push_back(aliasInstr);
  }

  
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  
  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = node->alias.empty() ? currentTemp : node->alias;
  printInstr.arg2 = "SCAN";
  instructions.push_back(printInstr);
}

void IRGenerator::visit(AnalyzeStmtNode *node) {
  currentTemp = newTemp("Analyze_" + node->analysisType);

  IRInstruction analyzeInstr;
  if (node->analysisType == "GC_CONTENT") {
    analyzeInstr.opcode = IROpCode::ANALYZE_GC;
  } else {
    analyzeInstr.opcode = IROpCode::ANALYZE_CPG;
  }
  
  analyzeInstr.arg1 = node->windowSize;
  analyzeInstr.arg2 = currentTemp;
  instructions.push_back(analyzeInstr);

  if (!node->alias.empty()) {
    IRInstruction aliasInstr;
    aliasInstr.opcode = IROpCode::SCAN_ALIAS; // We can reuse SCAN_ALIAS or FIND_ALIAS as they just assign names
    aliasInstr.arg1 = currentTemp;
    aliasInstr.arg2 = node->alias;
    instructions.push_back(aliasInstr);
  }

  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = node->alias.empty() ? currentTemp : node->alias;
  printInstr.arg2 = "ANALYZE";
  instructions.push_back(printInstr);
}

void IRGenerator::visit(IfStmtNode *node) {
  std::string prop, op, val;
  if (auto simple = dynamic_cast<SimpleConditionNode*>(node->condition.get())) {
    prop = simple->property;
    op = simple->op;
    val = simple->value;
  }
  IRInstruction ifBegin;
  ifBegin.opcode = IROpCode::IF_BEGIN;
  ifBegin.arg1 = prop;
  ifBegin.arg2 = op;
  ifBegin.arg3 = val;
  instructions.push_back(ifBegin);

  for (auto &stmt : node->thenStatements) {
    if (stmt) stmt->accept(*this);
  }

  if (!node->elseStatements.empty()) {
    IRInstruction ifElse;
    ifElse.opcode = IROpCode::IF_ELSE;
    instructions.push_back(ifElse);
    for (auto &stmt : node->elseStatements) {
      if (stmt) stmt->accept(*this);
    }
  }

  IRInstruction ifEnd;
  ifEnd.opcode = IROpCode::IF_END;
  instructions.push_back(ifEnd);
}

void IRGenerator::visit(ForeachStmtNode *node) {
  for (const auto &val : node->collection) {
    size_t startIdx = instructions.size();
    for (auto &stmt : node->bodyStatements) {
      if (stmt) stmt->accept(*this);
    }
    for (size_t i = startIdx; i < instructions.size(); ++i) {
      if (instructions[i].arg1 == node->iteratorVar) instructions[i].arg1 = val;
      if (instructions[i].arg2 == node->iteratorVar) instructions[i].arg2 = val;
      if (instructions[i].arg3 == node->iteratorVar) instructions[i].arg3 = val;
      if (instructions[i].arg4 == node->iteratorVar) instructions[i].arg4 = val;
      if (instructions[i].arg5 == node->iteratorVar) instructions[i].arg5 = val;
    }
  }
}
