#include "IRGenerator.h"

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
    default: return "UNKNOWN";
  }
}

IRGenerator::IRGenerator() : tempCounter(0) {}

std::string IRGenerator::newTemp() {
  return "t" + std::to_string(tempCounter++);
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
  currentTemp = newTemp();

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
  printInstr.arg1 = currentTemp;
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
  currentTemp = newTemp();

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
  currentTemp = newTemp();

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
  currentTemp = newTemp();

  // Strand option
  if (!node->strandFilter.empty()) {
    IRInstruction strandInstr;
    strandInstr.opcode = IROpCode::SCAN_OPT_STRAND;
    strandInstr.arg1 = node->strandFilter;
    instructions.push_back(strandInstr);
  }

  // Threshold option
  if (!node->threshold.empty()) {
    IRInstruction threshInstr;
    threshInstr.opcode = IROpCode::SCAN_OPT_THRESHOLD;
    threshInstr.arg1 = node->threshold;
    instructions.push_back(threshInstr);
  }

  // Execute scan
  IRInstruction scanInstr;
  scanInstr.opcode = IROpCode::SCAN_EXEC;
  scanInstr.arg1 = node->matrixAlias;
  scanInstr.arg2 = currentTemp;
  instructions.push_back(scanInstr);

  // Alias
  if (!node->alias.empty()) {
    IRInstruction aliasInstr;
    aliasInstr.opcode = IROpCode::SCAN_ALIAS;
    aliasInstr.arg1 = currentTemp;
    aliasInstr.arg2 = node->alias;
    instructions.push_back(aliasInstr);
  }

  // Where clause
  if (node->whereClause) {
    node->whereClause->accept(*this);
  }

  // Print results
  IRInstruction printInstr;
  printInstr.opcode = IROpCode::PRINT_RESULTS;
  printInstr.arg1 = currentTemp;
  printInstr.arg2 = "SCAN";
  instructions.push_back(printInstr);
}
