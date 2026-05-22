#include "AST.h"
#include <iostream>

void SimpleConditionNode::print(int indent) const {
    printIndent(indent);
    std::cout << "SimpleCondition(" << property << " " << op << " " << value << ")" << std::endl;
}

void BinaryConditionNode::print(int indent) const {
    printIndent(indent);
    std::cout << "BinaryCondition(" << op << ")" << std::endl;
    if (left) left->print(indent + 1);
    if (right) right->print(indent + 1);
}

void NotConditionNode::print(int indent) const {
    printIndent(indent);
    std::cout << "NotCondition" << std::endl;
    if (condition) condition->print(indent + 1);
}

void LoadStmtNode::print(int indent) const {
    printIndent(indent);
    std::cout << "LoadStmt(File: " << filename << ", Alias: " << alias << ")" << std::endl;
}

void FindOptNode::print(int indent) const {
    printIndent(indent);
    std::cout << "FindOpt(" << type;
    if (!value1.empty()) std::cout << ", " << value1;
    if (!value2.empty()) std::cout << " " << value2;
    if (!value3.empty()) std::cout << " " << value3;
    if (!value4.empty()) std::cout << " " << value4;
    if (!value5.empty()) std::cout << " " << value5;
    std::cout << ")" << std::endl;
}

void FindStmtNode::print(int indent) const {
    printIndent(indent);
    std::cout << "FindStmt(Motif: " << motif << ")" << std::endl;
    for (const auto& opt : opts) {
        opt->print(indent + 1);
    }
}

void ExtractStmtNode::print(int indent) const {
    printIndent(indent);
    std::cout << "ExtractStmt(Entity: " << entity << ")" << std::endl;
    if (whereClause) {
        printIndent(indent + 1);
        std::cout << "Where:" << std::endl;
        whereClause->print(indent + 2);
    }
}

void SetOpStmtNode::print(int indent) const {
    printIndent(indent);
    std::cout << "SetOperationStmt(" << op << " " << entity1 << " AND " << entity2 << ")" << std::endl;
    if (whereClause) {
        printIndent(indent + 1);
        std::cout << "Where:" << std::endl;
        whereClause->print(indent + 2);
    }
}

void ProgramNode::print(int indent) const {
    printIndent(indent);
    std::cout << "Program" << std::endl;
    for (const auto& stmt : statements) {
        stmt->print(indent + 1);
    }
}
