#include "AST.h"
#include <iostream>

void SimpleConditionNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "SimpleCondition(" << property << " " << op << " " << value << ")" << std::endl;
}

void BinaryConditionNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "BinaryCondition(" << op << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (left) left->print(childPrefix, right == nullptr);
    if (right) right->print(childPrefix, true);
}

void NotConditionNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "NotCondition" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (condition) condition->print(childPrefix, true);
}

void LoadStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "LoadStmt(File: " << filename << ", Alias: " << alias << ")" << std::endl;
}

void FindOptNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "FindOpt(" << type;
    if (!value1.empty()) std::cout << ", " << value1;
    if (!value2.empty()) std::cout << " " << value2;
    if (!value3.empty()) std::cout << " " << value3;
    if (!value4.empty()) std::cout << " " << value4;
    if (!value5.empty()) std::cout << " " << value5;
    std::cout << ")" << std::endl;
}

void FindStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "FindStmt(Motif: " << motif << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    for (size_t i = 0; i < opts.size(); ++i) {
        opts[i]->print(childPrefix, i == opts.size() - 1);
    }
}

void ExtractStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "ExtractStmt(Entity: " << entity << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void SetOpStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "SetOperationStmt(" << op << " " << entity1 << " AND " << entity2 << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void ProgramNode::print(std::string prefix, bool isLast) const {
    if (prefix.empty()) {
        std::cout << "Program" << std::endl;
    } else {
        std::cout << prefix << (isLast ? "└── " : "├── ") << "Program" << std::endl;
    }
    std::string childPrefix = prefix.empty() ? "" : prefix + (isLast ? "    " : "│   ");
    for (size_t i = 0; i < statements.size(); ++i) {
        statements[i]->print(childPrefix, i == statements.size() - 1);
    }
}

void SimpleConditionNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void BinaryConditionNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void NotConditionNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void LoadStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void FindOptNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void FindStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ExtractStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void SetOpStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ProgramNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
