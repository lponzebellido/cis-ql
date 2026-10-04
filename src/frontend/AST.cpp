#include "AST.h"
#include <iostream>

void SimpleConditionNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "SimpleCondition(";
    if (!sliceStart.empty()) std::cout << "SLICE ";
    std::cout << property;
    if (!reference.empty()) std::cout << (property == "ATTRIBUTE" ? " " : " TO ") << reference;
    if (!modifier.empty()) std::cout << " " << modifier << " " << modifierValue;
    if (!sliceStart.empty()) std::cout << " FROM " << sliceStart << " TO " << sliceEnd;
    std::cout << " " << op << " " << value << ")" << std::endl;
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
    std::cout << prefix << (isLast ? "└── " : "├── ") << "LoadStmt(" << loadType << " File: " << filename;
    if (!format.empty()) std::cout << ", Format: " << format;
    if (!evidenceClass.empty()) std::cout << ", Evidence: " << evidenceClass;
    if (!assay.empty()) std::cout << ", Assay: " << assay;
    if (!sample.empty()) std::cout << ", Sample: " << sample;
    if (!condition.empty()) std::cout << ", Condition: " << condition;
    if (!replicate.empty()) std::cout << ", Replicate: " << replicate;
    if (!control.empty()) std::cout << ", Control: " << control;
    std::cout << ", Alias: " << alias << ")" << std::endl;
}

void UseStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "UseStmt(" << datasetType << ": " << alias << ")"
              << std::endl;
}

void ExportStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "ExportStmt(" << alias << " TO " << filename
              << " FORMAT " << format << ")" << std::endl;
}

void DefinePromotersStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "DefinePromotersStmt(Source: " << source
              << ", From: TSS, Upstream: " << upstreamValue << " "
              << upstreamUnit << ", Downstream: " << downstreamValue << " "
              << downstreamUnit << ", AS: " << alias << ")" << std::endl;
}

void DefineModuleStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "DefineModuleStmt(From: " << firstSet << " WITH "
              << secondSet << ", Spacing: " << minimumSpacingValue << " "
              << minimumSpacingUnit << " TO " << maximumSpacingValue << " "
              << maximumSpacingUnit << ", Order: " << orderPolicy
              << ", Orientation: " << orientationPolicy << ", AS: "
              << alias << ")" << std::endl;
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
    std::cout << prefix << (isLast ? "└── " : "├── ") << "FindStmt(Motif: " << motif;
    if (!alias.empty()) std::cout << ", AS: " << alias;
    std::cout << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    bool hasWhere = whereClause != nullptr;
    for (size_t i = 0; i < opts.size(); ++i) {
        opts[i]->print(childPrefix, !hasWhere && i == opts.size() - 1);
    }
    if (hasWhere) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void ExtractStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "ExtractStmt(";
    if (relation.empty())
        std::cout << "Entity: " << entity;
    else {
        std::cout << "Relation: " << relation << ", Source: " << source;
        if (relation == "TRANSCRIPTS") {
            std::cout << ", Types: [";
            for (size_t index = 0; index < transcriptTypes.size(); ++index) {
                if (index > 0) std::cout << ", ";
                std::cout << transcriptTypes[index];
            }
            std::cout << "], Policy: " << transcriptPolicy;
            if (!transcriptPolicyKey.empty())
                std::cout << " " << transcriptPolicyKey;
            if (!transcriptPolicyValue.empty())
                std::cout << " " << transcriptPolicyValue;
        }
    }
    if (!alias.empty()) std::cout << ", AS: " << alias;
    std::cout << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void SetOpStmtNode::print(std::string prefix, bool isLast) const {
    if (op == "CONSENSUS") {
        std::cout << prefix << (isLast ? "└── " : "├── ")
                  << "ConsensusStmt(From: [";
        for (size_t index = 0; index < entities.size(); ++index) {
            if (index > 0) std::cout << ", ";
            std::cout << entities[index];
        }
        std::cout << "], Anchor: " << anchor
                  << ", Minimum support: " << minimumSupport;
        if (!minimumReciprocalOverlap.empty())
            std::cout << ", Minimum reciprocal overlap: "
                      << minimumReciprocalOverlap << "%";
        if (!maximumSummitDistanceValue.empty())
            std::cout << ", Maximum summit distance: "
                      << maximumSummitDistanceValue << " "
                      << maximumSummitDistanceUnit;
        std::cout << ", AS: " << alias << ")" << std::endl;
        if (whereClause) {
            std::string childPrefix = prefix + (isLast ? "    " : "│   ");
            std::cout << childPrefix << "└── Where:" << std::endl;
            whereClause->print(childPrefix + "    ", true);
        }
        return;
    }
    const std::string separator = op == "EXCEPT" ? " FROM "
                                  : op == "OVERLAPS" ? " WITH "
                                  : op == "NEAR" ? " TO " : " AND ";
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "SetOperationStmt(" << op << " " << entity1
              << separator << entity2;
    if (op == "NEAR")
        std::cout << " WITHIN " << distanceValue << " " << distanceUnit;
    if (!alias.empty()) std::cout << ", AS: " << alias;
    std::cout << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void ValidateStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "ValidateStmt(Annotation: " << annotationAlias
              << ", AS: " << alias << ")" << std::endl;
}

void GroupStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "GroupStmt(Source: " << sourceAlias
              << ", By: ID, AS: " << alias << ")" << std::endl;
}

void TranslateStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "TranslateStmt(Source: " << sourceAlias
              << ", Code: " << geneticCode << ", Frame: " << frame
              << ", AS: " << alias << ")" << std::endl;
}

void CountStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ")
              << "CountStmt(" << countedEntity << " IN " << containerEntity
              << ", AS: " << alias << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void ScanStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "ScanStmt(Matrix: " << matrixAlias;
    if (!target.empty()) std::cout << ", IN: " << target;
    if (!strandFilter.empty()) std::cout << ", Strand: " << strandFilter;
    if (!threshold.empty()) std::cout << ", Threshold: " << threshold;
    if (!significanceMetric.empty()) {
        std::cout << ", " << significanceMetric << " "
                  << significanceOperator << " " << significanceThreshold;
    }
    if (!backgroundMode.empty()) {
        std::cout << ", Background: " << backgroundMode;
        if (!backgroundSource.empty()) std::cout << " " << backgroundSource;
    }
    if (!alias.empty()) std::cout << ", AS: " << alias;
    std::cout << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void AnalyzeStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "AnalyzeStmt(" << analysisType;
    if (!windowSize.empty()) std::cout << ", Window: " << windowSize;
    if (!alias.empty()) std::cout << ", AS: " << alias;
    std::cout << ")" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (whereClause) {
        std::cout << childPrefix << "└── Where:" << std::endl;
        whereClause->print(childPrefix + "    ", true);
    }
}

void IfStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "IfStmt" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    if (condition) condition->print(childPrefix, false);
    for (size_t i = 0; i < thenStatements.size(); ++i) {
        thenStatements[i]->print(childPrefix, elseStatements.empty() && i == thenStatements.size() - 1);
    }
    for (size_t i = 0; i < elseStatements.size(); ++i) {
        elseStatements[i]->print(childPrefix, i == elseStatements.size() - 1);
    }
}

void ForeachStmtNode::print(std::string prefix, bool isLast) const {
    std::cout << prefix << (isLast ? "└── " : "├── ") << "ForeachStmt(" << iteratorVar << " IN [";
    for (size_t i = 0; i < collection.size(); ++i) {
        std::cout << collection[i] << (i + 1 < collection.size() ? ", " : "");
    }
    std::cout << "])" << std::endl;
    std::string childPrefix = prefix + (isLast ? "    " : "│   ");
    for (size_t i = 0; i < bodyStatements.size(); ++i) {
        bodyStatements[i]->print(childPrefix, i == bodyStatements.size() - 1);
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
void UseStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ExportStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void DefinePromotersStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void DefineModuleStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void FindOptNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void FindStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ExtractStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void SetOpStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void CountStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ScanStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void AnalyzeStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ValidateStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void GroupStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void TranslateStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void IfStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ForeachStmtNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
void ProgramNode::accept(ASTVisitor& visitor) { visitor.visit(this); }
