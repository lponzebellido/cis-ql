#include "Parser.h"
#include <iostream>
#include <stdexcept>

Parser::Parser(const std::vector<Token>& tokenList) : tokens(tokenList), current(0), hasError(false) {}

const Token& Parser::peek() const {
    return tokens[current];
}

const Token& Parser::previous() const {
    return tokens[current - 1];
}

bool Parser::isAtEnd() const {
    return peek().type == TokenType::END_OF_FILE;
}

Token Parser::advance() {
    if (!isAtEnd()) current++;
    return previous();
}

bool Parser::check(TokenType type) const {
    if (isAtEnd()) return false;
    return peek().type == type;
}

bool Parser::match(TokenType type) {
    if (check(type)) {
        advance();
        return true;
    }
    return false;
}

void Parser::consume(TokenType type, const std::string& message) {
    if (check(type)) {
        advance();
        return;
    }
    reportError(peek(), message);
    throw std::runtime_error("Parse error");
}

void Parser::reportError(const Token& token, const std::string& message) {
    hasError = true;
    std::cerr << "Error Sintactico en L" << token.line << ":C" << token.column << " a las '" << token.lexeme << "': " << message << std::endl;
}

void Parser::synchronize() {
    advance();
    while (!isAtEnd()) {
        if (previous().type == TokenType::SEMICOLON) return;
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
        } catch (std::runtime_error&) {
            synchronize();
        }
    }
    return program;
}

std::unique_ptr<StatementNode> Parser::parseStatement() {
    if (match(TokenType::LOAD)) return parseLoad();
    if (match(TokenType::FIND)) return parseFind();
    if (match(TokenType::EXTRACT)) return parseExtract();
    if (match(TokenType::INTERSECT) || match(TokenType::UNION) || match(TokenType::EXCEPT)) {
        current--; // revert match
        return parseSetOperation();
    }
    
    // Ignore ERROR_TOKEN and keep parsing
    if (match(TokenType::ERROR_TOKEN)) {
        return nullptr;
    }

    reportError(peek(), "Se esperaba el inicio de una sentencia (LOAD, FIND, EXTRACT, INTERSECT, UNION, EXCEPT)");
    throw std::runtime_error("Parse error");
}

std::unique_ptr<LoadStmtNode> Parser::parseLoad() {
    consume(TokenType::SEQUENCE, "Se esperaba 'SEQUENCE' despues de LOAD.");
    consume(TokenType::STRING, "Se esperaba un nombre de archivo (cadena).");
    std::string file = previous().lexeme;
    consume(TokenType::AS, "Se esperaba 'AS' despues del nombre del archivo.");
    consume(TokenType::ID, "Se esperaba un identificador alias.");
    std::string alias = previous().lexeme;
    consume(TokenType::SEMICOLON, "Se esperaba ';' al final de la instruccion LOAD.");
    return std::unique_ptr<LoadStmtNode>(new LoadStmtNode(file, alias));
}

std::unique_ptr<FindStmtNode> Parser::parseFind() {
    consume(TokenType::MOTIF, "Se esperaba 'MOTIF' despues de FIND.");
    consume(TokenType::STRING, "Se esperaba el nombre del motivo (cadena).");
    auto node = std::unique_ptr<FindStmtNode>(new FindStmtNode(previous().lexeme));
    
    while (!check(TokenType::SEMICOLON) && !isAtEnd()) {
        auto opt = std::unique_ptr<FindOptNode>(new FindOptNode());
        if (match(TokenType::WITHIN)) {
            opt->type = "WITHIN";
            consume(TokenType::NUM, "Se esperaba un numero para WITHIN.");
            opt->value1 = previous().lexeme;
            
            if (match(TokenType::BP) || match(TokenType::KB) || match(TokenType::MB)) {
                opt->value2 = previous().lexeme;
            }
            
            if (match(TokenType::UPSTREAM) || match(TokenType::DOWNSTREAM)) {
                opt->value3 = previous().lexeme;
            } else {
                reportError(peek(), "Se esperaba UPSTREAM o DOWNSTREAM.");
                throw std::runtime_error("Parse error");
            }
            
            consume(TokenType::FROM, "Se esperaba 'FROM' despues de la direccion.");
            
            if (match(TokenType::GENE) || match(TokenType::PROMOTER) || match(TokenType::ENHANCER) || 
                match(TokenType::EXON) || match(TokenType::INTRON) || match(TokenType::UTR) || match(TokenType::TSS) || match(TokenType::REGION)) {
                opt->value4 = previous().lexeme;
            } else {
                reportError(peek(), "Se esperaba una entidad biologica.");
                throw std::runtime_error("Parse error");
            }
            
            consume(TokenType::STRING, "Se esperaba el nombre de la entidad.");
            opt->value5 = previous().lexeme;
            
        } else if (match(TokenType::STRAND)) {
            opt->type = "STRAND";
            if (match(TokenType::POSITIVE) || match(TokenType::NEGATIVE)) {
                opt->value1 = previous().lexeme;
            } else {
                reportError(peek(), "Se esperaba POSITIVE o NEGATIVE.");
                throw std::runtime_error("Parse error");
            }
        } else if (match(TokenType::CHR)) {
            opt->type = "CHR";
            consume(TokenType::STRING, "Se esperaba el nombre del cromosoma.");
            opt->value1 = previous().lexeme;
        } else {
            reportError(peek(), "Opcion no reconocida en FIND.");
            throw std::runtime_error("Parse error");
        }
        node->opts.push_back(std::move(opt));
    }
    consume(TokenType::SEMICOLON, "Se esperaba ';' al final de FIND.");
    return node;
}

std::unique_ptr<ExtractStmtNode> Parser::parseExtract() {
    if (match(TokenType::GENE) || match(TokenType::PROMOTER) || match(TokenType::ENHANCER) || 
        match(TokenType::EXON) || match(TokenType::INTRON) || match(TokenType::UTR) || match(TokenType::TSS) || match(TokenType::REGION)) {
        std::string entity = previous().lexeme;
        auto whereClause = parseWhereClause();
        consume(TokenType::SEMICOLON, "Se esperaba ';' al final de EXTRACT.");
        return std::unique_ptr<ExtractStmtNode>(new ExtractStmtNode(entity, std::move(whereClause)));
    }
    reportError(peek(), "Se esperaba una entidad para EXTRACT.");
    throw std::runtime_error("Parse error");
}

std::unique_ptr<SetOpStmtNode> Parser::parseSetOperation() {
    std::string op;
    if (match(TokenType::INTERSECT)) op = "INTERSECT";
    else if (match(TokenType::UNION)) op = "UNION";
    else if (match(TokenType::EXCEPT)) op = "EXCEPT";
    
    if (match(TokenType::GENE) || match(TokenType::PROMOTER) || match(TokenType::ENHANCER) || 
        match(TokenType::EXON) || match(TokenType::INTRON) || match(TokenType::UTR) || match(TokenType::TSS) || match(TokenType::REGION)) {
        std::string e1 = previous().lexeme;
        
        std::string sepError = "Se esperaba 'AND' (para INTERSECT/UNION) o 'FROM' (para EXCEPT).";
        if (op == "EXCEPT") consume(TokenType::FROM, sepError);
        else consume(TokenType::AND, sepError);
        
        if (match(TokenType::GENE) || match(TokenType::PROMOTER) || match(TokenType::ENHANCER) || 
            match(TokenType::EXON) || match(TokenType::INTRON) || match(TokenType::UTR) || match(TokenType::TSS) || match(TokenType::REGION)) {
            std::string e2 = previous().lexeme;
            auto whereClause = parseWhereClause();
            consume(TokenType::SEMICOLON, "Se esperaba ';' al final de la operacion de conjuntos.");
            return std::unique_ptr<SetOpStmtNode>(new SetOpStmtNode(op, e1, e2, std::move(whereClause)));
        }
    }
    reportError(peek(), "Se esperaba una entidad en la operacion de conjuntos.");
    throw std::runtime_error("Parse error");
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

std::unique_ptr<ConditionNode> Parser::parseConditionPrime(std::unique_ptr<ConditionNode> left) {
    if (match(TokenType::OR)) {
        auto rightTerm = parseTerm();
        auto newLeft = std::unique_ptr<BinaryConditionNode>(new BinaryConditionNode("OR", std::move(left), std::move(rightTerm)));
        return parseConditionPrime(std::move(newLeft));
    }
    return left;
}

std::unique_ptr<ConditionNode> Parser::parseTerm() {
    auto factor = parseFactor();
    return parseTermPrime(std::move(factor));
}

std::unique_ptr<ConditionNode> Parser::parseTermPrime(std::unique_ptr<ConditionNode> left) {
    if (match(TokenType::AND)) {
        auto rightFactor = parseFactor();
        auto newLeft = std::unique_ptr<BinaryConditionNode>(new BinaryConditionNode("AND", std::move(left), std::move(rightFactor)));
        return parseTermPrime(std::move(newLeft));
    }
    return left;
}

std::unique_ptr<ConditionNode> Parser::parseFactor() {
    if (match(TokenType::NOT)) {
        return std::unique_ptr<NotConditionNode>(new NotConditionNode(parseFactor()));
    }
    if (match(TokenType::LPAREN)) {
        auto cond = parseCondition();
        consume(TokenType::RPAREN, "Se esperaba ')' despues de la condicion.");
        return cond;
    }
    return parseSimpleCondition();
}

std::unique_ptr<SimpleConditionNode> Parser::parseSimpleCondition() {
    std::string prop;
    if (match(TokenType::LENGTH) || match(TokenType::SIMILARITY)) {
        prop = previous().lexeme;
    } else {
        reportError(peek(), "Se esperaba 'LENGTH' o 'SIMILARITY'.");
        throw std::runtime_error("Parse error");
    }
    
    std::string op;
    if (match(TokenType::GREATER) || match(TokenType::LESS) || match(TokenType::GREATER_EQ) || match(TokenType::LESS_EQ) || match(TokenType::ASSIGN)) {
        op = previous().lexeme;
    } else {
        reportError(peek(), "Se esperaba un operador relacional (> < >= <= =).");
        throw std::runtime_error("Parse error");
    }
    
    std::string val;
    if (match(TokenType::NUM) || match(TokenType::FLOAT) || match(TokenType::STRING)) {
        val = previous().lexeme;
        if (match(TokenType::BP) || match(TokenType::KB) || match(TokenType::MB) || match(TokenType::PERCENT)) {
            val += " " + previous().lexeme;
        }
    } else {
        reportError(peek(), "Se esperaba un valor (NUM, FLOAT, STRING).");
        throw std::runtime_error("Parse error");
    }
    
    return std::unique_ptr<SimpleConditionNode>(new SimpleConditionNode(prop, op, val));
}
