#ifndef LEXER_H
#define LEXER_H

#include <string>
#include <vector>
#include <unordered_map>
#include "Token.h"
#include "SymbolTable.h"

class Lexer {
private:
    std::string source;
    size_t pos;
    int line;
    int column;
    std::unordered_map<std::string, TokenType> keywords;
    SymbolTable* symbolTable;

    void initKeywords();
    char advance();
    char peek();
    char peekNext();
    void skipWhitespace();
    
    Token createToken(TokenType type, const std::string& lexeme);
    Token identifierOrKeyword();
    Token number();
    Token stringLiteral();

public:
    Lexer(const std::string& sourceCode, SymbolTable* st);
    std::vector<Token> tokenize();
};

#endif
