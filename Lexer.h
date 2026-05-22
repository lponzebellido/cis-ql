#ifndef LEXER_H
#define LEXER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "Token.h"
#include "SymbolTable.h"

#define MAXLENBUF 1024

class Lexer {
private:
    std::istream& in;
    char buffer[MAXLENBUF];
    int colIndex;
    int line;
    int bufferLen;
    bool eofReached;
    std::unordered_map<std::string, TokenType> keywords;
    SymbolTable* symbolTable;

    void initKeywords();
    char getChar();
    void ungetChar();
    void skipWhitespaceAndComments();
    
    Token createToken(TokenType type, const std::string& lexeme, int startLine, int startCol);
    Token identifierOrKeyword(char firstChar, int startLine, int startCol);
    Token number(char firstChar, int startLine, int startCol);
    Token stringLiteral(int startLine, int startCol);

public:
    Lexer(std::istream& inputStream, SymbolTable* st);
    std::vector<Token> tokenize();
};

#endif
