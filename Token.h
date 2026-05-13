#ifndef TOKEN_H
#define TOKEN_H

#include <string>

enum class TokenType {
    FIND, SEARCH, EXTRACT, MATCH, COUNT, LOAD, SAVE,
    MOTIF, GENE, SEQUENCE, REGION, PROMOTER, ENHANCER,
    UPSTREAM, DOWNSTREAM, WITHIN, FROM, AT,
    BP, KB,
    AND, OR, NOT, WHERE, LENGTH, AS,

    ID,
    NUM,
    STRING,

    SEMICOLON, // ;
    COMMA,     // ,
    ASSIGN,    // =
    GREATER,   // >
    LESS,      // <
    GREATER_EQ,// >=
    LESS_EQ,   // <=
    LPAREN,    // (
    RPAREN,    // )

    END_OF_FILE,
    UNKNOWN
};

struct Token {
    TokenType type;
    std::string lexeme;
    int line;
    int column;
};

std::string tokenTypeToString(TokenType type);

#endif
