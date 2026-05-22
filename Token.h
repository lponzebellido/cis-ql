#ifndef TOKEN_H
#define TOKEN_H

#include <string>

enum class TokenType {
    FIND, SEARCH, EXTRACT, MATCH, COUNT, LOAD, SAVE,
    MOTIF, GENE, SEQUENCE, REGION, PROMOTER, ENHANCER,
    UPSTREAM, DOWNSTREAM, WITHIN, FROM, AT,
    BP, KB,
    AND, OR, NOT, WHERE, LENGTH, AS,
    EXON, INTRON, UTR, TSS,
    CHR, CHROMOSOME, STRAND, POSITIVE, NEGATIVE,
    INTERSECT, UNION, EXCEPT, OVERLAPS,
    SIMILARITY, REVERSE_COMPLEMENT,

    ID,
    NUM,       
    FLOAT,     
    STRING,

    // Symbols
    SEMICOLON, // ;
    COMMA,     // ,
    ASSIGN,    // =
    GREATER,   // >
    LESS,      // <
    GREATER_EQ,// >=
    LESS_EQ,   // <=
    LPAREN,    // (
    RPAREN,    // )
    PERCENT,   // %

    END_OF_FILE,
    UNKNOWN,
    ERROR_TOKEN
};

struct Token {
    TokenType type;
    std::string lexeme;
    int line;
    int column;
};

std::string tokenTypeToString(TokenType type);

#endif
