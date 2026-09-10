#ifndef TOKEN_H
#define TOKEN_H

#include <string>

enum class TokenType {
    FIND, SEARCH, EXTRACT, MATCH, COUNT, LOAD, SAVE, USE, EXPORT, DEFINE,
    MOTIF, GENE, SEQUENCE, ANNOTATION, REGION, PROMOTER, PROMOTERS, ENHANCER,
    OF,
    UPSTREAM, DOWNSTREAM, WITHIN, FROM, AT, TO, FORMAT,
    BED, GFF3, TSV,
    BP, KB, MB,
    AND, OR, NOT, WHERE, LENGTH, AS,
    EXON, INTRON, UTR, TSS, CDS,
    CHR, CHROMOSOME, STRAND, POSITIVE, NEGATIVE,
    INTERSECT, UNION, EXCEPT, OVERLAPS,
    SIMILARITY, REVERSE_COMPLEMENT,
    MATRIX, SCAN, THRESHOLD,
    ANALYZE, GC_CONTENT, CPG_ISLANDS, WINDOW,
    IF, THEN, ELSE, ENDIF,
    FOREACH, IN, DO, ENDFOR,

    ID,
    NUM,       
    FLOAT,     
    STRING,

    SEMICOLON, 
    COMMA,     
    ASSIGN,    
    GREATER,   
    LESS,      
    GREATER_EQ,
    LESS_EQ,   
    LPAREN,    
    RPAREN,    
    LBRACKET,
    RBRACKET,
    PERCENT,   

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
