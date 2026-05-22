#include "Lexer.h"
#include <cctype>
#include <iostream>

Lexer::Lexer(std::istream& inputStream, SymbolTable* st) 
    : in(inputStream), colIndex(0), line(1), bufferLen(0), eofReached(false), symbolTable(st) {
    initKeywords();
}

void Lexer::initKeywords() {
    keywords["FIND"] = TokenType::FIND;
    keywords["SEARCH"] = TokenType::SEARCH;
    keywords["EXTRACT"] = TokenType::EXTRACT;
    keywords["MATCH"] = TokenType::MATCH;
    keywords["COUNT"] = TokenType::COUNT;
    keywords["LOAD"] = TokenType::LOAD;
    keywords["SAVE"] = TokenType::SAVE;
    keywords["MOTIF"] = TokenType::MOTIF;
    keywords["GENE"] = TokenType::GENE;
    keywords["SEQUENCE"] = TokenType::SEQUENCE;
    keywords["REGION"] = TokenType::REGION;
    keywords["PROMOTER"] = TokenType::PROMOTER;
    keywords["ENHANCER"] = TokenType::ENHANCER;
    keywords["UPSTREAM"] = TokenType::UPSTREAM;
    keywords["DOWNSTREAM"] = TokenType::DOWNSTREAM;
    keywords["WITHIN"] = TokenType::WITHIN;
    keywords["FROM"] = TokenType::FROM;
    keywords["AT"] = TokenType::AT;
    keywords["BP"] = TokenType::BP;
    keywords["KB"] = TokenType::KB;
    keywords["MB"] = TokenType::MB;
    keywords["AND"] = TokenType::AND;
    keywords["OR"] = TokenType::OR;
    keywords["NOT"] = TokenType::NOT;
    keywords["WHERE"] = TokenType::WHERE;
    keywords["LENGTH"] = TokenType::LENGTH;
    keywords["AS"] = TokenType::AS;

    keywords["EXON"] = TokenType::EXON;
    keywords["INTRON"] = TokenType::INTRON;
    keywords["UTR"] = TokenType::UTR;
    keywords["TSS"] = TokenType::TSS;

    keywords["CHR"] = TokenType::CHR;
    keywords["CHROMOSOME"] = TokenType::CHROMOSOME;
    keywords["STRAND"] = TokenType::STRAND;
    keywords["POSITIVE"] = TokenType::POSITIVE;
    keywords["NEGATIVE"] = TokenType::NEGATIVE;

    keywords["INTERSECT"] = TokenType::INTERSECT;
    keywords["UNION"] = TokenType::UNION;
    keywords["EXCEPT"] = TokenType::EXCEPT;
    keywords["OVERLAPS"] = TokenType::OVERLAPS;

    keywords["SIMILARITY"] = TokenType::SIMILARITY;
    keywords["REVERSE_COMPLEMENT"] = TokenType::REVERSE_COMPLEMENT;
}

char Lexer::getChar() {
    if (eofReached) return EOF;
    
    if (colIndex >= bufferLen) {
        if (!in) {
            eofReached = true;
            return EOF;
        }
        in.getline(buffer, MAXLENBUF);
        if (in.fail() && !in.eof()) {
            in.clear(); 
        }
        bufferLen = std::char_traits<char>::length(buffer);
        if (!in.eof() || bufferLen > 0) {
            if (!in.eof()) {
                buffer[bufferLen++] = '\n';
                buffer[bufferLen] = '\0';
            }
        }
        colIndex = 0;
        if (bufferLen == 0) {
            eofReached = true;
            return EOF;
        }
    }
    
    char c = buffer[colIndex++];
    if (c == '\n') {
        line++;
    }
    return c;
}

void Lexer::ungetChar() {
    if (colIndex > 0) {
        colIndex--;
        if (buffer[colIndex] == '\n') {
            line--;
        }
    }
}

void Lexer::skipWhitespaceAndComments() {
    while (true) {
        char c = getChar();
        if (c == EOF) break;
        
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            continue;
        } else if (c == '/') {
            char next = getChar();
            if (next == '/') {
                while ((c = getChar()) != '\n' && c != EOF);
            } else if (next == '*') {
                while (true) {
                    c = getChar();
                    if (c == EOF) break;
                    if (c == '*') {
                        if (getChar() == '/') break;
                        ungetChar();
                    }
                }
            } else {
                ungetChar(); 
                ungetChar();
            }
        } else {
            ungetChar();
            break;
        }
    }
}

Token Lexer::createToken(TokenType type, const std::string& lexeme, int startLine, int startCol) {
    return {type, lexeme, startLine, startCol};
}

Token Lexer::identifierOrKeyword(char firstChar, int startLine, int startCol) {
    std::string lexeme = "";
    lexeme += firstChar;
    while (true) {
        char c = getChar();
        if (c == EOF) break;
        if (std::isalnum(c) || c == '_') {
            lexeme += c;
        } else {
            ungetChar();
            break;
        }
    }

    if (keywords.find(lexeme) != keywords.end()) {
        return createToken(keywords[lexeme], lexeme, startLine, startCol);
    }
    
    symbolTable->insert(lexeme, "ID");
    return createToken(TokenType::ID, lexeme, startLine, startCol);
}

Token Lexer::number(char firstChar, int startLine, int startCol) {
    std::string lexeme = "";
    lexeme += firstChar;
    bool isFloat = false;
    
    while (true) {
        char c = getChar();
        if (c == EOF) break;
        if (std::isdigit(c)) {
            lexeme += c;
        } else if (c == '.' && !isFloat) {
            char next = getChar();
            if (std::isdigit(next)) {
                isFloat = true;
                lexeme += c;
                lexeme += next;
            } else {
                if (next != EOF) ungetChar();
                ungetChar(); 
                break;
            }
        } else {
            ungetChar();
            break;
        }
    }
    return createToken(isFloat ? TokenType::FLOAT : TokenType::NUM, lexeme, startLine, startCol);
}

Token Lexer::stringLiteral(int startLine, int startCol) {
    std::string lexeme = "";
    while (true) {
        char c = getChar();
        if (c == EOF || c == '"') {
            break;
        }
        lexeme += c;
    }
    symbolTable->insert(lexeme, "LITERAL_STRING");
    return createToken(TokenType::STRING, "\"" + lexeme + "\"", startLine, startCol);
}

std::vector<Token> Lexer::tokenize() {
    std::vector<Token> tokens;
    
    while (true) {
        skipWhitespaceAndComments();
        
        int startLine = line;
        char c = getChar();
        
        if (c == EOF) {
            tokens.push_back(createToken(TokenType::END_OF_FILE, "EOF", line, colIndex));
            break;
        }
        
        int startCol = colIndex;
        
        if (std::isalpha(c) || c == '_') {
            tokens.push_back(identifierOrKeyword(c, startLine, startCol));
        } else if (std::isdigit(c)) {
            tokens.push_back(number(c, startLine, startCol));
        } else if (c == '"') {
            tokens.push_back(stringLiteral(startLine, startCol));
        } else {
            std::string lexeme(1, c);
            switch (c) {
                case ';': tokens.push_back(createToken(TokenType::SEMICOLON, lexeme, startLine, startCol)); break;
                case ',': tokens.push_back(createToken(TokenType::COMMA, lexeme, startLine, startCol)); break;
                case '(': tokens.push_back(createToken(TokenType::LPAREN, lexeme, startLine, startCol)); break;
                case ')': tokens.push_back(createToken(TokenType::RPAREN, lexeme, startLine, startCol)); break;
                case '%': tokens.push_back(createToken(TokenType::PERCENT, lexeme, startLine, startCol)); break;
                case '=': 
                    tokens.push_back(createToken(TokenType::ASSIGN, lexeme, startLine, startCol)); 
                    break;
                case '>': {
                    char next = getChar();
                    if (next == '=') {
                        lexeme += next;
                        tokens.push_back(createToken(TokenType::GREATER_EQ, lexeme, startLine, startCol));
                    } else {
                        if (next != EOF) ungetChar();
                        tokens.push_back(createToken(TokenType::GREATER, lexeme, startLine, startCol));
                    }
                    break;
                }
                case '<': {
                    char next = getChar();
                    if (next == '=') {
                        lexeme += next;
                        tokens.push_back(createToken(TokenType::LESS_EQ, lexeme, startLine, startCol));
                    } else {
                        if (next != EOF) ungetChar();
                        tokens.push_back(createToken(TokenType::LESS, lexeme, startLine, startCol));
                    }
                    break;
                }
                default:
                    std::cerr << "Lexical Error at L" << startLine << ":C" << startCol 
                              << " - Invalid character '" << c << "'" << std::endl;
                    tokens.push_back(createToken(TokenType::ERROR_TOKEN, lexeme, startLine, startCol));
                    break;
            }
        }
    }
    
    return tokens;
}

std::string tokenTypeToString(TokenType type) {
    switch (type) {
        case TokenType::FIND: return "FIND";
        case TokenType::SEARCH: return "SEARCH";
        case TokenType::EXTRACT: return "EXTRACT";
        case TokenType::MATCH: return "MATCH";
        case TokenType::COUNT: return "COUNT";
        case TokenType::LOAD: return "LOAD";
        case TokenType::SAVE: return "SAVE";
        case TokenType::MOTIF: return "MOTIF";
        case TokenType::GENE: return "GENE";
        case TokenType::SEQUENCE: return "SEQUENCE";
        case TokenType::REGION: return "REGION";
        case TokenType::PROMOTER: return "PROMOTER";
        case TokenType::ENHANCER: return "ENHANCER";
        case TokenType::UPSTREAM: return "UPSTREAM";
        case TokenType::DOWNSTREAM: return "DOWNSTREAM";
        case TokenType::WITHIN: return "WITHIN";
        case TokenType::FROM: return "FROM";
        case TokenType::AT: return "AT";
        case TokenType::BP: return "BP";
        case TokenType::KB: return "KB";
        case TokenType::MB: return "MB";
        case TokenType::AND: return "AND";
        case TokenType::OR: return "OR";
        case TokenType::NOT: return "NOT";
        case TokenType::WHERE: return "WHERE";
        case TokenType::LENGTH: return "LENGTH";
        case TokenType::AS: return "AS";
        case TokenType::EXON: return "EXON";
        case TokenType::INTRON: return "INTRON";
        case TokenType::UTR: return "UTR";
        case TokenType::TSS: return "TSS";
        case TokenType::CHR: return "CHR";
        case TokenType::CHROMOSOME: return "CHROMOSOME";
        case TokenType::STRAND: return "STRAND";
        case TokenType::POSITIVE: return "POSITIVE";
        case TokenType::NEGATIVE: return "NEGATIVE";
        case TokenType::INTERSECT: return "INTERSECT";
        case TokenType::UNION: return "UNION";
        case TokenType::EXCEPT: return "EXCEPT";
        case TokenType::OVERLAPS: return "OVERLAPS";
        case TokenType::SIMILARITY: return "SIMILARITY";
        case TokenType::REVERSE_COMPLEMENT: return "REVERSE_COMPLEMENT";
        case TokenType::ID: return "ID";
        case TokenType::NUM: return "NUM";
        case TokenType::FLOAT: return "FLOAT";
        case TokenType::STRING: return "STRING";
        case TokenType::SEMICOLON: return "SEMICOLON";
        case TokenType::COMMA: return "COMMA";
        case TokenType::ASSIGN: return "ASSIGN";
        case TokenType::GREATER: return "GREATER";
        case TokenType::LESS: return "LESS";
        case TokenType::GREATER_EQ: return "GREATER_EQ";
        case TokenType::LESS_EQ: return "LESS_EQ";
        case TokenType::LPAREN: return "LPAREN";
        case TokenType::RPAREN: return "RPAREN";
        case TokenType::PERCENT: return "PERCENT";
        case TokenType::END_OF_FILE: return "EOF";
        case TokenType::ERROR_TOKEN: return "ERROR_TOKEN";
        default: return "UNKNOWN";
    }
}
