#include "Lexer.h"
#include <cctype>
#include <iostream>

Lexer::Lexer(std::istream& inputStream)
    : in(inputStream), colIndex(0), line(1), bufferLen(0), eofReached(false),
      hasError(false) {
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
    keywords["USE"] = TokenType::USE;
    keywords["EXPORT"] = TokenType::EXPORT;
    keywords["DEFINE"] = TokenType::DEFINE;
    keywords["MOTIF"] = TokenType::MOTIF;
    keywords["MODULE"] = TokenType::MODULE;
    keywords["GENE"] = TokenType::GENE;
    keywords["SEQUENCE"] = TokenType::SEQUENCE;
    keywords["ANNOTATION"] = TokenType::ANNOTATION;
    keywords["TRACK"] = TokenType::TRACK;
    keywords["REGION"] = TokenType::REGION;
    keywords["PROMOTER"] = TokenType::PROMOTER;
    keywords["PROMOTERS"] = TokenType::PROMOTERS;
    keywords["ENHANCER"] = TokenType::ENHANCER;
    keywords["OF"] = TokenType::OF;
    keywords["UPSTREAM"] = TokenType::UPSTREAM;
    keywords["DOWNSTREAM"] = TokenType::DOWNSTREAM;
    keywords["WITHIN"] = TokenType::WITHIN;
    keywords["WITH"] = TokenType::WITH;
    keywords["FROM"] = TokenType::FROM;
    keywords["AT"] = TokenType::AT;
    keywords["TO"] = TokenType::TO;
    keywords["FORMAT"] = TokenType::FORMAT;
    keywords["SPACING"] = TokenType::SPACING;
    keywords["ORDER"] = TokenType::ORDER;
    keywords["ORIENTATION"] = TokenType::ORIENTATION;
    keywords["ANY"] = TokenType::ANY;
    keywords["SAME"] = TokenType::SAME;
    keywords["OPPOSITE"] = TokenType::OPPOSITE;
    keywords["AS_WRITTEN"] = TokenType::AS_WRITTEN;
    keywords["EVIDENCE"] = TokenType::EVIDENCE;
    keywords["ACCESSIBILITY"] = TokenType::ACCESSIBILITY;
    keywords["BINDING"] = TokenType::BINDING;
    keywords["OTHER"] = TokenType::OTHER;
    keywords["ASSAY"] = TokenType::ASSAY;
    keywords["SAMPLE"] = TokenType::SAMPLE;
    keywords["CONDITION"] = TokenType::CONDITION;
    keywords["REPLICATE"] = TokenType::REPLICATE;
    keywords["CONTROL"] = TokenType::CONTROL;
    keywords["TRACK_SCORE"] = TokenType::TRACK_SCORE;
    keywords["SIGNAL_VALUE"] = TokenType::SIGNAL_VALUE;
    keywords["MINUS_LOG10_PVALUE"] = TokenType::MINUS_LOG10_PVALUE;
    keywords["MINUS_LOG10_QVALUE"] = TokenType::MINUS_LOG10_QVALUE;
    keywords["EVIDENCE_CLASS"] = TokenType::EVIDENCE_CLASS;
    keywords["SUPPORT_COUNT"] = TokenType::SUPPORT_COUNT;
    keywords["BED"] = TokenType::BED;
    keywords["NARROWPEAK"] = TokenType::NARROWPEAK;
    keywords["GFF3"] = TokenType::GFF3;
    keywords["TSV"] = TokenType::TSV;
    keywords["BP"] = TokenType::BP;
    keywords["KB"] = TokenType::KB;
    keywords["MB"] = TokenType::MB;
    keywords["AND"] = TokenType::AND;
    keywords["OR"] = TokenType::OR;
    keywords["NOT"] = TokenType::NOT;
    keywords["WHERE"] = TokenType::WHERE;
    keywords["LENGTH"] = TokenType::LENGTH;
    keywords["START"] = TokenType::START;
    keywords["END"] = TokenType::END;
    keywords["MOD"] = TokenType::MOD;
    keywords["AS"] = TokenType::AS;

    keywords["EXON"] = TokenType::EXON;
    keywords["INTRON"] = TokenType::INTRON;
    keywords["UTR"] = TokenType::UTR;
    keywords["TSS"] = TokenType::TSS;
    keywords["CDS"] = TokenType::CDS;

    keywords["CHR"] = TokenType::CHR;
    keywords["CHROMOSOME"] = TokenType::CHROMOSOME;
    keywords["STRAND"] = TokenType::STRAND;
    keywords["POSITIVE"] = TokenType::POSITIVE;
    keywords["NEGATIVE"] = TokenType::NEGATIVE;

    keywords["INTERSECT"] = TokenType::INTERSECT;
    keywords["UNION"] = TokenType::UNION;
    keywords["EXCEPT"] = TokenType::EXCEPT;
    keywords["OVERLAPS"] = TokenType::OVERLAPS;
    keywords["NEAR"] = TokenType::NEAR;
    keywords["CONSENSUS"] = TokenType::CONSENSUS;
    keywords["ANCHOR"] = TokenType::ANCHOR;
    keywords["MIN_SUPPORT"] = TokenType::MIN_SUPPORT;
    keywords["MIN_RECIPROCAL_OVERLAP"] =
        TokenType::MIN_RECIPROCAL_OVERLAP;
    keywords["MAX_SUMMIT_DISTANCE"] = TokenType::MAX_SUMMIT_DISTANCE;

    keywords["SIMILARITY"] = TokenType::SIMILARITY;
    keywords["REVERSE_COMPLEMENT"] = TokenType::REVERSE_COMPLEMENT;

    keywords["MATRIX"] = TokenType::MATRIX;
    keywords["SCAN"] = TokenType::SCAN;
    keywords["THRESHOLD"] = TokenType::THRESHOLD;
    keywords["PVALUE"] = TokenType::PVALUE;
    keywords["QVALUE"] = TokenType::QVALUE;
    keywords["BACKGROUND"] = TokenType::BACKGROUND;
    keywords["UNIFORM"] = TokenType::UNIFORM;

    keywords["ANALYZE"] = TokenType::ANALYZE;
    keywords["GC_CONTENT"] = TokenType::GC_CONTENT;
    keywords["CPG_ISLANDS"] = TokenType::CPG_ISLANDS;
    keywords["WINDOW"] = TokenType::WINDOW;

    keywords["IF"] = TokenType::IF;
    keywords["THEN"] = TokenType::THEN;
    keywords["ELSE"] = TokenType::ELSE;
    keywords["ENDIF"] = TokenType::ENDIF;
    keywords["FOREACH"] = TokenType::FOREACH;
    keywords["IN"] = TokenType::IN;
    keywords["DO"] = TokenType::DO;
    keywords["ENDFOR"] = TokenType::ENDFOR;
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

    char exponent = getChar();
    if (exponent == 'e' || exponent == 'E') {
        std::string exponentText(1, exponent);
        char next = getChar();
        if (next == '+' || next == '-') {
            exponentText += next;
            next = getChar();
        }
        if (!std::isdigit(static_cast<unsigned char>(next))) {
            if (next != EOF) ungetChar();
            hasError = true;
            std::cerr << "Lexical Error at L" << startLine << ":C"
                      << startCol << " - Invalid numeric exponent in '"
                      << lexeme + exponentText << "'" << std::endl;
            return createToken(TokenType::ERROR_TOKEN,
                               lexeme + exponentText,
                               startLine, startCol);
        }
        isFloat = true;
        exponentText += next;
        while (true) {
            next = getChar();
            if (!std::isdigit(static_cast<unsigned char>(next))) {
                if (next != EOF) ungetChar();
                break;
            }
            exponentText += next;
        }
        lexeme += exponentText;
    } else if (exponent != EOF) {
        ungetChar();
    }
    return createToken(isFloat ? TokenType::FLOAT : TokenType::NUM, lexeme, startLine, startCol);
}

Token Lexer::stringLiteral(int startLine, int startCol) {
    std::string lexeme = "";
    while (true) {
        char c = getChar();
        if (c == EOF) {
            hasError = true;
            std::cerr << "Lexical Error at L" << startLine << ":C" << startCol
                      << " - Unterminated string literal" << std::endl;
            break;
        }
        if (c == '"') {
            break;
        }
        lexeme += c;
    }
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
                case '[': tokens.push_back(createToken(TokenType::LBRACKET, lexeme, startLine, startCol)); break;
                case ']': tokens.push_back(createToken(TokenType::RBRACKET, lexeme, startLine, startCol)); break;
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
                    hasError = true;
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
        case TokenType::USE: return "USE";
        case TokenType::EXPORT: return "EXPORT";
        case TokenType::DEFINE: return "DEFINE";
        case TokenType::MOTIF: return "MOTIF";
        case TokenType::MODULE: return "MODULE";
        case TokenType::GENE: return "GENE";
        case TokenType::SEQUENCE: return "SEQUENCE";
        case TokenType::ANNOTATION: return "ANNOTATION";
        case TokenType::TRACK: return "TRACK";
        case TokenType::REGION: return "REGION";
        case TokenType::PROMOTER: return "PROMOTER";
        case TokenType::PROMOTERS: return "PROMOTERS";
        case TokenType::ENHANCER: return "ENHANCER";
        case TokenType::OF: return "OF";
        case TokenType::UPSTREAM: return "UPSTREAM";
        case TokenType::DOWNSTREAM: return "DOWNSTREAM";
        case TokenType::WITHIN: return "WITHIN";
        case TokenType::WITH: return "WITH";
        case TokenType::FROM: return "FROM";
        case TokenType::AT: return "AT";
        case TokenType::TO: return "TO";
        case TokenType::FORMAT: return "FORMAT";
        case TokenType::SPACING: return "SPACING";
        case TokenType::ORDER: return "ORDER";
        case TokenType::ORIENTATION: return "ORIENTATION";
        case TokenType::ANY: return "ANY";
        case TokenType::SAME: return "SAME";
        case TokenType::OPPOSITE: return "OPPOSITE";
        case TokenType::AS_WRITTEN: return "AS_WRITTEN";
        case TokenType::EVIDENCE: return "EVIDENCE";
        case TokenType::ACCESSIBILITY: return "ACCESSIBILITY";
        case TokenType::BINDING: return "BINDING";
        case TokenType::OTHER: return "OTHER";
        case TokenType::ASSAY: return "ASSAY";
        case TokenType::SAMPLE: return "SAMPLE";
        case TokenType::CONDITION: return "CONDITION";
        case TokenType::REPLICATE: return "REPLICATE";
        case TokenType::CONTROL: return "CONTROL";
        case TokenType::TRACK_SCORE: return "TRACK_SCORE";
        case TokenType::SIGNAL_VALUE: return "SIGNAL_VALUE";
        case TokenType::MINUS_LOG10_PVALUE: return "MINUS_LOG10_PVALUE";
        case TokenType::MINUS_LOG10_QVALUE: return "MINUS_LOG10_QVALUE";
        case TokenType::EVIDENCE_CLASS: return "EVIDENCE_CLASS";
        case TokenType::SUPPORT_COUNT: return "SUPPORT_COUNT";
        case TokenType::BED: return "BED";
        case TokenType::NARROWPEAK: return "NARROWPEAK";
        case TokenType::GFF3: return "GFF3";
        case TokenType::TSV: return "TSV";
        case TokenType::BP: return "BP";
        case TokenType::KB: return "KB";
        case TokenType::MB: return "MB";
        case TokenType::AND: return "AND";
        case TokenType::OR: return "OR";
        case TokenType::NOT: return "NOT";
        case TokenType::WHERE: return "WHERE";
        case TokenType::LENGTH: return "LENGTH";
        case TokenType::START: return "START";
        case TokenType::END: return "END";
        case TokenType::MOD: return "MOD";
        case TokenType::AS: return "AS";
        case TokenType::EXON: return "EXON";
        case TokenType::INTRON: return "INTRON";
        case TokenType::UTR: return "UTR";
        case TokenType::TSS: return "TSS";
        case TokenType::CDS: return "CDS";
        case TokenType::CHR: return "CHR";
        case TokenType::CHROMOSOME: return "CHROMOSOME";
        case TokenType::STRAND: return "STRAND";
        case TokenType::POSITIVE: return "POSITIVE";
        case TokenType::NEGATIVE: return "NEGATIVE";
        case TokenType::INTERSECT: return "INTERSECT";
        case TokenType::UNION: return "UNION";
        case TokenType::EXCEPT: return "EXCEPT";
        case TokenType::OVERLAPS: return "OVERLAPS";
        case TokenType::NEAR: return "NEAR";
        case TokenType::CONSENSUS: return "CONSENSUS";
        case TokenType::ANCHOR: return "ANCHOR";
        case TokenType::MIN_SUPPORT: return "MIN_SUPPORT";
        case TokenType::MIN_RECIPROCAL_OVERLAP:
            return "MIN_RECIPROCAL_OVERLAP";
        case TokenType::MAX_SUMMIT_DISTANCE: return "MAX_SUMMIT_DISTANCE";
        case TokenType::SIMILARITY: return "SIMILARITY";
        case TokenType::REVERSE_COMPLEMENT: return "REVERSE_COMPLEMENT";
        case TokenType::MATRIX: return "MATRIX";
        case TokenType::SCAN: return "SCAN";
        case TokenType::THRESHOLD: return "THRESHOLD";
        case TokenType::PVALUE: return "PVALUE";
        case TokenType::QVALUE: return "QVALUE";
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
        case TokenType::BACKGROUND: return "BACKGROUND";
        case TokenType::UNIFORM: return "UNIFORM";
        case TokenType::END_OF_FILE: return "EOF";
        case TokenType::ERROR_TOKEN: return "ERROR_TOKEN";
        default: return "UNKNOWN";
    }
}
