#include "Lexer.h"
#include <cctype>

Lexer::Lexer(const std::string& sourceCode, SymbolTable* st) 
    : source(sourceCode), pos(0), line(1), column(1), symbolTable(st) {
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
    keywords["AND"] = TokenType::AND;
    keywords["OR"] = TokenType::OR;
    keywords["NOT"] = TokenType::NOT;
    keywords["WHERE"] = TokenType::WHERE;
    keywords["LENGTH"] = TokenType::LENGTH;
    keywords["AS"] = TokenType::AS;
}

char Lexer::advance() {
    column++;
    return source[pos++];
}

char Lexer::peek() {
    if (pos >= source.length()) return '\0';
    return source[pos];
}

char Lexer::peekNext() {
    if (pos + 1 >= source.length()) return '\0';
    return source[pos + 1];
}

void Lexer::skipWhitespace() {
    while (true) {
        char c = peek();
        if (c == ' ' || c == '\r' || c == '\t') {
            advance();
        } else if (c == '\n') {
            line++;
            column = 1;
            advance();
        } else {
            break;
        }
    }
}

Token Lexer::createToken(TokenType type, const std::string& lexeme) {
    return {type, lexeme, line, column - (int)lexeme.length()};
}

Token Lexer::identifierOrKeyword() {
    std::string lexeme = "";
    while (std::isalnum(peek()) || peek() == '_') {
        lexeme += advance();
    }

    if (keywords.find(lexeme) != keywords.end()) {
        return createToken(keywords[lexeme], lexeme);
    }
    
    symbolTable->insert(lexeme, "ID");
    return createToken(TokenType::ID, lexeme);
}

Token Lexer::number() {
    std::string lexeme = "";
    while (std::isdigit(peek())) {
        lexeme += advance();
    }
    return createToken(TokenType::NUM, lexeme);
}

Token Lexer::stringLiteral() {
    std::string lexeme = "";
    advance();
    while (peek() != '"' && peek() != '\0') {
        if (peek() == '\n') line++;
        lexeme += advance();
    }
    
    if (peek() == '"') advance(); 
    
    symbolTable->insert(lexeme, "LITERAL_STRING");
    return createToken(TokenType::STRING, "\"" + lexeme + "\"");
}

std::vector<Token> Lexer::tokenize() {
    std::vector<Token> tokens;
    
    while (pos < source.length()) {
        skipWhitespace();
        if (pos >= source.length()) break;
        
        char c = peek();
        
        if (std::isalpha(c)) {
            tokens.push_back(identifierOrKeyword());
        } else if (std::isdigit(c)) {
            tokens.push_back(number());
        } else if (c == '"') {
            tokens.push_back(stringLiteral());
        } else {
            advance();
            std::string lexeme(1, c);
            switch (c) {
                case ';': tokens.push_back(createToken(TokenType::SEMICOLON, lexeme)); break;
                case ',': tokens.push_back(createToken(TokenType::COMMA, lexeme)); break;
                case '(': tokens.push_back(createToken(TokenType::LPAREN, lexeme)); break;
                case ')': tokens.push_back(createToken(TokenType::RPAREN, lexeme)); break;
                case '=': 
                    tokens.push_back(createToken(TokenType::ASSIGN, lexeme)); 
                    break;
                case '>':
                    if (peek() == '=') {
                        lexeme += advance();
                        tokens.push_back(createToken(TokenType::GREATER_EQ, lexeme));
                    } else {
                        tokens.push_back(createToken(TokenType::GREATER, lexeme));
                    }
                    break;
                case '<':
                    if (peek() == '=') {
                        lexeme += advance();
                        tokens.push_back(createToken(TokenType::LESS_EQ, lexeme));
                    } else {
                        tokens.push_back(createToken(TokenType::LESS, lexeme));
                    }
                    break;
                default:
                    tokens.push_back(createToken(TokenType::UNKNOWN, lexeme));
                    break;
            }
        }
    }
    
    tokens.push_back(createToken(TokenType::END_OF_FILE, ""));
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
        case TokenType::AND: return "AND";
        case TokenType::OR: return "OR";
        case TokenType::NOT: return "NOT";
        case TokenType::WHERE: return "WHERE";
        case TokenType::LENGTH: return "LENGTH";
        case TokenType::AS: return "AS";
        case TokenType::ID: return "ID";
        case TokenType::NUM: return "NUM";
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
        case TokenType::END_OF_FILE: return "EOF";
        default: return "UNKNOWN";
    }
}
