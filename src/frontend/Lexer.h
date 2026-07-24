#ifndef LEXER_H
#define LEXER_H

#include "Token.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#define MAXLENBUF 1024

class Lexer {
private:
  std::istream &in;
  char buffer[MAXLENBUF];
  int colIndex;
  int line;
  int bufferLen;
  bool eofReached;
  bool hasError;
  std::unordered_map<std::string, TokenType> keywords;
  void initKeywords();
  char getChar();
  void ungetChar();
  void skipWhitespaceAndComments();

  Token createToken(TokenType type, const std::string &lexeme, int startLine,
                    int startCol);
  Token identifierOrKeyword(char firstChar, int startLine, int startCol);
  Token number(char firstChar, int startLine, int startCol);
  Token stringLiteral(int startLine, int startCol);

public:
  explicit Lexer(std::istream &inputStream);
  std::vector<Token> tokenize();
  bool hadError() const { return hasError; }
};

#endif
