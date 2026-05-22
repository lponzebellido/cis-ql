#include "Lexer.h"
#include "Parser.h"
#include <fstream>
#include <iostream>

void printLexicalAnalysis(const std::vector<Token> &tokens) {
  std::cout << "\n\nLEXICAL ANALYSIS:\n\n";
  std::cout << "TYPE\t\t| LEXEME\t\t| POSITION" << std::endl;
  std::cout << "---------------------------------------------------"
            << std::endl;
  for (const Token &t : tokens) {
    if (t.type == TokenType::END_OF_FILE)
      break;

    std::string tipo = tokenTypeToString(t.type);
    std::cout << tipo;
    if (tipo.length() < 8)
      std::cout << "\t\t| ";
    else
      std::cout << "\t| ";

    std::cout << t.lexeme;
    if (t.lexeme.length() < 8)
      std::cout << "\t\t\t| ";
    else if (t.lexeme.length() < 16)
      std::cout << "\t\t| ";
    else
      std::cout << "\t| ";

    std::cout << "L" << t.line << ":C" << t.column << std::endl;
  }
}

void printSyntaxAnalysis(const Parser &parser,
                         const std::unique_ptr<ProgramNode> &ast) {
  std::cout << "\n\nSYNTAX ANALYSIS:\n\n";
  if (!parser.hadError()) {
    std::cout << "Syntax analysis completed without errors.\n\n";
    std::cout << "Abstract Syntax Tree (AST):\n";
    ast->print();
  } else {
    std::cout << "\nSyntax analysis finished WITH ERRORS.\n";
  }
}

int main() {
  std::string filename = "prueba.dsl";
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: Could not open " << filename << std::endl;
    return 1;
  }

  // 1. Lexical
  SymbolTable symbolTable;
  Lexer lexer(file, &symbolTable);
  std::vector<Token> tokens = lexer.tokenize();

  // 2. Syntax
  Parser parser(tokens);
  auto ast = parser.parse();

  printLexicalAnalysis(tokens);
  symbolTable.print();
  printSyntaxAnalysis(parser, ast);

  return 0;
}
