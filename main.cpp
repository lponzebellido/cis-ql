#include "Lexer.h"
#include "Parser.h"
#include <fstream>
#include <iostream>

int main() {
  std::string filename = "prueba.dsl";
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error: No se pudo abrir " << filename << std::endl;
    return 1;
  }

  SymbolTable symbolTable;
  Lexer lexer(file, &symbolTable);

  std::vector<Token> tokens = lexer.tokenize();

  std::cout << "\n\n### ANALISIS LEXICO ###\n\n";
  std::cout << "TIPO\t\t| LEXEMA\t\t| POSICION" << std::endl;
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

  // symbolTable.print(); // Opcional, comentar para no saturar la salida

  std::cout << "\n\n### ANALISIS SINTACTICO ###\n\n";
  Parser parser(tokens);
  auto ast = parser.parse();

  if (!parser.hadError()) {
      std::cout << "El analisis sintactico termino sin errores.\n\n";
      std::cout << "--- Arbol Sintactico Abstracto (AST) ---\n";
      ast->print();
  } else {
      std::cout << "\nEl analisis sintactico termino CON ERRORES.\n";
  }

  return 0;
}
