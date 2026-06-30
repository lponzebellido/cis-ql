#include "backend/IRGenerator.h"
#include "backend/Interpreter.h"
#include "backend/SemanticAnalyzer.h"
#include "frontend/Lexer.h"
#include "frontend/Parser.h"
#include <fstream>
#include <iomanip>
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

void printSemanticAnalysis(const SemanticAnalyzer &semantic) {
  std::cout << "\n\nSEMANTIC ANALYSIS:\n\n";
  if (!semantic.hadError()) {
    std::cout
        << "Semantic analysis completed without errors. Program is valid.\n";
  } else {
    std::cout << "\nSemantic analysis finished WITH ERRORS.\n";
  }
}

void printIRCode(const std::vector<IRInstruction> &ir) {
  std::cout << "\n\nINTERMEDIATE CODE:\n\n";
  for (size_t i = 0; i < ir.size(); i++) {
    std::cout << "[" << std::setw(3) << i << "] "
              << irOpcodeToString(ir[i].opcode);
    if (!ir[i].arg1.empty())
      std::cout << "\t" << ir[i].arg1;
    if (!ir[i].arg2.empty())
      std::cout << "\t" << ir[i].arg2;
    if (!ir[i].arg3.empty())
      std::cout << "\t" << ir[i].arg3;
    if (!ir[i].arg4.empty())
      std::cout << "\t" << ir[i].arg4;
    if (!ir[i].arg5.empty())
      std::cout << "\t" << ir[i].arg5;
    std::cout << std::endl;
  }
}

int main(int argc, char *argv[]) {
  std::string filename = "";
  bool debugMode = false;

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--debug" || arg == "-d") {
      debugMode = true;
    } else {
      filename = arg;
    }
  }

  if (filename.empty()) {
    std::cerr << "Usage: " << argv[0] << " <filename.cql> [--debug]"
              << std::endl;
    return 1;
  }

  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open " << filename << std::endl;
    return 1;
  }

  // 1. Lexical Analysis
  SymbolTable symbolTable;
  Lexer lexer(file, &symbolTable);
  std::vector<Token> tokens = lexer.tokenize();

  if (debugMode)
    printLexicalAnalysis(tokens);

  // 2. Syntax Analysis
  Parser parser(tokens);
  auto ast = parser.parse();

  if (!parser.hadError()) {
    if (debugMode)
      printSyntaxAnalysis(parser, ast);

    // 3. Semantic Analysis
    SemanticAnalyzer semantic(symbolTable);
    semantic.analyze(ast.get());

    if (debugMode)
      printSemanticAnalysis(semantic);

    if (!semantic.hadError()) {
      // 4. IR Generation
      IRGenerator irGen;
      const auto &ir = irGen.generate(ast.get());
      if (debugMode)
        printIRCode(ir);

      // 5. Execution
      Interpreter interpreter;
      interpreter.execute(ir, debugMode);
    }
  } else {
    if (debugMode)
      printSyntaxAnalysis(parser, ast);
  }

  if (debugMode)
    symbolTable.print();

  return 0;
}
