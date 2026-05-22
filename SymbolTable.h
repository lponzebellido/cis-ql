#ifndef SYMBOL_TABLE_H
#define SYMBOL_TABLE_H

#include <iostream>
#include <string>
#include <unordered_map>

struct SymbolInfo {
  std::string type;
};

class SymbolTable {
private:
  std::unordered_map<std::string, SymbolInfo> table;

public:
  void insert(const std::string &name, const std::string &type) {
    if (table.find(name) == table.end()) {
      table[name] = {type};
    }
  }

  bool lookup(const std::string &name) {
    return table.find(name) != table.end();
  }

  void print() const {
    std::cout << "\nSYMBOL TABLE" << std::endl;
    for (const auto &par : table) {
      std::cout << "Lexeme: " << par.first << "\t| Type: " << par.second.type
                << std::endl;
    }
    std::cout << "\n" << std::endl;
  }
};

#endif
