#ifndef SYMBOL_TABLE_H
#define SYMBOL_TABLE_H

#include <string>
#include <unordered_map>
#include <iostream>

struct SymbolInfo {
    std::string type;
};

class SymbolTable {
private:
    std::unordered_map<std::string, SymbolInfo> table;

public:
    void insert(const std::string& name, const std::string& type) {
        if (table.find(name) == table.end()) {
            table[name] = {type};
        }
    }

    bool lookup(const std::string& name) {
        return table.find(name) != table.end();
    }

    void print() const {
        std::cout << "\n ### TABLA DE SIMBOLOS ###" << std::endl;
        for (const auto& par : table) {
            std::cout << "Clave: " << par.first << "\t| Tipo: " << par.second.type << std::endl;
        }
        std::cout << "\n" << std::endl;
    }
};

#endif
