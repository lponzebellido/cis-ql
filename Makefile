CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -MMD -MP -I src/frontend -I src/backend -I src/bioinfo

TARGET = cisql
CORE_TEST_TARGET = tests/test_core
CORE_BENCHMARK_TARGET = benchmarks/benchmark_core

SRCS = src/main.cpp \
       src/frontend/Lexer.cpp src/frontend/Parser.cpp src/frontend/AST.cpp \
       src/backend/SemanticAnalyzer.cpp src/backend/IRGenerator.cpp src/backend/Interpreter.cpp \
       src/bioinfo/FastaReader.cpp src/bioinfo/GFFReader.cpp src/bioinfo/MotifFinder.cpp \
       src/bioinfo/SmithWaterman.cpp src/bioinfo/SetOperations.cpp src/bioinfo/PWMScanner.cpp \
       src/bioinfo/GCAnalyzer.cpp

OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(DEPS) $(TARGET) $(CORE_TEST_TARGET) $(CORE_BENCHMARK_TARGET)

$(CORE_TEST_TARGET): tests/test_core.cpp \
       src/bioinfo/MotifFinder.cpp src/bioinfo/SmithWaterman.cpp \
       src/bioinfo/SetOperations.cpp src/bioinfo/PWMScanner.cpp \
       src/bioinfo/GCAnalyzer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

test: $(TARGET) $(CORE_TEST_TARGET)
	./$(CORE_TEST_TARGET)
	python3 tests/test_language.py

benchmark: $(TARGET)
	python3 benchmarks/run_examples.py

$(CORE_BENCHMARK_TARGET): benchmarks/benchmark_core.cpp \
       src/bioinfo/MotifFinder.cpp src/bioinfo/SmithWaterman.cpp \
       src/bioinfo/SetOperations.cpp src/bioinfo/PWMScanner.cpp \
       src/bioinfo/GCAnalyzer.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

benchmark-core: $(CORE_BENCHMARK_TARGET)
	./$(CORE_BENCHMARK_TARGET)

.PHONY: all run clean test benchmark benchmark-core

-include $(DEPS)
