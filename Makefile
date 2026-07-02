CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 -I src/frontend -I src/backend -I src/bioinfo

TARGET = cisql

SRCS = src/main.cpp \
       src/frontend/Lexer.cpp src/frontend/Parser.cpp src/frontend/AST.cpp \
       src/backend/SemanticAnalyzer.cpp src/backend/IRGenerator.cpp src/backend/Interpreter.cpp \
       src/bioinfo/FastaReader.cpp src/bioinfo/GFFReader.cpp src/bioinfo/MotifFinder.cpp \
       src/bioinfo/SmithWaterman.cpp src/bioinfo/SetOperations.cpp src/bioinfo/PWMScanner.cpp

OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
