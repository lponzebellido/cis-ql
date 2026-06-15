CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
TARGET = cisql

SRCS = main.cpp Lexer.cpp Parser.cpp AST.cpp SemanticAnalyzer.cpp \
       IRGenerator.cpp Interpreter.cpp FastaReader.cpp GFFReader.cpp \
       MotifFinder.cpp SmithWaterman.cpp SetOperations.cpp
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
