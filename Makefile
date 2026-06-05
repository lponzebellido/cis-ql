CXX = g++
CXXFLAGS = -std=c++11 -Wall
TARGET = cisql

SRCS = main.cpp Lexer.cpp Parser.cpp AST.cpp SemanticAnalyzer.cpp
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
