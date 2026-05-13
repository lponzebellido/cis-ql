CXX = g++
CXXFLAGS = -std=c++11 -Wall
TARGET = lexer

SRCS = main.cpp Lexer.cpp
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
