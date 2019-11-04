# program executable, program folder, and library folder
TARGET = run
BUILD = build
LIB := ../../GFD

# c++ compiler and compiler flags
CXX = g++
CXXFLAGS := -g -Wall -Wextra -pedantic -std=c++11


# commands for insert and delete folder
MKDIR_P ?= mkdir -p
RMDIR_R ?= $(RM) -r

# source files for the program
SOURCES := LoadMsh.cpp $(filter-out $(LIB)/Types/MpiEasy.cpp, $(wildcard $(LIB)/Types/*.cpp)) $(wildcard $(LIB)/Mesh/*.cpp) $(wildcard $(LIB)/Output/*.cpp)
OBJECTS := $(addprefix $(BUILD)/,$(patsubst %.cpp,%.o,$(patsubst $(LIB)/%,%,$(SOURCES))))

# general part of makefile
all: $(BUILD)/$(TARGET)
$(BUILD)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@
$(BUILD)/%.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(BUILD)/%.o: $(LIB)/%.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
.PHONY: clean
clean:
	$(RMDIR_R) $(BUILD)
run: $(BUILD)/$(TARGET)
	cd $(BUILD); ./$(TARGET) ../test.msh test.jrm pic.svg stat.txt; cd ../

