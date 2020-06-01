## MIT License
## 
## Copyright (c) 2020 Alon Rashelbach
## 
## Permission is hereby granted, free of charge, to any person obtaining a copy
## of this software and associated documentation files (the "Software"), to deal
## in the Software without restriction, including without limitation the rights
## to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
## copies of the Software, and to permit persons to whom the Software is
## furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all
## copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
## IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
## FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
## AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
## LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
## OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
## SOFTWARE.

SRC_DIR		?=		src
BIN_DIR		?=		bin
TOOLS_DIR	?=		tools
VENDOR_DIR	?=		vendor
INCLUDE_DIR	?=		include

CXX			?=		g++
CXXFLAGS	?=		-std=c++11 -pthread -Wall -fPIC -DNOPYTHON
DBGFLAGS 	?= 		-g
SIMDFLAGS	?=		-mavx2 -mfma
INCLUDES	?= 		-I $(BIN_DIR) -I $(INCLUDE_DIR) -I tuplemerge -I vendor
LIBRARIES	?=		-L $(BIN_DIR)
RM			?=		rm -f

export SRC_DIR
export BIN_DIR
export INCLUDE_DIR

# Make the bin directory
$(shell mkdir -p $(BIN_DIR) 2> /dev/null )

# Generate a list of all source files to be compiles as objects
SOURCES			:=$(wildcard $(TOOLS_DIR)/*.cpp)
SOURCES			+=$(wildcard $(SRC_DIR)/*.cpp)

# Generate a list of objects
OBJECTS 		+=$(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))
OBJECTS 		+=$(patsubst $(VENDOR_DIR)/%.cpp,$(BIN_DIR)/%.o,$(wildcard $(VENDOR_DIR)/*.cpp))

# Generate a list of executable files
EXECUTABLES		:=$(patsubst $(TOOLS_DIR)/%.cpp,$(BIN_DIR)/%.exe,$(wildcard $(TOOLS_DIR)/*.cpp))

# Generate all Executables
release: $(EXECUTABLES) python
	@echo "Done making release"

debug: $(EXECUTABLES) python
	@echo "Done making debug"

# Rules to create all objects
include $(BIN_DIR)/objects.mk

# Any executable file requires all object files!
$(BIN_DIR)/%.exe: $(OBJECTS) $(BIN_DIR)/%.o
	$(CXX) $(CXXFLAGS) $(SIMDFLAGS) $(DBGFLAGS) $(OFFLAGS) $(INCLUDES) $(LIBRARIES) $+ -o $@ -ltuplemerge

# Create librqrmi.a for python extension
librqrmi.a: $(OBJECTS)
	@ar crf $(BIN_DIR)/librqrmi.a \
	$(BIN_DIR)/algorithms.o $(BIN_DIR)/argument_handler.o $(BIN_DIR)/cpu_core_tools.o \
	$(BIN_DIR)/logging.o $(BIN_DIR)/lookup.o $(BIN_DIR)/matrix_operations.o \
	$(BIN_DIR)/rqrmi_fast.o $(BIN_DIR)/rqrmi_model.o $(BIN_DIR)/rqrmi_tools.o \
	$(BIN_DIR)/object_io.o $(BIN_DIR)/python_library.o $(BIN_DIR)/simd_aux.o \
	$(BIN_DIR)/vector_list.o

# Python file
python: librqrmi.a
	@cp $(SRC_DIR)/*.py $(BIN_DIR)/
	python3 $(SRC_DIR)/setup.py build_ext -b $(BIN_DIR) --build-lib $(BIN_DIR) -t $(BIN_DIR)

# Target specific variables
release:	DBGFLAGS = -O2 -DNDEBUG
debug:		DBGFLAGS = -O0 -g

# Clean bin directory
clean:
	$(RM) $(BIN_DIR)/*.o $(BIN_DIR)/*.exe
