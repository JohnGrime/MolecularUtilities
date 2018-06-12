#
# Compiler information. The -pthread switch is needed for c++11 std::thread under some g++
#
CC     := clang++
CFLAGS := -Wall -Wextra -pedantic -std=c++11 -fPIC -O2 -pthread
LFLAGS := 

ifdef OMP
	CFLAGS += -fopenmp
endif

#
# Some default paths and file sets etc.
#
OBJ_DIR        := obj
BIN_DIR        := bin

UTIL_DIR       := Util
PROG_DIR       := Programs

UTIL_HEADERS   := $(wildcard $(UTIL_DIR)/*.h)
UTIL_FILES     := $(wildcard $(UTIL_DIR)/*.cpp)

PROG_FILES     := $(wildcard $(PROG_DIR)/*.cpp)

UTIL_OBJECTS   := $(UTIL_FILES:$(UTIL_DIR)/%.cpp=$(OBJ_DIR)/util_%.o)

PROGRAMS := $(PROG_FILES:$(PROG_DIR)/%.cpp=$(BIN_DIR)/%)

#
# Phony targets; do not produce output files, so will always execute where invoked
#
.PHONY: all clean


#
# First target = default, call it "all" as per usual. Explicitly add $(UTIL_OBJECTS)
# to dependencies; if left implicit, deleted automatically after programs built.
#
all: $(OBJ_DIR) $(BIN_DIR) $(UTIL_OBJECTS) $(PROGRAMS)
ifdef OMP
	@echo ""
	@echo "Attempted to build with OpenMP"
	@echo ""
else
	@echo ""
	@echo "You can build with OpenMP via e.g. OMP=yes on make command line (CC may need override)"
	@echo ""
endif


#
# Clean sweep. Ignore error messages if nothing to delete etc.
#
clean:
	rm -rf $(OBJ_DIR)/ $(BIN_DIR)/

#
# Enables "make print-VARIABLE" on command line to print $(VARIABLE) to stdout.
# This is handy for debugging the makefile if something weird is happening.
#
print-%:
	@echo $* = \"$($*)\"

#
# Required output directories
#
$(OBJ_DIR):
	mkdir $@
$(BIN_DIR):
	mkdir $@


#
# --- Util objects ---
#
# All util objects ( obj/util_x.o ) should have dependencies on $(UTIL_DIR)/x.cpp
# and all util headers.
#
# A general rule for this is given below. Note that dependency order is important:
# by placing the .cpp file as the FIRST dependency, only that file is passed to the
# compiler. This first dependency is referenced as "$<" ("$?" == all dependencies).
#
$(OBJ_DIR)/util_%.o: $(UTIL_DIR)/%.cpp $(UTIL_HEADERS)
	$(CC) -c $(CFLAGS) $< -o $@


#
# --- Programs ---
#
$(BIN_DIR)/%: $(PROG_DIR)/%.cpp $(UTIL_OBJECTS)
	$(CC)    $(CFLAGS) $(LFLAGS) $(UTIL_OBJECTS) $< -o $@
