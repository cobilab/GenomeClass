# Compiler and flags
CC = gcc

# Use -pg for gprof profiling
CFLAGS = -pg -g -O0 -fno-lto -march=native -funroll-loops -ffast-math -flto -Wall -Wextra -pthread
LDFLAGS = -pg -lm -pthread  # Add -pthread for linking pthread library

# Object and header dependencies
OBJ = genomeclass.o io.o mem.o buffer.o array.o hash.o misc.o pmodels.o alphabet.o context.o
# Separate utility in Samples/ that can be built as its own executable
DM_SRC = dataset_modifier.c
DM_OBJ = dataset_modifier.o
DEPS = genomeclass.h defs.h io.h rc.h msg.h mem.h buffer.h array.h hash.h misc.h pmodel.h alphabet.h context.h

# Default target: build both executables
.PHONY: all
all: genomeclass dataset_modifier

# Rule to compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Link final binary with profiling support
genomeclass: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

# Build dataset_modifier executable from Samples/dataset_modifier.c
dataset_modifier: $(DM_OBJ)
	$(CC) $(CFLAGS) $(DM_OBJ) -o $@ $(LDFLAGS)

# Rule to compile the dataset modifier source (path-aware)
$(DM_OBJ): $(DM_SRC)
	$(CC) $(CFLAGS) -c $< -o $@

# Run the program
run: genomeclass
	./genomeclass

# Profile with gprof
profile: run
	gprof genomeclass gmon.out > profile.txt
	@echo "Profile report saved to profile.txt"

# Clean build artifacts
.PHONY: clean
clean:
	rm -f *.o genomeclass dataset_modifier gmon.out profile.txt

