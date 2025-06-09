# Compiler and flags
CC = gcc

# Use -pg for gprof profiling
CFLAGS = -pg -g -O0 -fno-lto -march=native -funroll-loops -ffast-math -flto -Wall -Wextra
LDFLAGS = -pg -lm  # Include -pg here too for proper linking

# Object and header dependencies
OBJ = genomeclass.o io.o mem.o buffer.o array.o hash.o misc.o pmodels.o alphabet.o context.o
DEPS = genomeclass.h defs.h io.h rc.h msg.h mem.h buffer.h array.h hash.h misc.h pmodel.h alphabet.h context.h

# Rule to compile .c to .o
%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

# Link final binary with profiling support
genomeclass: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LDFLAGS)

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
	rm -f *.o genomeclass gmon.out profile.txt
