# compiler name (gcc seems to use smarter tricks since the executable is faster when using gcc instead of clang on an arch linux machine)
CC = gcc

SRC_DIR := src/
HEADERS_DIR := src/headers/

# build debug as default
OBJ_DIR = obj/debug/
BIN_DIR = bin/debug/
CFLAGS = -Wall -g -mavx2 -Isrc/headers
LDFLAGS = -lm

# condition to check value passed
ifeq ($(MODE),exec)
OBJ_DIR = obj/exec/
BIN_DIR = bin/exec/
CFLAGS = -O3 -ftree-loop-im -ffast-math -mavx2 -march=native -mtune=native -Isrc/headers
endif

SOURCE_NAMES = main.c paralhullIO.c quickhull.c

HEADER_NAMES = paralhull.h

# files list
HEADER_FILES := $(HEADER_NAMES:%=$(HEADERS_DIR)%)

OBJ_FILES := $(SOURCE_NAMES:%.c=$(OBJ_DIR)%.o)

# command used to check variables value and "debug" the makefile
print:
	@echo HEADER_FILES = $(HEADER_FILES)
	@echo SOURCE_NAMES = $(SOURCE_NAMES:%=$(SRC_DIR)%)
	@echo OBJ_FILES = $(OBJ_FILES)

# build options when debugging
build: $(BIN_DIR)main

$(BIN_DIR)main: $(OBJ_FILES)
	$(CC) $(CFLAGS) $(OBJ_FILES) -o $(BIN_DIR)main $(LDFLAGS)

$(OBJ_DIR)%.o: $(SRC_DIR)%.c $(HEADER_FILES)
	$(CC) -c $(CFLAGS) $(SRC_DIR)$(*F).c -o $@

SRC_FILES_PATH := $(SOURCE_NAMES:%=$(SRC_DIR)%)

final:
	$(CC) -O3 -ftree-loop-im -mavx2 -march=native -mtune=native -Isrc/headers $(SRC_FILES_PATH) -o bin/exec/main $(LDFLAGS)

# delete all gcc output files
clean:
	rm -f bin/debug/main bin/exec/main
	rm -f obj/debug/*.o obj/exec/*.o
