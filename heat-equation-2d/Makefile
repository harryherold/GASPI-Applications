TOPROOT = $(PWD)
GASPI_DIR = $(GPI2_ROOT)

BIN_DIR = $(PWD)/bin

CC =  gcc

#gcc
CFLAGS = -Wall -c -std=gnu99 -O0

SOURCES = main.c main_orig.c

LIBS    = -L$(GASPI_DIR)/lib -libverbs -lpthread -lrt -lGPI2  -Wl,-rpath=$(GASPI_DIR)/lib:
INCLUDE = -I$(GASPI_DIR)/include -I./

OBJS = $(SOURCES:%.c=%.o)

main_orig:main_orig.o
	$(CC) -o $@ main_orig.o $(LIBS)
	mkdir -p $(BIN_DIR)
	mv main_orig $(BIN_DIR)/main_orig

main:main.o
	$(CC) -o $@ main.o $(LIBS)
	mkdir -p $(BIN_DIR)
	mv main $(BIN_DIR)/main

%.o:%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $<

all : $(OBJS)

clean:
	rm -rf *~ \#_* *.o
