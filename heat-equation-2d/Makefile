TOPROOT = $(PWD)
GASPI_DIR = $(GPI2_ROOT)

BIN_DIR = $(PWD)/bin

CC =  gcc

#gcc
CFLAGS = -Wall -c -std=gnu99

SOURCES = main.c test.c

LIBS    = -L$(GASPI_DIR)/lib -libverbs -lpthread -lrt -lGPI2  -Wl,-rpath=$(GASPI_DIR)/lib:
INCLUDE = -I$(GASPI_DIR)/include -I./

OBJS = $(SOURCES:%.c=%.o)

main:$(OBJS)
	$(CC) -o $@ $(OBJS) $(LIBS)
	mkdir -p $(BIN_DIR)
	mv main $(BIN_DIR)/main

%.o:%.c
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $<

all : $(OBJS)

clean:
	rm -rf *~ \#_* *.o
