TOPROOT = $(PWD)
GASPI_DIR = /home/harry/GPI-2

BIN_DIR = $(PWD)/bin

CC =  gcc

#gcc
CFLAGS = -Wall -c -std=gnu99

SOURCES = main.c

LIBS    = -L$(GASPI_DIR)/lib64 -libverbs -lpthread -lrt -lGPI2  -Wl,-rpath=$(GASPI_DIR)/lib64:
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
