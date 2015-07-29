GASPI_DIR = $(GPI2_ROOT)

CC = gcc

CFLAGS = -std=gnu99 -Wall -c -g3 -O0

LIBS    = -L$(GASPI_DIR)/lib -libverbs -lpthread -lGPI2 -Wl,-rpath=$(GASPI_DIR)/lib:

INCLUDE = -I$(GASPI_DIR)/include

OBJS = test_utils.o gaspi_utils.o

test-utils:$(OBJS)
	$(CC) -o $@ $(OBJS) $(LIBS)

%.o:%.c
	$(CC) $(CFLAGS) $(INCLUDE) $<

clean:
	rm -rf *~ \#_* *.o
	rm test-utils
