CFLAGS = -Wall -g -DINLINE
MFLAGS = $(CFLAGS) -lm
CC = gcc $(CFLAGS)
MCC = gcc $(MFLAGS)

nrlmsise-test    :     nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o
	$(MCC) -o nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o \
                     nrlmsise-00_data.o

nrlmsise-00.o :		nrlmsise-00.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00.c

nrlmsise-00_test.o :	nrlmsise-00_test.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_test.c

nrlmsise-00_data.o :	nrlmsise-00_data.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_data.c

clean   :
	rm -rf  nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o \
                   nrlmsise-00_data.o
