SHELL=/bin/sh

CFLAGS= -c -g -Wall -I../headers
SRCS = ${wildcard *.c}
OBJS= ${SRCS:%.c=%.o}

all: graph2manifold.o graph_complement.o casson.o test_parse_util test_parse_orb test_read_snappea

test_parse_orb: test_parse_orb.o parse_orb.o casson.o parse_util.o ../code/liborb.a
	gcc test_parse_orb.o parse_orb.o parse_util.o ../code/*.o casson.o  -lm -o test_parse_orb

test_read_snappea: test_read_snappea.o parse_orb.o casson.o parse_util.o ../code/liborb.a
	gcc test_read_snappea.o parse_orb.o parse_util.o ../code/*.o casson.o  -lm -o test_read_snappea

test_parse_util: test_parse_util.o parse_util.o
	gcc test_parse_util.o parse_util.o -o test_parse_util

%.o: %.c
	gcc ${CFLAGS} $*.c -o $@

depend:
	@echo making dependency file Makefile.dependencies ...
	gcc -MM -I../headers *.c > Makefile.dependencies

include Makefile.dependencies
