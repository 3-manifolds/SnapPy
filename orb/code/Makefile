SHELL=/bin/sh

CFLAGS= -c -g -Wall -I../headers
SRCS = ${wildcard *.c}
OBJS= ${SRCS:%.c=%.o}

../orb_solution_for_snappea_finite_triangulation: ../orb_solution_for_snappea_finite_triangulation.c liborb.a
	gcc -I../headers ../orb_solution_for_snappea_finite_triangulation.c liborb.a -lm -o../orb_solution_for_snappea_finite_triangulation

lib: liborb.a

liborb.a: ${OBJS}
	ar r liborb.a ${OBJS}
	ranlib liborb.a

%.o: %.c
	gcc ${CFLAGS} $*.c -o $@

depend:
	@echo making dependency file Makefile.dependencies ...
	gcc -MM -I../headers *.c > Makefile.dependencies

clean:
	rm -f liborb.a ../orb_solution_for_snappea_finite_triangulation *.o

include Makefile.dependencies
