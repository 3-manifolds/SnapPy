CFLAGS   = -O3 -Iqd/include -I../../headers
LFLAGS   = -Lqd/lib
LIBS     = -lqd -lm
SOURCES  = hp_sl2c_matrices.cpp hp_matrix_conversion.cpp hp_o31_matrices.cpp \
hp_transcendentals.cpp hp_Dirichlet.cpp hp_Dirichlet_construction.cpp \
hp_Dirichlet_basepoint.cpp hp_Dirichlet_extras.cpp hp_complex_length.cpp \
#hp_length_spectrum.cpp
HEADERS  = hp_Dirichlet.h
OBJECTS  = hp_sl2c_matrices.o hp_matrix_conversion.o hp_o31_matrices.o \
hp_transcendentals.o hp_Dirichlet_extras.o hp_Dirichlet_construction.o \
hp_Dirichlet_basepoint.o hp_Dirichlet.o hp_complex_length.o \
#hp_length_spectrum.o
SNAPPEA = $(wildcard ../../build/temp.*/kernel_code/*.o) fakeUI.o
CC = gcc
CPP = g++


all: qd $(OBJECTS) hp_test fakeUI.o

clean:
	rm -f *.o hp_test

qd: qd-2.3.14.tar.gz
	tar xvfz qd-2.3.14.tar.gz
	cd qd-2.3.14 ;\
	./configure --prefix=`pwd`/../qd ;\
	make ; \
	make install

fake_UI.o: fakeUI.c
	$(CC) -c $(CFLAGS) fakeUI.c

hp_test.o: hp_test.cpp
	$(CPP) -c $(CFLAGS) hp_test.cpp

hp_test: qd $(OBJECTS) hp_test.o fakeUI.o
	$(CPP) -o hp_test $(CFLAGS) $(LFLAGS) hp_test.o $(OBJECTS) $(SNAPPEA) $(LIBS)

%.o : %.c $(HEADERS)
	$(CC) -c $(CFLAGS) $<

%.o : %.cpp $(HEADERS)
	$(CPP) -c $(CFLAGS) $<
