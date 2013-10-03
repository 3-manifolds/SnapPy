INCLUDES = -Iqd/include -I../../headers
LIBS     =  -Lqd/lib hp_test.cpp
LFLAGS   = -lqd
SOURCES  = hp_o31_matrices.c hp_transcendentals.c hp_Dirichlet.c \
hp_Dirichlet_construction.c hp_Dirichlet_basepoint.c hp_Dirichlet_extras.c
HEADERS  = hp_Dirichlet.h
OBJECTS  = hp_o31_matrices.o hp_transcendentals.o hp_Dirichlet.o \
hp_Dirichlet_construction.o hp_Dirichlet_basepoint.o hp_Dirichlet_extras.o

all: qd $(OBJECTS) hp_test.cpp

clean:
	rm -f *.o hp_test

qd: qd-2.3.14.tar.gz
	tar xvfz qd-2.3.14.tar.gz
	cd qd-2.3.14 ;\
	./configure --prefix=`pwd`/../qd ;\
	make ;\
	make install

hp_test: qd $(OBJECTS) hp_test.cpp
	 g++ -o hp_test $(INCLUDES) $(LIBS) $(LFLAGS)

%.o : %.c $(HEADERS)
	g++ -c $< $(INCLUDES)