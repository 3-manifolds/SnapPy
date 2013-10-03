INCLUDES = -Iqd/include -I../../headers
OBJECTS = hp_o31_matrices.o hp_Dirichlet.o hp_Dirichlet_construction.o \
hp_Dirichlet_basepoint.o hp_Dirichlet_extras.o

all: qd $(OBJECTS)

clean:
	rm *.o

qd: qd-2.3.14.tar.gz
	tar xvfz qd-2.3.14.tar.gz
	cd qd-2.3.14 ;\
	./configure --prefix=`pwd`/../qd ;\
	make ;\
	make install

%.o : %.c hp_Dirichlet.h
	g++ -c $< $(INCLUDES)