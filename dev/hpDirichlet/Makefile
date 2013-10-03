INCLUDES = -Iqd/include -I../../headers
OBJECTS = hp_o31_matrices.o
all: qd $(OBJECTS)

clean:
	rm *.o

qd: qd-2.3.14.tar.gz
	tar xvfz qd-2.3.14.tar.gz
	cd qd-2.3.14 ;\
	./configure --prefix=`pwd`/../qd ;\
	make ;\
	make install

%.o : %.c
	g++ -c $< $(INCLUDES)