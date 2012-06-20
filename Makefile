COPTS= -O3
#COPTS= -g3
COPTS += -Wall -fopenmp
files=*.cpp *.h Makefile2
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

OBJ_LIST= init_functions.o distanze.o partizioni.o rand55.o adj_handler.o rand_mersenne.o

ALL: distanze_generiche #ising #distanze_lineari

distanze_generiche: general_distance.o ${OBJ_LIST}
	g++ -o distanze_generiche general_distance.o ${OBJ_LIST} -lm -lgomp 

distanze_lineari: linear_distance.o ${OBJ_LIST}
	g++ -o distanze_lineari linear_distance.o ${OBJ_LIST} -lm -lgomp 

rand55.o: rand55.cpp rand55.h
	g++ ${COPTS} -c rand55.cpp

ising: ising_simulation.cpp adj_handler.o adj_handler.h rand_marsenne.o 
	g++ ${COPTS} -o ising -DSTANDALONE ising_simulation.cpp  adj_handler.o rand_marsenne.o

general_distance.o: general_distance.cpp strutture.h 
	g++ ${COPTS} -c general_distance.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

init_functions.o: strutture.h init_functions.cpp
	g++ ${COPTS} -c init_functions.cpp

distanze.o: strutture.h distanze.cpp
	g++ ${COPTS} -c distanze.cpp

partizioni.o: strutture.h partizioni.cpp adj_handler.h
	g++ ${COPTS} -c partizioni.cpp

adj_handler.o: adj_handler.cpp adj_handler.h
	g++ ${COPTS} -c adj_handler.cpp
	
rand_mersenne.o: rand_mersenne.cpp rand_mersenne.h
	g++ ${COPTS} -c rand_mersenne.cpp
	
simple_partition.o: strutture.h simple_partition.cpp
	g++ ${COPTS} -c simple_partition.cpp

clean:
	rm -f *.o

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}
