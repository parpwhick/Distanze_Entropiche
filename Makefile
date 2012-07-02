COPTS= -O3 -march=native
#COPTS= -g3
COPTS += -Wall #-fopenmp
LINK_OPTS = -lm -lgomp
files=*.cpp *.h Makefile
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

OBJ_LIST= init_functions.o distanze.o partizioni.o rand55.o adj_handler.o rand_mersenne.o 
DIST_GEN_OBJ = general_distance.o ising_simulation.o
ALL: distanze_generiche ising distanze_lineari sierpinski

distanze_generiche: ${DIST_GEN_OBJ} ${OBJ_LIST}
	g++ ${COPTS} -o distanze_generiche ${DIST_GEN_OBJ} ${OBJ_LIST} ${LINK_OPTS}

distanze_lineari: sequence_partitions.o ${OBJ_LIST}
	g++ -o distanze_lineari sequence_partitions.o ${OBJ_LIST} ${LINK_OPTS}

ising: ising_simulation.cpp adj_handler.o adj_handler.h rand_mersenne.o 
	g++ ${COPTS} -o ising -DSTANDALONE ising_simulation.cpp adj_handler.o rand_mersenne.o 

sierpinski: sierpinski.cpp
	g++ -o sierpinski -O3 -Wall sierpinski.cpp

general_distance.o: general_distance.cpp strutture.h 
	g++ ${COPTS} -c general_distance.cpp

sequence_partitions.o: sequence_partitions.cpp strutture.h 
	g++ ${COPTS} -c sequence_partitions.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

ising_simulation.o: ising_simulation.cpp ising_simulation.h adj_handler.h
	g++ ${COPTS} -c ising_simulation.cpp


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

rand55.o: rand55.cpp rand55.h
	g++ ${COPTS} -c rand55.cpp

clean:
	rm -f *.o

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}

arch: ${files}
	mkdir distanze_entropiche
	cp ${files} distanze_entropiche
	tar -cvzf prog_distanze.tar.gz distanze_entropiche
	rm -fr distanze_entropiche
