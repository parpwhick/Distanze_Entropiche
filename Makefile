COPTS= -O3 -march=native
# debugging symbols
#COPTS+=-g2
# debugging the STL library
#COPTS+=-D_GLIBCXX_DEBUG
# profiling
#COPTS+= -pg
# STL profiling
#COPTS+= -D_GLIBCXX_PROFILE
COPTS +=-Wall -fopenmp
#COPTS +=-std=c++0x
LINK_OPTS = -lm -lgomp
files=*.cpp *.h Makefile Doxyfile
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

OBJ_LIST= init_functions.o distance.o partizioni.o adj_handler.o rand_mersenne.o 
DIST_GEN_OBJ = general_distance.o ising_simulation.o

all: distanze_generiche distanze_lineari #sierpinski ising

distanze_generiche: ${DIST_GEN_OBJ} ${OBJ_LIST}
	g++ ${COPTS} -o distanze_generiche ${DIST_GEN_OBJ} ${OBJ_LIST} ${LINK_OPTS}

distanze_lineari: sequence_partitions.o ${OBJ_LIST}
	g++ -o distanze_lineari sequence_partitions.o ${OBJ_LIST} ${LINK_OPTS}

ising: ising_simulation.cpp adj_handler.o adj_handler.h rand_mersenne.o 
	g++ ${COPTS} -o ising -DSOLO_SIMULAZIONE ising_simulation.cpp adj_handler.o rand_mersenne.o 

sierpinski: sierpinski.cpp smart_data_types.h
	g++ -o sierpinski -O3 -Wall sierpinski.cpp

general_distance.o: general_distance.cpp strutture.h partizioni.h distance.h adj_handler.h
	g++ ${COPTS} -c general_distance.cpp

sequence_partitions.o: sequence_partitions.cpp strutture.h partizioni.h
	g++ ${COPTS} -c sequence_partitions.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

ising_simulation.o: ising_simulation.cpp ising_simulation.h adj_handler.h distance.h partizioni.h smart_data_types.h
	g++ ${COPTS} -c ising_simulation.cpp


init_functions.o: strutture.h init_functions.cpp
	g++ ${COPTS} -c init_functions.cpp

distance.o: strutture.h distance.cpp distance.h
	g++ ${COPTS} -c distance.cpp

partizioni.o: strutture.h partizioni.cpp adj_handler.h partizioni.h
	g++ ${COPTS} -c partizioni.cpp

adj_handler.o: adj_handler.cpp adj_handler.h
	g++ ${COPTS} -c adj_handler.cpp
	
rand_mersenne.o: rand_mersenne.cpp rand_mersenne.h
	g++ ${COPTS} -c rand_mersenne.cpp

rand55.o: rand55.cpp rand55.h
	g++ ${COPTS} -c rand55.cpp

clean:
	rm -f *.o
	rm -vf *.bin

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}

arch: ${files}
	mkdir distanze_entropiche
	cp ${files} distanze_entropiche
	tar -cvzf prog_distanze.tar.gz distanze_entropiche
	rm -fr distanze_entropiche
