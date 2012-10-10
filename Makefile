COPTS= -O3 -march=native
# debugging symbols
#COPTS+=-g2
# debugging the STL library
#COPTS+=-D_GLIBCXX_DEBUG
# profiling
#COPTS+= -pg
# STL profiling
#COPTS+= -D_GLIBCXX_PROFILE
ifeq ($(shell hostname),Atomowy-linux)
   CC=g++
else
   CC=icc
endif

COPTS +=-Wall

ifeq (${CC},g++)
    COPTS+=-fopenmp
    LINK_OPTS = -lgomp
else
    COPTS+=-xHOST -ipo -no-prec-div
endif
COPTS +=-std=c++0x

files=*.cpp *.h Makefile Doxyfile
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

OBJ_LIST= init_functions.o distance.o partizioni.o adj_handler.o rand_mersenne.o 
DIST_GEN_OBJ = general_distance.o ising_simulation.o

all: distanze_generiche distanze_lineari #sierpinski ising

distanze_generiche: ${DIST_GEN_OBJ} ${OBJ_LIST}
	${CC} -o distanze_generiche ${DIST_GEN_OBJ} ${OBJ_LIST} ${LINK_OPTS}

distanze_lineari: sequence_partitions.o ${OBJ_LIST}
	${CC} -o distanze_lineari sequence_partitions.o ${OBJ_LIST} ${LINK_OPTS}

ising: ising_simulation.cpp adj_handler.o adj_handler.h rand_mersenne.o 
	${CC} ${COPTS} -o ising -DSOLO_SIMULAZIONE ising_simulation.cpp adj_handler.o rand_mersenne.o 

sierpinski: sierpinski.cpp smart_data_types.h
	${CC} -o sierpinski -O3 -Wall sierpinski.cpp

general_distance.o: general_distance.cpp strutture.h partizioni.h distance.h adj_handler.h
	${CC} ${COPTS} -c general_distance.cpp

sequence_partitions.o: sequence_partitions.cpp strutture.h partizioni.h
	${CC} ${COPTS} -c sequence_partitions.cpp

translation.o: strutture.h translation.cpp
	${CC} ${COPTS} -c translation.cpp

ising_simulation.o: ising_simulation.cpp ising_simulation.h adj_handler.h distance.h partizioni.h smart_data_types.h
	${CC} ${COPTS} -c ising_simulation.cpp

nagaoka_simulation.o: nagaoka_simulation.cpp smart_data_types.h
	${CC} ${COPTS} -c nagaoka_simulation.cpp

init_functions.o: strutture.h init_functions.cpp
	${CC} ${COPTS} -c init_functions.cpp

distance.o: strutture.h distance.cpp distance.h
	${CC} ${COPTS} -c distance.cpp

partizioni.o: strutture.h partizioni.cpp adj_handler.h partizioni.h
	${CC} ${COPTS} -c partizioni.cpp

adj_handler.o: adj_handler.cpp adj_handler.h
	${CC} ${COPTS} -c adj_handler.cpp
	
rand_mersenne.o: rand_mersenne.cpp rand_mersenne.h
	${CC} ${COPTS} -c rand_mersenne.cpp

dopon_problem.o: dopon_problem.h dopon_problem.cpp adj_handler.h
	${CC} ${COPTS} ${INTEL_DEFINES} -I /usr/include/eigen3/ -I/usr/include/arpack++/ -c dopon_problem.cpp

rand55.o: rand55.cpp rand55.h
	${CC} ${COPTS} -c rand55.cpp

clean: clean_temp_files
	rm -f *.o

clean_temp_files:
	rm -f *.bin
	rm -f *.txt

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}

arch: ${files}
	mkdir distanze_entropiche
	cp ${files} distanze_entropiche
	tar -cvzf prog_distanze.tar.gz distanze_entropiche
	rm -fr distanze_entropiche
