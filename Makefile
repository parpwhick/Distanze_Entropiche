COPTS= -Ofast
# debugging symbols
#COPTS+=-g2
# debugging the STL library
#COPTS+=-D_GLIBCXX_DEBUG
# profiling
#COPTS+= -pg
#LINK_OPTS += -pg
# STL profiling
#COPTS+= -D_GLIBCXX_PROFILE
ifeq ($(shell hostname),marcin)
   CC=icc
#else ifeq ($(shell hostname),PiorunBeskidow)
#   CC=icc
else
   CC=g++
endif

ifeq (${CC},g++)
    COPTS+=-fopenmp -Wall -march=native
    LINK_OPTS += -lgomp
    COPTS +=-std=c++0x
else
    COPTS+=-xHOST
    COPTS+=-ipo -no-prec-div 
    COPTS+=-parallel -ansi-alias -fargument-noalias
    COPTS+=-mkl -w2
    LINK_OPTS +=-mkl
#    INTEL_DEFINES =-DEIGEN_USE_MKL_ALL -I /usr/include/eigen3/ 
endif

files=*.cpp *.h Makefile Doxyfile
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

OBJ_LIST= init_functions.o distance.o partizioni.o adj_handler.o rand_mersenne.o 
DIST_GEN_OBJ = general_distance.o ising_simulation.o dopon_problem.o nagaoka_simulation.o

all: distanze_generiche distanze_lineari #sierpinski ising

distanze_generiche: ${DIST_GEN_OBJ} ${OBJ_LIST}
	${CC} -o distanze_generiche ${DIST_GEN_OBJ} ${OBJ_LIST} ${LINK_OPTS}

distanze_lineari: sequence_partitions.o ${OBJ_LIST}
	${CC} -o distanze_lineari sequence_partitions.o ${OBJ_LIST} ${LINK_OPTS}

ising: ising_simulation.cpp adj_handler.o adj_handler.h rand_mersenne.o 
	${CC} ${COPTS} -o ising -DSOLO_SIMULAZIONE ising_simulation.cpp adj_handler.o rand_mersenne.o 

sierpinski: sierpinski.cpp smart_data_types.h
	${CC} -o sierpinski -O3 -Wall sierpinski.cpp

energy_test: energy_test.o dopon_problem.o adj_handler.o
	${CC} -o energy_test energy_test.o init_functions.o adj_handler.o dopon_problem.o rand_mersenne.o ${LINK_OPTS}

general_distance.o: general_distance.cpp strutture.h partizioni.h distance.h adj_handler.h
	${CC} ${COPTS} -c general_distance.cpp

sequence_partitions.o: sequence_partitions.cpp strutture.h partizioni.h
	${CC} ${COPTS} -c sequence_partitions.cpp

translation.o: strutture.h translation.cpp
	${CC} ${COPTS} -c translation.cpp

ising_simulation.o: ising_simulation.cpp ising_simulation.h adj_handler.h distance.h partizioni.h smart_data_types.h
	${CC} ${COPTS} -c ising_simulation.cpp

nagaoka_simulation.o: nagaoka_simulation.cpp smart_data_types.h dopon_problem.h
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
	${CC} ${COPTS} ${INTEL_DEFINES} -c dopon_problem.cpp

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
	rm -f prog_distanze.tar.gz
	mkdir distanze_entropiche
	cp ${files} distanze_entropiche
	tar -cvzf prog_distanze.tar.gz distanze_entropiche
	rm -fr distanze_entropiche
