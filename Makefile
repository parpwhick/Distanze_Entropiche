COPTS= -O3 -Wall -fopenmp
files=*.cpp *.h Makefile2
file_supporto=./utility/carica* ./utility/cluster* ./utility/comandi*

distanze: main.o init_functions.o distanze.o partizioni.o rand55.o 
	g++ -o distanze main.o init_functions.o distanze.o partizioni.o rand55.o -lm -lgomp 

ising: ising.o distanze.o partizioni.o init_functions.o translation.o rand55.o
	g++ -o ising ising.o distanze.o partizioni.o init_functions.o translation.o rand55.o

rand55.o: rand55.cpp rand55.h
	g++ ${COPTS} -c rand55.cpp

ising.o: ising.cpp strutture.h
	g++ ${COPTS} -c ising.cpp

main.o: main.cpp strutture.h 
	g++ ${COPTS} -c main.cpp

translation.o: strutture.h translation.cpp
	g++ ${COPTS} -c translation.cpp

init_functions.o: strutture.h init_functions.cpp
	g++ ${COPTS} -c init_functions.cpp

distanze.o: strutture.h distanze.cpp
	g++ ${COPTS} -c distanze.cpp

partizioni.o: strutture.h partizioni.cpp
	g++ ${COPTS} -c partizioni.cpp

simple_partition.o: strutture.h simple_partition.cpp
	g++ ${COPTS} -c simple_partition.cpp

clean:
	rm -f *.o

zip: ${files} ${file_supporto}
	zip -9 prog_distanze.zip ${files} ${file_supporto}
