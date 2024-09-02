
all: main.o tdigest.o
	gcc -o a.out main.o tdigest.o -lstdc++

main.o: main.cc tdigest.h
	gcc -c main.cc -Wall -Wpedantic -Wconversion

tdigest.o: tdigest.cc tdigest.h
	gcc -c tdigest.cc -Wall -Wpedantic -Wconversion

clean:
	rm -f *.o

