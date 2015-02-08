CC=g++
CFLAGS=
LDFLAGS=
SOURCES=src/*
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=bin/N-TR0PY

make: src/*
	g++ -O3 src/* -o bin/Riemann1D
