CC=g++
CFLAGS=-pg -ggdb
LFLAGS=-lm -pg

HEADERS=core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h
OBJECTS=core.o data.o optimizer.o psychometric.o sigmoid.o bootstrap.o mclist.o special.o
TESTS=tests_all

testfit: $(OBJECTS) testfit.o
	$(CC) $(OBJECTS) testfit.o $(LFLAGS) -o testfit

tests: tests_all

testfit.o: testfit.cc $(HEADERS)
	$(CC) -c $(CFLAGS) testfit.cc
core.o: core.cc $(HEADERS)
	$(CC) -c $(CFLAGS) core.cc
data.o: data.cc $(HEADERS)
	$(CC) -c $(CFLAGS) data.cc
optimizer.o: optimizer.cc $(HEADERS)
	$(CC) -c $(CFLAGS) optimizer.cc
psychometric.o: psychometric.cc $(HEADERS)
	$(CC) -c $(CFLAGS) psychometric.cc
sigmoid.o: sigmoid.cc $(HEADERS)
	$(CC) -c $(CFLAGS) sigmoid.cc
mclist.o: mclist.cc $(HEADERS)
	$(CC) -c $(CFLAGS) mclist.cc
bootstrap.o: bootstrap.cc $(HEADERS)
	$(CC) -c $(CFLAGS) bootstrap.cc
special.o: special.cc $(HEADERS)
	$(CC) -c $(CFLAGS) special.cc

tests_all: tests_all.cc $(HEADERS) $(OBJECTS) testing.h
	$(CC) -c $(CFLAGS) tests_all.cc
	$(CC) $(OBJECTS) tests_all.o $(LFLAGS) -o tests_all

clean:
	rm $(OBJECTS)
