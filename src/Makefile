CC=g++
CFLAGS=-pg -ggdb -Wall
LFLAGS=-lm -pg

BUILD=build
HEADERS=core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h mcmc.h rng.h linalg.h
OBJECTS= $(addprefix $(BUILD)/, core.o data.o optimizer.o psychometric.o sigmoid.o bootstrap.o mclist.o special.o mcmc.o rng.o linalg.o prior.o)
TESTS=tests_all

tests: $(BUILD) tests_all.log

play: $(OBJECTS) play.cc $(HEADERS)
	$(CC) -c $(CFLAGS) play.cc -o $(BUILD)/play.o
	$(CC) $(OBJECTS) $(BUILD)/play.o $(LFLAGS) -o play

$(BUILD):
	mkdir $(BUILD)

$(BUILD)/core.o: core.cc $(HEADERS)
	$(CC) -c $(CFLAGS) core.cc -o $(BUILD)/core.o
$(BUILD)/data.o: data.cc $(HEADERS)
	$(CC) -c $(CFLAGS) data.cc -o $(BUILD)/data.o
$(BUILD)/optimizer.o: optimizer.cc $(HEADERS)
	$(CC) -c $(CFLAGS) optimizer.cc -o $(BUILD)/optimizer.o
$(BUILD)/psychometric.o: psychometric.cc $(HEADERS)
	$(CC) -c $(CFLAGS) psychometric.cc -o $(BUILD)/psychometric.o
$(BUILD)/sigmoid.o: sigmoid.cc $(HEADERS)
	$(CC) -c $(CFLAGS) sigmoid.cc -o $(BUILD)/sigmoid.o
$(BUILD)/mclist.o: mclist.cc $(HEADERS)
	$(CC) -c $(CFLAGS) mclist.cc -o $(BUILD)/mclist.o
$(BUILD)/bootstrap.o: bootstrap.cc $(HEADERS)
	$(CC) -c $(CFLAGS) bootstrap.cc -o $(BUILD)/bootstrap.o
$(BUILD)/special.o: special.cc $(HEADERS)
	$(CC) -c $(CFLAGS) special.cc -o $(BUILD)/special.o
$(BUILD)/mcmc.o: mcmc.cc $(HEADERS)
	$(CC) -c $(CFLAGS) mcmc.cc -o $(BUILD)/mcmc.o
$(BUILD)/rng.o: rng.cc $(HEADERS)
	$(CC) -c $(CFLAGS) rng.cc -o $(BUILD)/rng.o
$(BUILD)/linalg.o: linalg.cc $(HEADERS)
	$(CC) -c $(CFLAGS) linalg.cc -o $(BUILD)/linalg.o
$(BUILD)/prior.o: prior.cc $(HEADERS)
	$(CC) -c $(CFLAGS) prior.cc -o $(BUILD)/prior.o

tests_all: tests_all.cc $(HEADERS) $(OBJECTS) testing.h
	$(CC) -c $(CFLAGS) tests_all.cc -o $(BUILD)/tests_all.o
	$(CC) $(OBJECTS) $(BUILD)/tests_all.o $(LFLAGS) -o tests_all
tests_all.log: tests_all
	./tests_all

clean:
	-rm -rf $(BUILD)
	-rm tests_all
