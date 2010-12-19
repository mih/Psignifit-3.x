######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

# Makefile for profile_mcmc


BUILD=build
SRC=../src/

CC=g++
CFLAGS=-pg -ggdb -Wall -I$(SRC)
LFLAGS=-lm -pg

HEADERS= $(addprefix $(SRC), psipp.h core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h mcmc.h rng.h linalg.h getstart.h)
OBJECTS= $(addprefix $(BUILD)/, core.o data.o optimizer.o psychometric.o sigmoid.o bootstrap.o mclist.o special.o mcmc.o rng.o linalg.o prior.o getstart.o)

profile: profile_mcmc

profile_mcmc: $(BUILD)/profile_mcmc.o $(OBJECTS)
	$(CC) $(OBJECTS) $(BUILD)/profile_mcmc.o -o profile_mcmc $(LFLAGS)

$(BUILD)/profile_mcmc.o: profile_mcmc.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) profile_mcmc.cc  -o $(BUILD)/profile_mcmc.o

$(BUILD):
	mkdir $(BUILD)

$(BUILD)/core.o: $(SRC)core.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)core.cc -o $(BUILD)/core.o
$(BUILD)/data.o: $(SRC)data.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)data.cc -o $(BUILD)/data.o
$(BUILD)/optimizer.o: $(SRC)optimizer.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)optimizer.cc -o $(BUILD)/optimizer.o
$(BUILD)/psychometric.o: $(SRC)psychometric.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)psychometric.cc -o $(BUILD)/psychometric.o
$(BUILD)/sigmoid.o: $(SRC)sigmoid.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)sigmoid.cc -o $(BUILD)/sigmoid.o
$(BUILD)/mclist.o: $(SRC)mclist.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)mclist.cc -o $(BUILD)/mclist.o
$(BUILD)/bootstrap.o: $(SRC)bootstrap.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)bootstrap.cc -o $(BUILD)/bootstrap.o
$(BUILD)/special.o: $(SRC)special.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)special.cc -o $(BUILD)/special.o
$(BUILD)/mcmc.o: $(SRC)mcmc.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)mcmc.cc -o $(BUILD)/mcmc.o
$(BUILD)/rng.o: $(SRC)rng.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)rng.cc -o $(BUILD)/rng.o
$(BUILD)/linalg.o: $(SRC)linalg.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)linalg.cc -o $(BUILD)/linalg.o
$(BUILD)/prior.o: $(SRC)prior.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)prior.cc -o $(BUILD)/prior.o
$(BUILD)/getstart.o: $(SRC)getstart.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)getstart.cc -o $(BUILD)/getstart.o

clean:
	-rm -rf $(BUILD)
	-rm profile_mcmc
