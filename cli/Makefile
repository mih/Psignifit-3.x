CC=g++
CFLAGS=-pg -ggdb -Wall
LFLAGS=-lm -lpsipp -pg
TODAY=`date +%d-%m-%G`
LONGTODAY=`date +%G-%m-%d`
GIT_DESCRIPTION=`git describe --tags`
CLI_VERSION_HEADER=cli_version.h

BUILD=build

SRC=../src
HEADERS= $(addprefix $(SRC)/, core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h mcmc.h rng.h linalg.h getstart.h )
CLI_H= cli.h cli_utilities.h
CLI_O= $(addprefix $(BUILD)/, cli.o cli_utilities.o)

compile: cli-version psignifit-mapestimate psignifit-bootstrap psignifit-mcmc psignifit-diagnostics

clean:
	-rm -r $(BUILD)

cli-version:
	if git rev-parse &> /dev/null ; then \
		echo "#ifndef CLI_VERSION_H" > $(CLI_VERSION_HEADER) ; \
		echo "#define CLI_VERSION_H" >> $(CLI_VERSION_HEADER) ; \
		echo "#define VERSION \""$(GIT_DESCRIPTION)"\"" >> $(CLI_VERSION_HEADER) ; \
		echo "#endif" >> $(CLI_VERSION_HEADER) ; \
	fi

psignifit-mapestimate: $(HEADERS) $(BUILD)/psignifit-mapestimate.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-mapestimate $(OBJECTS) $(BUILD)/psignifit-mapestimate.o $(CLI_O) $(LFLAGS)
psignifit-bootstrap: $(HEADERS) $(BUILD)/psignifit-bootstrap.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-bootstrap $(OBJECTS) $(BUILD)/psignifit-bootstrap.o $(CLI_O) $(LFLAGS)
psignifit-mcmc: $(HEADERS) $(BUILD)/psignifit-mcmc.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-mcmc $(OBJECTS) $(BUILD)/psignifit-mcmc.o $(CLI_O) $(LFLAGS)
psignifit-diagnostics: $(HEADERS) $(BUILD)/psignifit-diagnostics.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-diagnostics $(OBJECTS) $(BUILD)/psignifit-diagnostics.o $(CLI_O) $(LFLAGS)

$(BUILD)/psignifit-mapestimate.o: psignifit-mapestimate.cc $(HEADERS) $(CLI_H) $(BUILD)
	$(CC) -c $(CFLAGS) psignifit-mapestimate.cc -o $(BUILD)/psignifit-mapestimate.o
$(BUILD)/psignifit-bootstrap.o: psignifit-bootstrap.cc $(HEADERS) $(CLI_H) $(BUILD)
	$(CC) -c $(CFLAGS) psignifit-bootstrap.cc -o $(BUILD)/psignifit-bootstrap.o
$(BUILD)/psignifit-mcmc.o: psignifit-mcmc.cc $(HEADERS) $(CLI_H) $(BUILD)
	$(CC) -c $(CFLAGS) psignifit-mcmc.cc -o $(BUILD)/psignifit-mcmc.o
$(BUILD)/psignifit-diagnostics.o: psignifit-diagnostics.cc $(HEADERS) $(CLI_H) $(BUILD)
	$(CC) -c $(CFLAGS) psignifit-diagnostics.cc -o $(BUILD)/psignifit-diagnostics.o

$(BUILD)/cli.o: cli.cc cli.h
	$(CC) -c $(CFLAGS) cli.cc -o $(BUILD)/cli.o
$(BUILD)/cli_utilities.o: cli_utilities.cc cli_utilities.h
	$(CC) -c $(CFLAGS) cli_utilities.cc -o $(BUILD)/cli_utilities.o

$(BUILD):
	mkdir $(BUILD)


