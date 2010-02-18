# Global Makefile for psignifit2.0
# The following commands are currently available:
# *python*:
#     build the python extension
# *python-doc*:
#     use sphinx to build a html documentation of the python module
#     in the following folder:
DOCOUT=doc-html
# *build*:
#     build everything that can be build


################# GROUPING FILES ######################
PYTHONFILES=pypsignifit/__init__.py pypsignifit/psignidata.py pypsignifit/psignierrors.py pypsignifit/psigniplot.py pypsignifit/psigobservers.py pypsignifit/pygibbsit.py
CFILES_LIB=src/bootstrap.cc src/core.cc src/data.cc src/linalg.cc src/mclist.cc src/mcmc.cc src/optimizer.cc src/psychometric.cc src/rng.cc src/sigmoid.cc src/special.cc
HFILES_LIB=src/bootstrap.h  src/core.h  src/data.h  src/errors.h src/linalg.h src/mclist.h src/mcmc.h src/optimizer.h src/prior.h src/psychometric.h src/rng.h src/sigmoid.h src/special.h src/psipp.h
CHFILES_INTERFACE=src/psipy.cc src/psipy_doc.h src/pytools.h
DOCFILES=doc-src/API.rst doc-src/index.rst doc-src/TUTORIAL.rst doc-src/*.png

################ COMMAND DEFINITIONS #################

install: python-install python-doc

build: python-build python-doc

python-install: $(PYTHONFILES) $(CFILES) $(HFILES) setup.py
	echo "Installing python extension"
	python setup.py install

python-build: $(PYTHONFILES) $(CFILES) $(HFILES) setup.py
	echo "Building python extension"
	python setup.py build_ext
	printf "The module can be used if you set\nPYTHONPATH=%s/src/\n" `pwd`

clean-python-build:
	echo "clean python build"
	rm -rv build

python-doc: $(DOCFILES) $(PYTHONFILES) python-build
	echo "building sphinx documentation"
	PYTHONPATH=build/`ls -1 build | grep lib` sphinx-build doc-src $(DOCOUT)

clean-python-doc:
	echo "clean sphinx documentation"
	rm -rv $(DOCOUT)

test-cpp:
	cd src


clean: clean-python-doc clean-python-build
