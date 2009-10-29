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
PYTHONFILES=src/pypsignifit.py src/psignidata.py src/psignierrors.py src/psigniplot.py src/psigobservers.py src/pygibbsit.py
CFILES_LIB=src/bootstrap.cc src/core.cc src/data.cc src/linalg.cc src/mclist.cc src/mcmc.cc src/optimizer.cc src/psychometric.cc src/rng.cc src/sigmoid.cc src/special.cc
HFILES_LIB=src/bootstrap.h  src/core.h  src/data.h  src/errors.h src/linalg.h src/mclist.h src/mcmc.h src/optimizer.h src/prior.h src/psychometric.h src/rng.h src/sigmoid.h src/special.h src/psipp.h
CHFILES_INTERFACE=src/psipy.cc src/psipy_doc.h src/pytools.h
DOCFILES=doc-src/API.rst doc-src/index.rst doc-src/TUTORIAL.rst doc-src/*.png

################ COMMAND DEFINITIONS #################

build: python python-doc

python: $(PYTHONFILES) $(CFILES) $(HFILES) src/setup.py
	echo "Building python extension"
	cd src/ ; python setup.py build_ext

python-doc: $(DOCFILES) $(PYTHONFILES)
	echo "building sphinx documentation"
	PYTHONPATH=src/ sphinx-build doc-src $(DOCOUT)

