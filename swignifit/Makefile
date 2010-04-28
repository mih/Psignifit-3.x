######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

build: swignifit_wrap.cxx swignifit.py
	python setup.py build

swignifit_wrap.cxx: swignifit.i
	swig -c++ -python -v -I../src -I../psipy swignifit.i
swignifit.py: swignifit.i
	swig -c++ -python -v -I../src -I../psipy swignifit.i

clean:
	rm -rf build swignifit.py *.pyc swignifit_wrap.cxx

test: test-swignifit test-interface

test-swignifit:
	PYTHONPATH=build/`ls -1 build | grep lib` python swignifit_test.py
test-interface:
	PYTHONPATH=build/`ls -1 build | grep lib` python interface_test.py
