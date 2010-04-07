######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

swig:
	swig -c++ -python -v -I../src swignifit.i
	python setup.py build

clean:
	rm -rf build swignifit.py *.pyc swignifit_wrap.cxx
test:
	PYTHONPATH=build/`ls -1 build | grep lib` python swignifit_test.py
