swig:
	swig -c++ -python -I../src swignifit.i
	python setup.py build

clean:
	rm -rf build swignifit.py swignifit.pyc swignifit_wrap.cxx 
test:
	PYTHONPATH=build/`ls -1 build | grep lib` python test.py
