#include <Python.h>
#include <Numeric/arrayobject.h>
#include <vector>
#include <cstdio>

#include "psipp.h"

static char psibootstrap_doc [] =
"Bootstrap\n";

static PyObject * psibootstrap ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int i,j;
	int Nsamples ( 2000 );             // Number of bootstrap samples
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pydata;                  // python object holding the data
	PyObject *pystart (Py_None);       // python object holding the starting values
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization
	PyObject *pypriors (Py_None);      // prior specs

	PyObject *pyblock;       // python object holding data from a single block

	/************************************************************
	 * Parse command line
	 */
	// TODO: command line parsing should go to a separate function
	static char *kwlist[] = {
		"data",
		"start",
		"nsamples",
		"nafc",
		"sigmoid",
		"core",
		"priors",
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|OiissO",
				kwlist,
				&pydata,&pystart,&Nsamples,&Nafc,&sigmoidname,&corename,&pypriors ) )
		return NULL;

	/************************************************************
	 * Check data
	 */
	if ( !PySequence_Check ( pydata ) ) {
		// data are no sequence ~> Error
		PyErr_Format ( PyExc_ValueError, "data should be a sequence" );
		return NULL;
	}

	// Get the data
	int Nblocks ( PySequence_Size ( pydata ) );
	std::vector<double> x ( Nblocks );
	std::vector<int> k ( Nblocks );
	std::vector<int> n ( Nblocks );
	for ( i=0; i<Nblocks; i++ ) {
		pyblock = PySequence_GetItem ( pydata, i );
		if ( PySequence_Size ( pyblock ) != 3 ) {
			PyErr_Format ( PyExc_ValueError, "data in block %d do not have 3 entries", i );
			return NULL;
		}
		x[i] = PyFloat_AsDouble ( PySequence_GetItem ( pyblock, 0 ) );
		k[i] = PyInt_AsLong ( PySequence_GetItem ( pyblock, 1 ) );
		n[i] = PyInt_AsLong ( PySequence_GetItem ( pyblock, 2 ) );
		std::cerr << i << " " << x[i] << " " << k[i] << " " << n[i] << "\n";
	}

	/************************************************************
	 * Now set up the model and go
	 */
	PsiData * data = new PsiData (x,n,k,Nafc);

	// Determine sigmoid
	PsiSigmoid * sigmoid;
	if ( !strcmp(sigmoidname,"logistic") ) {
		std::cerr << "Using logistic sigmoid\n";
		sigmoid = new PsiLogistic;
	} else if ( !strcmp(sigmoidname,"gauss") ) {
		std::cerr << "Using gaussian cdf sigmoid\n";
		sigmoid = new PsiGauss;
	} else if ( !strcmp(sigmoidname,"gumbel_l") || !strcmp(sigmoidname,"lgumbel") ) {
		std::cerr << "Using gumbelL sigmoid\n";
		sigmoid = new PsiGumbelL;
	} else if ( !strcmp(sigmoidname,"gumbel_r") || !strcmp(sigmoidname,"rgumbel") ) {
		std::cerr << "Using gumbelR sigmoid\n";
		sigmoid = new PsiGumbelR;
	} else {
		PyErr_Format ( PyExc_ValueError, "invalid sigmoid type" );
		return NULL;
	}

	// Determine core
	PsiCore * core;
	if ( !strcmp(corename,"ab") ) {
		std::cerr << "Using core ab\n";
		core = new abCore;
	} else if ( !strncmp(corename,"mw",2) ) {
		double alpha;
		std::cerr << corename << "\n";
		if ( sscanf ( corename, "mw%lf", &alpha )==0 )
			alpha = 0.1;
		std::cerr << "Using core mw with parameter " << alpha << "\n";
		core = new mwCore ( sigmoid->getcode(), alpha );
	} else if ( !strcmp(corename,"linear") ) {
		std::cerr << "Using linear core\n";
		core = new linearCore;
	} else if ( !strcmp(corename,"log") || !strcmp(corename,"logarithmic") ) {
		std::cerr << "Using logarithmic core\n";
		core = new logCore ( data );
	} else {
		PyErr_Format ( PyExc_ValueError, "invalid core type" );
		return NULL;
	}

	PsiPsychometric * pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	int Nparams = pmf->getNparams ();
	std::vector<double> cuts (1, .5);

	// Set priors
	double priorpars[10];
	if ( pypriors == Py_None ) {
		std::cerr << "WARNING: No priors imposed! This might lead to strange results for guessing rate.\n";
	} else if ( PySequence_Check ( pypriors ) ) {
		// Priors are given as a sequence
		for ( i=0; i<Nparams; i++ ) {
			pyblock = PySequence_GetItem ( pypriors, i );
			if ( !strncmp ( PyString_AsString(pyblock), "Uniform", 7 ) ) {
				sscanf ( PyString_AsString(pyblock), "Uniform(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new UniformPrior ( priorpars[0], priorpars[1] ) );
				std::cerr << "Using Uniform Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(pyblock), "Gauss", 5 ) ) {
				sscanf ( PyString_AsString(pyblock), "Gauss(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new GaussPrior ( priorpars[0], priorpars[1] ) );
				std::cerr << "Using Gauss Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(pyblock), "Beta", 4 ) ) {
				sscanf ( PyString_AsString(pyblock), "Beta(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new BetaPrior ( priorpars[0], priorpars[1] ) );
				std::cerr << "Using Beta Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(pyblock), "Gamma", 6 ) ) {
				sscanf ( PyString_AsString(pyblock), "Gamma(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new GammaPrior ( priorpars[0], priorpars[1] ) );
				std::cerr << "Using Gamma Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else {
				std::cerr << "Imposing no constraints on parameter " << i << "\n";
			}
		}
	} else {
		PyErr_Format ( PyExc_ValueError, "priors should be given as dictionary or sequence" );
		return NULL;
	}

	std::vector<double> *start = new std::vector<double> (Nparams);
	if ( pystart!=Py_None ) {
		for ( i=0; i<Nparams; i++ )
			(*start)[i] = PyFloat_AsDouble ( PySequence_GetItem ( pystart, i ) );
	} else
		start = NULL;

	// Perform the real work
	BootstrapList boots = parametricbootstrap ( Nsamples, data, pmf, cuts, start );
	JackKnifeList jack  = jackknifedata       ( data, pmf );

	/************************************************************
	 * Prepare output
	 */
	// Data and samples
	PyArrayObject *pysamples;
	PyArrayObject *pyestimates;
	PyArrayObject *pydeviance;
	int samplesdim[2]   = {Nsamples, Nblocks};
	int estimatesdim[2] = {Nsamples, Nparams};
	pysamples   = (PyArrayObject*) PyArray_FromDims ( 2, samplesdim, PyArray_INT );
	pyestimates = (PyArrayObject*) PyArray_FromDims ( 2, estimatesdim, PyArray_DOUBLE );
	pydeviance = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	for ( i=0; i<Nsamples; i++ ) {
		k = boots.getData ( i );
		for ( j=0; j<Nblocks; j++ ) {
			((int*)pysamples->data)[i*Nblocks+j] = k[j];
		}
		for ( j=0; j<Nparams; j++ ) {
			((double*)pyestimates->data)[i*Nparams+j] = boots.getEst ( i, j );
		}
		((double*)pydeviance->data)[i] = boots.getdeviance ( i );
	}

	/************************************************************
	 * Return
	 */
	pyblock = Py_BuildValue ( "(OOO)", pysamples, pyestimates, pydeviance );
	Py_DECREF ( pysamples );
	Py_DECREF ( pyestimates );
	Py_DECREF ( pydeviance );

	delete data;
	delete pmf;    // also deletes core and sigmoid

	return pyblock;
}

static PyMethodDef psipy_methods[] = {
	{"psibootstrap", (PyCFunction) psibootstrap, METH_VARARGS | METH_KEYWORDS, psibootstrap_doc },
	{NULL,NULL}
};

extern "C" {
	void init_psipy() {
		(void) Py_InitModule("_psipy",psipy_methods);
		import_array();
	}
}
