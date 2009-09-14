#include <Python.h>
#include <Numeric/arrayobject.h>
#include <vector>
#include <cstdio>

#include "psipp.h"
#include "pytools.h"

static char psibootstrap_doc [] =
"bootstrap ( data, start=None, nsamples=2000, nafc=2, sigmoid=\"logistic\", core=\"ab\", priors=None )\n"
"parametric bootstrap of a psychometric function\n"
"\n"
":Parameters:\n"
"\n"
"  data      A list of lists or an array of data. The first column should be stimulus intensity,\n"
"            the second column should be number of correct responses (in 2AFC) or number of yes-\n"
"            responses (in Yes/No), the third column should be number of trials\n"
"  start     generating values for the bootstrap samples. If this is None, the generating value\n"
"            will be the MAP estimate\n"
"  nsamples  number of bootstrap samples to be drawn\n"
"  nafc      number of alternatives for nAFC tasks. If nafc==1 a Yes/No task is assumed\n"
"  sigmoid   name of the sigmoid to be fitted. Valid sigmoids are currently\n"
"               logistic\n"
"               gauss\n"
"               gumbel_l\n"
"               gumbel_r\n"
"  core      \"core\"-type of the psychometric function. Valid choices are currently\n"
"               ab       (x-a)/b\n"
"               mw(%g)   midpoint + width\n"
"               linear   a+bx\n"
"               log      a+b log(x)\n"
"  priors    constraints on the likelihood estimation. These are expressed in the form of a list of\n"
"            prior names. Valid prior are currently\n"
"               Uniform(%g,%g)\n"
"               Gauss(%g,%g)\n"
"               Beta(%g,%g)\n"
"               Gamma(%g,%g)\n"
"            if an invalid prior is selected, no constraints are imposed at all.\n"
"\n"
":Output:\n"
"  samples,estimates,deviance\n"
"\n"
"  samples   a nsamplesXnblocks array of the bootstrap sampled data\n"
"  estimates a nsamplesXnparameters array of estimated parameters associated with the data sets\n"
"  deviance  a nsamples array of the associated deviances\n"
"\n"
"\n"
":Example:\n"
">>> x = [float(2*k) for k in xrange(6)]\n"
">>> k = [34,32,40,48,50,48]\n"
">>> n = [50]*6\n"
">>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]\n"
">>> priors = ('flat','flat','Uniform(0,0.1)')\n"
">>> samples,estimates,deviance = bootstrap(d,nsamples=2000,priors=priors)\n"
">>> mean(estimates[:,0])\n"
"2.7762481672120902\n"
">>> mean(estimates[:,1])\n"
"1.4243919674602623\n"
"\n";

static PyObject * psibootstrap ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int i,j;
	int Nsamples ( 2000 );             // Number of bootstrap samples
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pydata;                  // python object holding the data
	PyObject *pystart (Py_None);       // python object holding the starting values
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization
	PyObject *pypriors (Py_None);      // prior specs

	int Nblocks;
	PsiSigmoid * sigmoid;
	PsiCore * core;
	PsiData *data;
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
	 * prepare data and sigmoid
	 */
	try {
		data = create_dataset ( pydata, Nafc, &Nblocks );       // prepare data
		sigmoid = getsigmoid ( sigmoidname );                   // prepare sigmoid
		core = getcore ( corename, sigmoid->getcode(), data );  // prepare core object
	} catch (std::string message) {
		PyErr_Format ( PyExc_ValueError, message.c_str() );
		return NULL;
	}

	PsiPsychometric * pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	int Nparams = pmf->getNparams ();
	std::vector<double> cuts (1, .5);

	try {
		setpriors ( pypriors, pmf );
	} catch ( std::string msg ) {
		PyErr_Format ( PyExc_ValueError, msg.c_str() );
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
	std::vector<int> k (Nblocks);
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
	{"bootstrap", (PyCFunction) psibootstrap, METH_VARARGS | METH_KEYWORDS, psibootstrap_doc },
	{NULL,NULL}
};

extern "C" {
	void init_psipy() {
		(void) Py_InitModule("_psipy",psipy_methods);
		import_array();
	}
}
