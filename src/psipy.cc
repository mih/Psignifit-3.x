#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <cstdio>

#include "psipp.h"
#include "pytools.h"
#include "psipy_doc.h"

static PyObject * psibootstrap ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int i,j;
	int Nsamples ( 2000 );             // Number of bootstrap samples
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pydata;                  // python object holding the data
	PyObject *pystart (Py_None);       // python object holding the starting values
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization
	PyObject *pypriors (Py_None);      // prior specs
	PyObject *pycuts (Py_None);        // specify cuts
	int parametric (1);            // perform parametric bootstrap?

	// local variables
	int Nblocks,Nparams;
	PsiSigmoid * sigmoid;
	PsiCore * core;
	PsiData *data;
	PsiPsychometric * pmf;
	PyObject *pynumber;

	/************************************************************
	 * Parse command line
	 */
	static char *kwlist[] = {
		"data",
		"start",
		"nsamples",
		"nafc",
		"sigmoid",
		"core",
		"priors",
		"cuts",
		"parametric",
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|OiissOOi",
				kwlist,
				&pydata,&pystart,&Nsamples,&Nafc,&sigmoidname,&corename,&pypriors,&pycuts, &parametric ) )
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

	pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	Nparams = pmf->getNparams ();

	std::vector<double> *cuts;
	int Ncuts;
	try {
		cuts = getcuts ( pycuts, &Ncuts );
		setpriors ( pypriors, pmf );
	} catch ( std::string msg ) {
		PyErr_Format ( PyExc_ValueError, msg.c_str() );
		return NULL;
	}

	std::vector<double> *start = new std::vector<double> (Nparams);
	if ( pystart!=Py_None ) {
		for ( i=0; i<Nparams; i++ ) {
			pynumber = PySequence_GetItem ( pystart, i );
			(*start)[i] = PyFloat_AsDouble ( pynumber );
			Py_DECREF ( pynumber );
		}
	} else
		start = NULL;

	// Perform the real work
	BootstrapList boots = bootstrap ( Nsamples, data, pmf, *cuts, start, true, bool(parametric) );
	JackKnifeList jack  = jackknifedata       ( data, pmf );

	/************************************************************
	 * Prepare output
	 */
	// Data and samples
	PyArrayObject *pysamples;
	PyArrayObject *pyestimates;
	PyArrayObject *pydeviance;
	PyArrayObject *pythres;
	PyArrayObject *pyRpd;
	PyArrayObject *pyRkd;
	PyArrayObject *pyoutliers;
	PyArrayObject *pyinfluential;
	PyArrayObject *pybias;
	PyArrayObject *pyacc;
	std::vector<int> k (Nblocks);
	npy_intp samplesdim[2]   = {Nsamples, Nblocks};
	npy_intp estimatesdim[2] = {Nsamples, Nparams};
	npy_intp thresdim[2]     = {Nsamples, Ncuts};
	npy_intp nNsamples (Nsamples);
	npy_intp nNblocks (Nblocks);
	npy_intp nNcuts (Ncuts);
	/*
	pysamples   = (PyArrayObject*) PyArray_FromDims ( 2, samplesdim, PyArray_INT );
	pyestimates = (PyArrayObject*) PyArray_FromDims ( 2, estimatesdim, PyArray_DOUBLE );
	pydeviance  = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	pythres     = (PyArrayObject*) PyArray_FromDims ( 2, thresdim, PyArray_DOUBLE );
	pyRpd       = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	pyRkd       = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	pyoutliers  = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks,  PyArray_INT );
	pyinfluential = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks,  PyArray_INT );
	pybias      = (PyArrayObject*) PyArray_FromDims ( 1, &Ncuts, PyArray_DOUBLE );
	pyacc       = (PyArrayObject*) PyArray_FromDims ( 1, &Ncuts, PyArray_DOUBLE );
	*/
	pysamples     = (PyArrayObject*) PyArray_SimpleNew ( 2, samplesdim, NPY_INT );
	pyestimates   = (PyArrayObject*) PyArray_SimpleNew ( 2, estimatesdim, NPY_DOUBLE );
	pydeviance    = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pythres       = (PyArrayObject*) PyArray_SimpleNew ( 2, thresdim, NPY_DOUBLE );
	pyRpd         = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pyRkd         = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pyoutliers    = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNblocks,  NPY_BOOL );
	pyinfluential = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNblocks,  NPY_DOUBLE );
	pybias        = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNcuts, NPY_DOUBLE );
	pyacc         = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNcuts, NPY_DOUBLE );
	for ( i=0; i<Nsamples; i++ ) {
		k = boots.getData ( i );
		for ( j=0; j<Nblocks; j++ ) {
			((int*)pysamples->data)[i*Nblocks+j] = k[j];
		}
		for ( j=0; j<Nparams; j++ ) {
			((double*)pyestimates->data)[i*Nparams+j] = boots.getEst ( i, j );
		}
		((double*)pydeviance->data)[i] = boots.getdeviance ( i );
		for ( j=0; j<Ncuts; j++ )
			((double*)pythres->data)[i*Ncuts+j]    = boots.getThres_byPos ( i, j );
		((double*)pyRpd->data)[i]      = boots.getRpd(i);
		((double*)pyRkd->data)[i]      = boots.getRkd(i);
	}
	std::vector<double> *ci_lower = new std::vector<double> ( Nparams );
	std::vector<double> *ci_upper = new std::vector<double> ( Nparams );
	for ( i=0; i<Nparams; i++ ) {
		(*ci_lower)[i] = boots.getPercentile(0.025,i);
		(*ci_upper)[i] = boots.getPercentile(0.975,i);
	}
	for ( i=0; i<Nblocks; i++ ) {
		((bool*)pyoutliers->data)[i]   = bool(jack.outlier ( i ));
		((double*)pyinfluential->data)[i]= double(jack.influential ( i, *ci_lower, *ci_upper ));
	}
	delete ci_lower;
	delete ci_upper;

	// BCa
	for ( i=0; i<Ncuts; i++ ) {
		((double*)pybias->data)[i] = boots.getBias(i);
		((double*)pyacc->data)[i]  = boots.getAcc(i);
	}

	/************************************************************
	 * Return
	 */
	pynumber = Py_BuildValue ( "(OOOOOOOOOO)", pysamples, pyestimates, pydeviance, pythres, pybias, pyacc, pyRpd, pyRkd, pyoutliers, pyinfluential );
	Py_DECREF ( pysamples );
	Py_DECREF ( pyestimates );
	Py_DECREF ( pydeviance );
	Py_DECREF ( pythres );
	Py_DECREF ( pyRpd );
	Py_DECREF ( pyRkd );
	Py_DECREF ( pyoutliers );
	Py_DECREF ( pyinfluential );

	delete data;
	delete pmf;    // also deletes core and sigmoid
	delete start;

	return pynumber;
}

static PyObject * psimcmc ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int Nsamples ( 10000 );             // Number of bootstrap samples
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pydata;                  // python object holding the data
	PyObject *pystart (Py_None);       // python object holding the starting values
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "mw0.1";          // name of the parameterization
	PyObject *pypriors (Py_None);      // prior specs
	PyObject *pysteps (Py_None);       // stepwidths of the proposal distribution

	PyObject * pynumber;
	int i,j, Nparams, Nblocks;
	PsiData * data;
	PsiPsychometric * pmf;
	PsiCore * core;
	PsiSigmoid * sigmoid;

	static char *kwlist[] = {
		"data",
		"start",
		"nsamples",
		"nafc",
		"sigmoid",
		"core",
		"priors",
		"stepwidths",
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|OiissOO",
				kwlist,
				&pydata,&pystart,&Nsamples,&Nafc,&sigmoidname,&corename,&pypriors,&pysteps ) )
		return NULL;

	try {
		data = create_dataset ( pydata, Nafc, &Nblocks );       // prepare data
		sigmoid = getsigmoid ( sigmoidname );                   // prepare sigmoid
		core = getcore ( corename, sigmoid->getcode(), data );  // prepare core object
	} catch (std::string message) {
		PyErr_Format ( PyExc_ValueError, message.c_str() );
		return NULL;
	}

	pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	Nparams = pmf->getNparams ();

	try {
		setpriors ( pypriors, pmf );
	} catch ( std::string msg ) {
		PyErr_Format ( PyExc_ValueError, msg.c_str() );
		return NULL;
	}

	std::vector<double> *start = new std::vector<double> (Nparams);
	if ( pystart!=Py_None ) {
		for ( i=0; i<Nparams; i++ ) {
			pynumber = PySequence_GetItem ( pystart, i );
			(*start)[i] = PyFloat_AsDouble ( pynumber );
			Py_DECREF ( pynumber );
		}
	} else {
		// If we have no explicit starting value, we start with the MAP estimate
		std::cerr << "No starting value for chain specified -- using MAP estimate\n";
		PsiOptimizer * opt = new PsiOptimizer ( pmf, data );
		*start = opt->optimize ( pmf, data );
		delete opt;
	}

	// Perform sampling
	MetropolisHastings * S = new MetropolisHastings ( pmf, data, new GaussRandom () );
	S->setTheta ( *start );
	try {
		setstepwidths ( pysteps, S );
	} catch ( std::string msg ) {
		PyErr_Format ( PyExc_ValueError, msg.c_str() );
		return NULL;
	}
	MCMCList post ( S->sample(Nsamples) );
	delete S;

	// Data and samples
	PyArrayObject *pyestimates;
	PyArrayObject *pydeviance;
	PyArrayObject *pyposterior_predictive_data;
	PyArrayObject *pyposterior_predictive_deviances;
	PyArrayObject *pyposterior_predictive_Rpd;
	PyArrayObject *pyposterior_predictive_Rkd;
	PyArrayObject *pylogposterior_ratios;
	npy_intp estimatesdim[2] = {Nsamples, Nparams};
	npy_intp datadim[2] = {Nsamples, Nblocks};
	npy_intp nNsamples ( Nsamples );
	/*
	pyestimates = (PyArrayObject*) PyArray_FromDims ( 2, estimatesdim, PyArray_DOUBLE );
	pydeviance  = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	*/
	pyestimates                      = (PyArrayObject*) PyArray_SimpleNew ( 2, estimatesdim, NPY_DOUBLE );
	pydeviance                       = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pyposterior_predictive_data      = (PyArrayObject*) PyArray_SimpleNew ( 2, datadim, NPY_INT );
	pyposterior_predictive_deviances = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pyposterior_predictive_Rpd       = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pyposterior_predictive_Rkd       = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNsamples, NPY_DOUBLE );
	pylogposterior_ratios            = (PyArrayObject*) PyArray_SimpleNew ( 2, datadim, NPY_DOUBLE );
	for ( i=0; i<Nsamples; i++ ) {
		for ( j=0; j<Nparams; j++ ) {
			((double*)pyestimates->data)[i*Nparams+j] = post.getEst ( i, j );
		}
		((double*)pydeviance->data)[i] = post.getdeviance ( i );
		for ( j=0; j<Nblocks; j++ ) {
			((int*)pyposterior_predictive_data->data)[i*Nblocks+j] = post.getppData ( i, j );
			((double*)pylogposterior_ratios->data)[i*Nblocks+j]    = post.getlogratio ( i, j );
		}
		((double*)pyposterior_predictive_deviances->data)[i] = post.getppDeviance ( i );
		((double*)pyposterior_predictive_Rpd->data)[i]       = post.getppRpd ( i );
		((double*)pyposterior_predictive_Rkd->data)[i]       = post.getppRkd ( i );
	}

	pynumber = Py_BuildValue ( "(OOOOOOO)", pyestimates, pydeviance,
			pyposterior_predictive_data, pyposterior_predictive_deviances, pyposterior_predictive_Rpd, pyposterior_predictive_Rkd,pylogposterior_ratios );

	Py_DECREF ( pyestimates );
	Py_DECREF ( pydeviance );
	Py_DECREF ( pyposterior_predictive_data );
	Py_DECREF ( pyposterior_predictive_deviances );
	Py_DECREF ( pyposterior_predictive_Rpd );
	Py_DECREF ( pyposterior_predictive_Rkd );
	Py_DECREF ( pylogposterior_ratios );

	delete pmf;
	delete data;
	delete start;

	return pynumber;
}

static PyObject * psimapestimate ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pydata;                  // python object holding the data
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization
	PyObject *pypriors (Py_None);      // prior specs
	PyObject *pycuts (Py_None);        // cuts
	PyObject *pystart (Py_None);       // starting value for the optimizer

	PyObject * pyout;
	int i, Nparams, Nblocks;
	PsiData * data;
	PsiPsychometric * pmf;
	PsiCore * core;
	PsiSigmoid * sigmoid;

	static char *kwlist[] = {
		"data",
		"nafc",
		"sigmoid",
		"core",
		"priors",
		"cuts",
		"start",
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|issOOO",
				kwlist,
				&pydata,&Nafc,&sigmoidname,&corename,&pypriors,&pycuts,&pystart ) )
		return NULL;

	try {
		data = create_dataset ( pydata, Nafc, &Nblocks );       // prepare data
		sigmoid = getsigmoid ( sigmoidname );                   // prepare sigmoid
		core = getcore ( corename, sigmoid->getcode(), data );  // prepare core object
	} catch (std::string message) {
		PyErr_Format ( PyExc_ValueError, message.c_str() );
		return NULL;
	}

	pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	Nparams = pmf->getNparams ();

	int Ncuts;
	std::vector<double> *cuts;
	try {
		cuts = getcuts ( pycuts, &Ncuts );
		setpriors ( pypriors, pmf );
	} catch ( std::string msg ) {
		PyErr_Format ( PyExc_ValueError, msg.c_str() );
		return NULL;
	}

	std::vector<double> *estimate = new std::vector<double> (Nparams);
	std::vector<double> *start;
	PsiOptimizer * opt = new PsiOptimizer ( pmf, data );
	if ( pystart == Py_None ) {
		*estimate = opt->optimize ( pmf, data );
	} else {
		start = getcuts ( pystart, &Nparams );
		*estimate = opt->optimize ( pmf, data, start );
		delete start;
	}
	delete opt;

	Matrix * H = pmf->ddnegllikeli ( *estimate, data );
	Matrix * I = H->inverse();

	PyArrayObject *pyestimate;
	PyArrayObject *pythres;
	PyArrayObject *pyvars;
	npy_intp nNparams (Nparams);
	npy_intp nNcuts (Ncuts);
	/*
	pyestimate = (PyArrayObject*) PyArray_FromDims ( 1, &Nparams, PyArray_DOUBLE );
	pythres    = (PyArrayObject*) PyArray_FromDims ( 1, &Ncuts, PyArray_DOUBLE );
	*/
	pyestimate = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNparams, NPY_DOUBLE );
	pythres    = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNcuts, NPY_DOUBLE );
	pyvars     = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNparams, NPY_DOUBLE );
	for (i=0; i<Nparams; i++) {
		((double*)pyestimate->data)[i] = (*estimate)[i];
		((double*)pyvars->data)[i] = (*I)(i,i);
	}

	for (i=0; i<Ncuts; i++)
		((double*)pythres->data)[i] = pmf->getThres ( *estimate, (*cuts)[i] );

	pyout = Py_BuildValue ( "OOOd", pyestimate, pyvars, pythres, pmf->deviance ( *estimate, data ) );

	delete estimate;
	delete data;
	delete pmf;
	delete cuts;
	delete H;
	delete I;
	Py_DECREF ( pyestimate );
	Py_DECREF ( pythres );
	Py_DECREF ( pyvars );

	return pyout;
}

static PyObject * psidiagnostics ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pyparams;                // estimated parameters
	PyObject *pydata;                  // python object holding the data
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization
	PyObject *pycuts (Py_None);         // cuts

	PyObject * pyout;
	int intensityonly(-1);    // This variable is a bit tricky: if it remains -1, we have "real" data, if it becomes 1, we only have intensities and can safe some steps
	int i, Nparams, Nblocks;
	PsiData * data;
	PsiPsychometric * pmf;
	PsiCore * core;
	PsiSigmoid * sigmoid;
	std::vector<double> *params;
	std::vector<double> *devianceresiduals;

	static char *kwlist[] = {
		"data",
		"params",
		"nafc",
		"sigmoid",
		"core",
		"cuts",
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "OO|issO",
				kwlist,
				&pydata,&pyparams,&Nafc,&sigmoidname,&corename,&pycuts ) )
		return NULL;

	try {
		data = create_dataset ( pydata, Nafc, &Nblocks, &intensityonly );   // prepare data
		sigmoid = getsigmoid ( sigmoidname );                   // prepare sigmoid
		core = getcore ( corename, sigmoid->getcode(), data );  // prepare core object
	} catch (std::string message) {
		PyErr_Format ( PyExc_ValueError, message.c_str() );
		return NULL;
	}

	pmf = new PsiPsychometric ( Nafc, core, sigmoid );
	Nparams = pmf->getNparams ();

	std::vector<double> *cuts = NULL;
	int Ncuts (0);
	try {
		cuts = getcuts ( pycuts, &Ncuts );
		params = getparams ( pyparams, Nparams );
	} catch (std::string message) {
		PyErr_Format ( PyExc_ValueError, message.c_str() );
		return NULL;
	}

	if ( intensityonly==-1)
		devianceresiduals = new std::vector<double> (pmf->getDevianceResiduals ( *params, data ) );

	PyArrayObject *pydevianceresiduals;
	PyArrayObject *pypredicted;
	PyArrayObject *pythres;
	npy_intp nNblocks (Nblocks);
	npy_intp nNcuts (Ncuts);
	if ( intensityonly==-1 )
		// pydevianceresiduals = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks, PyArray_DOUBLE );
		pydevianceresiduals = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNblocks, NPY_DOUBLE );
	/*
	pypredicted = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks, PyArray_DOUBLE );
	pythres     = (PyArrayObject*) PyArray_FromDims ( 1, &Ncuts,   PyArray_DOUBLE );
	*/
	pypredicted = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNblocks, NPY_DOUBLE );
	pythres     = (PyArrayObject*) PyArray_SimpleNew ( 1, &nNcuts,   NPY_DOUBLE );
	for (i=0; i<Nblocks; i++) {
		if ( intensityonly==-1 )
			((double*)pydevianceresiduals->data)[i] = (*devianceresiduals)[i];
		((double*)pypredicted->data)[i] = pmf->evaluate ( data->getIntensity ( i ), *params );
	}
	for (i=0; i<Ncuts; i++) {
		((double*)pythres->data)[i] = pmf->getThres ( *params, (*cuts)[i] );
	}

	if ( intensityonly==-1 )
		pyout = Py_BuildValue ( "(OOdOdd)", pypredicted, pydevianceresiduals, pmf->deviance ( *params, data ), pythres,
				pmf->getRpd ( *devianceresiduals, *params, data ),
				pmf->getRkd ( *devianceresiduals, data ) );
	else
		pyout = Py_BuildValue ( "O", pypredicted );

	delete data;
	delete pmf;
	delete params;
	if ( intensityonly==-1 ) {
		delete devianceresiduals;
		Py_DECREF ( pydevianceresiduals );
	}
	Py_DECREF ( pypredicted );
	Py_DECREF ( pythres );

	return pyout;
}

static PyMethodDef psipy_methods[] = {
	{"bootstrap", (PyCFunction) psibootstrap, METH_VARARGS | METH_KEYWORDS, psibootstrap_doc },
	{"mcmc", (PyCFunction) psimcmc, METH_VARARGS | METH_KEYWORDS, psimcmc_doc },
	{"mapestimate", (PyCFunction) psimapestimate, METH_VARARGS | METH_KEYWORDS, psimapestimate_doc },
	{"diagnostics", (PyCFunction) psidiagnostics, METH_VARARGS | METH_KEYWORDS, psidiagnostics_doc },
	{NULL,NULL}
};

extern "C" {
	void init_psipy() {
		(void) Py_InitModule("_psipy",psipy_methods);
		import_array();
	}
}
