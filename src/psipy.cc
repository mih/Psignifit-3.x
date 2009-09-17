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
"               mw%g     midpoint and width\n"
"               linear   a+bx\n"
"               log      a+b log(x)\n"
"  priors    constraints on the likelihood estimation. These are expressed in the form of a list of\n"
"            prior names. Valid prior are currently\n"
"               Uniform(%g,%g)\n"
"               Gauss(%g,%g)\n"
"               Beta(%g,%g)\n"
"               Gamma(%g,%g)\n"
"            if an invalid prior is selected, no constraints are imposed at all.\n"
"  cuts      a single number or a sequence of 'cuts' indicating the performances that should be\n"
"            considered 'threshold' performances. This means that in a 2AFC task, cuts==0.5 the\n"
"            'threshold' is somewhere around 75%% correct performance, depending on the lapse rate\n"
"\n"
":Output:\n"
"  samples,estimates,deviance,threshold,bias,acceleration,Rkd,Rpd,outliers,influential\n"
"\n"
"  samples   a nsamplesXnblocks array of the bootstrap sampled data\n"
"  estimates a nsamplesXnparameters array of estimated parameters associated with the data sets\n"
"  deviance  a nsamples array of the associated deviances\n"
"  threshold a nsamples array of thresholds\n"
"  bias      a ncuts array of the bias term associated with the threshold\n"
"  acc       a ncuts array of the acceleration constant associated with the threshold\n"
"  Rkd       a nsamples array of correlations between block index and deviance residuals\n"
"  Rpd       a nsamples array of correlations between model prediction and deviance residuals\n"
"  outliers  a nblocks array indicating points that are outliers\n"
"  influential a nblocks array indicating points that are influential observations\n"
"\n"
"\n"
":Example:\n"
">>> x = [float(2*k) for k in xrange(6)]\n"
">>> k = [34,32,40,48,50,48]\n"
">>> n = [50]*6\n"
">>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]\n"
">>> priors = ('flat','flat','Uniform(0,0.1)')\n"
">>> samples,est,D,thres,bias,acc,Rkd,Rpd,out,influ = bootstrap(d,nsamples=2000,priors=priors)\n"
">>> mean(est[:,0])\n"
"2.7762481672120902\n"
">>> mean(est[:,1])\n"
"1.4243919674602623\n"
"\n";

static char psimcmc_doc [] =
"mcmc ( data, start=None, nsamples=10000, nafc=2, sigmoid='logistic', core='ab', priors=None, stepwidths=None )\n"
"\n"
"Markov Chain Monte Carlo sampling for a psychometric function.\n"
"\n"
":Parameters:\n"
"\n"
"  data      A list of lists or an array of data. The first column should be stimulus intensity,\n"
"            the second column should be number of correct responses in (nAFC) or number of yes-\n"
"            responses (in Yes/No), the third column should be number of trials\n"
"  start     starting values for the markov chain. If this is None, the MAP estimate will be used\n"
"  nsamples  number of samples to be taken from the posterior (note that due to suboptimal sampling,\n"
"            this number may be much lower than the effective number of samples.\n"
"  nafc      number of responses alternatives for nAFC tasks. If nafc==1 a Yes/No task is assumed\n"
"  sigmoid   name of the sigmoid to be fitted. Valid sigmoids are currently\n"
"               logistic\n"
"               gauss\n"
"               gumbel_l\n"
"               gumbel_r\n"
"  core      \"core\"-type of the psychometric function. Valid choices are currently\n"
"               ab     (x-a)/b\n"
"               mw%g   midpoint and width\n"
"               linear a+bx\n"
"               log    a+b log(x)\n"
"  priors    prior distributions on the parameters of the psychometric function. These are expressed\n"
"            in the form of a list of prior names. Valid prior names are currently\n"
"               Uniform(%g,%g)\n"
"               Gauss(%g,%g)\n"
"               Beta(%g,%g)\n"
"               Gamma(%g,%g)\n"
"            if an invalid prior is selected, no constraints are imposed on that parameter resulting\n"
"            in an improper prior distribution.\n"
"  stepwidths standard deviations of the proposal distribution. The best choice is sometimes a bit\n"
"            tricky here. However, as a rule of thumb we can state: if the stepwidths are too small,\n"
"            the samples might not cover the whole posterior, if the stepwidths are too large, most\n"
"            steps will leave the area of high posterior density and will therefore be rejected.\n"
"            Thus, in general stepwidths should be somewhere in the middle.\n"
"\n"
":Output:\n"
"  estimates,deviances\n"
"\n"
"  estimates a nsamplesXnparameters array of parameters sampled from the posterior density\n"
"            of parameters given the data\n"
"  deviances a nsamples array of associated deviances\n"
"\n"
"\n"
":Example:\n"
">>> x = [float(2*k) for k in xrange(6)]\n"
">>> k = [34,32,40,48,50,48]\n"
">>> n = [50]*6\n"
">>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]\n"
">>> priors = ('Gauss(0,1000)','Gauss(0,1000)','Beta(3,100)')\n"
">>> stepwidths = (1.,1.,0.01)\n"
">>> estimates,deviances = mcmc(d,nsamples=10000),priors=priors,stepwidths=stepwidths)\n"
">>> mean(estimates[:,0])\n"
"2.4881405291765764\n"
">>> mean(estimates[:,1])\n"
"1.6920641469955848\n"
"\n";

static char psimapestimate_doc [] =
"mapestimate ( data, nafc=2, sigmoid='logistic', core='ab', priors=None )\n"
"\n"
"MAP or constrained maximum likelihood estimation for a psychometric function\n"
"\n"
":Parameters:\n"
"  data       an array or a list with three columns, the first of which contains stimulus intensity, the second\n"
"             contains the number of correct responses (in nAFC) or the number of 'Yes' responses (in Yes/No),\n"
"             the third column contains the number of trials.\n"
"  nafc       number of response alternatives in an nAFC task. For Yes/No tasks, nafc should be 1\n"
"  sigmoid    type of the sigmoid to be fitted. Valid choices are:\n"
"                'logistic'\n"
"                'gauss'\n"
"                'gumbel_l'\n"
"                'gumbel_r'\n"
"  core       term inside the sigmoid. Valid choices are:\n"
"                'ab'       (x-a)/b\n"
"                'mw%g'     midpoint and width\n"
"                'linear'   a+bx\n"
"                'log'      a+b log(x)\n"
"  priors     constraints to be imposed on the fitting procedure. This should be a list of strings indicating\n"
"             the imposed constraints. Valid choices are:\n"
"                'Uniform(%g,%g)'     Uniform prior on interval\n"
"                'Gauss(%g,%g)'       Gaussian with mean and standard deviation\n"
"                'Beta(%g,%g)'        Beta distribution\n"
"                'Gamma(%g,%g)'       Gamma distribution\n"
"             If an invalid prior is selected, no constraints are imposed on that parameter, resulting in\n"
"             an improper prior distribution.\n"
"\n"
":Output:\n"
"  estimate,deviance\n"
"  estimate  a nparameters array of estimated parameters\n"
"  deviance  the deviance associated with the estimated parameters\n"
"\n"
":Example:\n"
">>> x = [float(2*k) for k in xrange(6)]\n"
">>> k = [34,32,40,48,50,48]\n"
">>> n = [50]*6\n"
">>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]\n"
">>> priors = ('flat','flat','Uniform(0,0.1)')\n"
">>> estimate,deviance = mapestimate ( d, priors=priors )\n"
">>> estimate\n"
"array([2.751768597693296, 1.4572372412562276, 0.015556356934318862], 'd')\n"
">>> deviance\n"
"8.071331367479198\n";


static char psidiagnostics_doc [] =
"diagnostics ( data, params, nafc=2, sigmoid='logistics', core='ab' )\n"
"\n"
"Some diagnostic statistics for a psychometric function fit\n"
"\n"
":Parameters:\n"
"  data     a list of lists or an array containing stimulus intensity in the first column, number of\n"
"           correct responses (for nAFC) or number of Yes-responses (for Yes/No) in the second column,\n"
"           and total number of trials in the third column.\n"
"           Alternatively, only the intensities can be given. In that case, we only obtain the predicted\n"
"           values\n"
"  nafc     number of response alternatives in nAFC tasks. For Yes/No tasks, this should be 1.\n"
"  sigmoid  type of the sigmoid that was used for fitting. Valid choices are:\n"
"             'logistic'\n"
"             'gauss'\n"
"             'gumbel_l'\n"
"             'gumbel_r'\n"
"  core     term inside the sigmoid, i.e. in most cases the parameterization of the sigmoid. Valid\n"
"           choices are:\n"
"              'ab'      (x-a)/b\n"
"              'mw%g'    midpoint and width\n"
"              'linear'  a+bx\n"
"              'log'     a+b log(x)\n"
"\n"
":Output:\n"
"  predicted(,devianceresiduals,deviance,Rpd,Rkd)\n"
"  predicted          predicted values associated with the respective stimulus intensities\n"
"  devianceresiduals  deviance residuals of the data\n"
"  deviance           deviance of the data\n"
"  Rpd                correlation between predicted performance and deviance residuals\n"
"  Rkd                correlation between block index and deviance residuals\n"
"\n"
":Example:\n"
">>> x = [float(2*k) for k in xrange(6)]\n"
">>> k = [34,32,40,48,50,48]\n"
">>> n = [50]*6\n"
">>> d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]\n"
">>> prm = [2.75, 1.45, 0.015]\n"
">>> pred,di,D,Rpd,Rkd = _psipy.diagnostics(d,prm)\n"
">>> D\n"
"8.0748485860836254\n"
">>> di[0]\n"
"1.6893279652591433\n"
">>> Rpd\n"
"-0.19344675783032761\n";

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
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|OiissOO",
				kwlist,
				&pydata,&pystart,&Nsamples,&Nafc,&sigmoidname,&corename,&pypriors,&pycuts ) )
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
	BootstrapList boots = parametricbootstrap ( Nsamples, data, pmf, *cuts, start );
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
	int samplesdim[2]   = {Nsamples, Nblocks};
	int estimatesdim[2] = {Nsamples, Nparams};
	int thresdim[2]     = {Nsamples, Ncuts};
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
		((bool*)pyoutliers->data)[i]   = jack.outlier ( i );
		((bool*)pyinfluential->data)[i]= jack.influential ( i, *ci_lower, *ci_upper );
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
	PsiMClist post ( S->sample(Nsamples) );
	delete S;

	// Data and samples
	PyArrayObject *pyestimates;
	PyArrayObject *pydeviance;
	int estimatesdim[2] = {Nsamples, Nparams};
	pyestimates = (PyArrayObject*) PyArray_FromDims ( 2, estimatesdim, PyArray_DOUBLE );
	pydeviance  = (PyArrayObject*) PyArray_FromDims ( 1, &Nsamples, PyArray_DOUBLE );
	for ( i=0; i<Nsamples; i++ ) {
		for ( j=0; j<Nparams; j++ ) {
			((double*)pyestimates->data)[i*Nparams+j] = post.getEst ( i, j );
		}
		((double*)pydeviance->data)[i] = post.getdeviance ( i );
	}

	pynumber = Py_BuildValue ( "(OO)", pyestimates, pydeviance );

	Py_DECREF ( pyestimates );
	Py_DECREF ( pydeviance );

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
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "O|issOO",
				kwlist,
				&pydata,&Nafc,&sigmoidname,&corename,&pypriors,&pycuts ) )
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
	PsiOptimizer * opt = new PsiOptimizer ( pmf, data );
	*estimate = opt->optimize ( pmf, data );
	delete opt;

	PyArrayObject *pyestimate;
	PyArrayObject *pythres;
	pyestimate = (PyArrayObject*) PyArray_FromDims ( 1, &Nparams, PyArray_DOUBLE );
	pythres    = (PyArrayObject*) PyArray_FromDims ( 1, &Ncuts, PyArray_DOUBLE );
	for (i=0; i<Nparams; i++)
		((double*)pyestimate->data)[i] = (*estimate)[i];

	for (i=0; i<Ncuts; i++)
		((double*)pythres->data)[i] = pmf->getThres ( *estimate, (*cuts)[i] );

	pyout = Py_BuildValue ( "OOd", pyestimate, pythres, pmf->deviance ( *estimate, data ) );

	delete estimate;
	delete data;
	delete pmf;
	delete cuts;
	Py_DECREF ( pyestimate );
	Py_DECREF ( pythres );

	return pyout;
}

static PyObject * psidiagnostics ( PyObject * self, PyObject * args, PyObject * kwargs ) {
	int Nafc ( 2 );                    // Number of response alternatives
	PyObject *pyparams;                // estimated parameters
	PyObject *pydata;                  // python object holding the data
	char *sigmoidname = "logistic";    // name of the sigmoid
	char *corename    = "ab";          // name of the parameterization

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
		NULL };
	if ( !PyArg_ParseTupleAndKeywords ( args, kwargs, "OO|iss",
				kwlist,
				&pydata,&pyparams,&Nafc,&sigmoidname,&corename ) )
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

	params = getparams ( pyparams, Nparams );
	if ( intensityonly==-1)
		devianceresiduals = new std::vector<double> (pmf->getDevianceResiduals ( *params, data ) );

	PyArrayObject *pydevianceresiduals;
	PyArrayObject *pypredicted;
	if ( intensityonly==-1 )
		pydevianceresiduals = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks, PyArray_DOUBLE );
	pypredicted = (PyArrayObject*) PyArray_FromDims ( 1, &Nblocks, PyArray_DOUBLE );
	for (i=0; i<Nblocks; i++) {
		if ( intensityonly==-1 )
			((double*)pydevianceresiduals->data)[i] = (*devianceresiduals)[i];
		((double*)pypredicted->data)[i] = pmf->evaluate ( data->getIntensity ( i ), *params );
	}

	if ( intensityonly==-1 )
		pyout = Py_BuildValue ( "(OOddd)", pypredicted, pydevianceresiduals, pmf->deviance ( *params, data ),
				pmf->getRpd ( *devianceresiduals, *params, data ),
				pmf->getRkd ( *devianceresiduals ) );
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
