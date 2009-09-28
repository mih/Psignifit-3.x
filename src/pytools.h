#ifndef PYTOOLS_H
#define PYTOOLS_H

#include <vector>
#include <cstdlib>
#include <string>

PsiData * create_dataset ( PyObject * pydata, int Nafc, int *nblocks, int *allow1d=NULL ) {
	if ( !PySequence_Check ( pydata ) )
		throw std::string ( "data should be a sequence" );

	PyObject * pyblock, * pynumber;

	int Nblocks ( PySequence_Size ( pydata ) );
	*nblocks = Nblocks;
	std::vector<double> x ( Nblocks );
	std::vector<int> k ( Nblocks );
	std::vector<int> n ( Nblocks );

	for ( int i=0; i<Nblocks; i++ ) {
		pyblock = PySequence_GetItem ( pydata, i );
		if ( PySequence_Check ( pyblock ) && (allow1d==NULL || *allow1d==-1) ) {
			// The data are a sequence, we are interested in "real" data
			if ( PySequence_Size ( pyblock ) != 3 ) {
				char msg[50];
				sprintf ( msg,"data in block %d do not have 3 entries", i );
				Py_DECREF ( pyblock );
				throw std::string ( msg );
			}
			pynumber = PySequence_GetItem ( pyblock, 0 );   x[i] = PyFloat_AsDouble ( pynumber );    Py_DECREF ( pynumber );
			pynumber = PySequence_GetItem ( pyblock, 1 );   k[i] = PyInt_AsLong ( pynumber );        Py_DECREF ( pynumber );
			pynumber = PySequence_GetItem ( pyblock, 2 );   n[i] = PyInt_AsLong ( pynumber );        Py_DECREF ( pynumber );
			// std::cerr << i << " " << x[i] << " " << k[i] << " " << n[i] << "\n";
		} else if ( PyNumber_Check ( pyblock ) && allow1d!=NULL ) {     // We can only get to the intensity stuff, if allow1d is not NULL
			// The data are just a number, we are interested in the intensity only
			pynumber = PyNumber_Float(pyblock);
			*allow1d = 1;                          // From now on only intensities alone are allowed
			x[i] = PyFloat_AsDouble ( pynumber );
			Py_DECREF ( pynumber );
			k[i] = 0;
			n[i] = 0;
		} else {
			char msg[50];
			sprintf ( msg, "data in block %d are strange, I don't know what to do with them", i );
			Py_DECREF ( pyblock );
			throw std::string ( msg );
		}
		Py_DECREF ( pyblock );
	}

	return new PsiData ( x, n, k, Nafc );
}

PsiSigmoid * getsigmoid ( const char * sigmoidname ) {
	if ( !strcmp(sigmoidname,"logistic") ) {
		// std::cerr << "Using logistic sigmoid\n";
		return new PsiLogistic;
	} else if ( !strcmp(sigmoidname,"gauss") ) {
		// std::cerr << "Using gaussian cdf sigmoid\n";
		return new PsiGauss;
	} else if ( !strcmp(sigmoidname,"gumbel_l") || !strcmp(sigmoidname,"lgumbel") ) {
		// std::cerr << "Using gumbelL sigmoid\n";
		return new PsiGumbelL;
	} else if ( !strcmp(sigmoidname,"gumbel_r") || !strcmp(sigmoidname,"rgumbel") ) {
		// std::cerr << "Using gumbelR sigmoid\n";
		return new PsiGumbelR;
	} else if ( !strcmp(sigmoidname,"cauchy") ) {
		return new PsiCauchy;
	} else {
		throw std::string ( "invalid sigmoid type" );
	}
}

PsiCore * getcore ( const char * corename, int sigmoidcode, const PsiData * data ) {
	if ( !strcmp(corename,"ab") ) {
		// std::cerr << "Using core ab\n";
		return new abCore;
	} else if ( !strncmp(corename,"mw",2) ) {
		double alpha(0.1);
		if ( sscanf ( corename, "mw%lf", &alpha )==0 )
			alpha = 0.1;
		if (alpha<=0 || alpha>=1)
			throw std::string ( "alpha should be >0 and <1, but isn't." );
		return new mwCore ( sigmoidcode, alpha );
	} else if ( !strcmp(corename,"linear") ) {
		// std::cerr << "Using linear core\n";
		return new linearCore;
	} else if ( !strcmp(corename,"log") || !strcmp(corename,"logarithmic") ) {
		// std::cerr << "Using logarithmic core\n";
		return new logCore ( data );
	} else {
		throw std::string ( "invalid core type" );
	}
}

void setpriors ( PyObject * pypriors, PsiPsychometric * pmf ) {
	double priorpars[10];
	int i, Nparams ( pmf->getNparams() );
	PyObject * singleprior;
	if ( pypriors == Py_None ) {
		std::cerr << "WARNING: No priors imposed! This might lead to strange results for guessing rate.\n";
	} else if ( PySequence_Check ( pypriors ) ) {
		// Priors are given as a sequence
		for ( i=0; i<Nparams; i++ ) {
			singleprior = PySequence_GetItem ( pypriors, i );
			if ( !strncmp ( PyString_AsString(singleprior), "Uniform", 7 ) ) {
				sscanf ( PyString_AsString(singleprior), "Uniform(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new UniformPrior ( priorpars[0], priorpars[1] ) );
				// std::cerr << "Using Uniform Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(singleprior), "Gauss", 5 ) ) {
				sscanf ( PyString_AsString(singleprior), "Gauss(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new GaussPrior ( priorpars[0], priorpars[1] ) );
				// std::cerr << "Using Gauss Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(singleprior), "Beta", 4 ) ) {
				sscanf ( PyString_AsString(singleprior), "Beta(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new BetaPrior ( priorpars[0], priorpars[1] ) );
				// std::cerr << "Using Beta Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else if ( !strncmp ( PyString_AsString(singleprior), "Gamma", 6 ) ) {
				sscanf ( PyString_AsString(singleprior), "Gamma(%lf,%lf)", priorpars,priorpars+1 );
				pmf->setPrior ( i, new GammaPrior ( priorpars[0], priorpars[1] ) );
				// std::cerr << "Using Gamma Prior with params " << priorpars[0] << " " << priorpars[1] << " for parameter " << i << "\n";
			} else {
				// std::cerr << "Imposing no constraints on parameter " << i << "\n";
			}
			Py_DECREF ( singleprior );
		}
	} else {
		throw std::string ( "priors should be given as a sequence" );
	}
}

void setstepwidths ( PyObject * pysteps, MetropolisHastings * S ) {
	int i, Nparams ( S->getNparams() );
	PyObject * singlestep;

	if ( pysteps == Py_None ) {
		// std::cerr << "Warning: stepwidths were not touched! This might lead to bad convergence of the markov chains.\n";
	} else if ( PySequence_Check ( pysteps ) ) {
		for ( i=0; i<Nparams; i++ ) {
			singlestep = PySequence_GetItem ( pysteps, i );
			S->setstepsize ( PyFloat_AsDouble ( singlestep ), i );
			Py_DECREF ( singlestep );
		}
	} else {
		throw std::string ( "stepwidths should be given as a sequence" );
	}
}

std::vector<double> *getparams ( PyObject * pyparams, int Nparams ) {
	if ( !PySequence_Check ( pyparams ) )
		throw std::string ( "parameters should be a sequence" );

	if ( PySequence_Size ( pyparams ) != Nparams )
		throw std::string ( "parameters have incorrect lengths" );

	int i;
	std::vector<double> *out = new std::vector<double> ( Nparams );
	PyObject * pynumber;
	for ( i=0; i<Nparams; i++ ) {
		pynumber = PySequence_GetItem ( pyparams, i );
		(*out)[i] = PyFloat_AsDouble ( pynumber );
		Py_DECREF ( pynumber );
	}

	return out;
}

std::vector<double> *getcuts ( PyObject * pycuts, int * Ncuts ) {
	if ( pycuts == Py_None ) {
		*Ncuts = 1;
		return new std::vector<double> (1,.5);
	} else if ( PySequence_Check ( pycuts ) ) {
		int i,len ( PySequence_Size ( pycuts ) );
		*Ncuts = len;
		std::vector<double>* out = new std::vector<double> (len);
		PyObject *pynumber;
		for ( i=0; i<len; i++ ) {
			pynumber = PySequence_GetItem ( pycuts, i );
			(*out)[i] = PyFloat_AsDouble ( pynumber );
			Py_DECREF ( pynumber );
		}
		return out;
	} else if ( PyNumber_Check ( pycuts ) ) {
		*Ncuts = 1;
		return new std::vector<double> (1, PyFloat_AsDouble ( pycuts ) );
	} else {
		throw std::string ( "cuts should be a sequence or a single number" );
	}
}
#endif
