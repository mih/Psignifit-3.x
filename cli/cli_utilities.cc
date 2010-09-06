#include "cli_utilities.h"
#include <cstring>

PsiData * allocateDataFromFile ( std::string fname, int nafc ) {
	std::fstream infile ( fname.c_str() );
	std::list<std::string> inputlines;
	char line[100];
	unsigned int i;
	std::list<std::string>::iterator inputlines_i;
	
	while ( !infile.eof() ) {
		infile.getline ( line, 100 );
		if (strlen(line)>0)
			inputlines.push_back ( std::string ( line ) );
	}

	std::vector<double> x ( inputlines.size() );
	std::vector<int>    k ( inputlines.size() );
	std::vector<int>    n ( inputlines.size() );
	double xx,kk,nn;

	for ( i=0, inputlines_i=inputlines.begin(); inputlines_i != inputlines.end(); inputlines_i++, i++ ) {
		sscanf ( (*inputlines_i).c_str(), "%lf %lf %lf\n", &xx,&kk,&nn );
		if ( kk<1 ) kk *= nn;
		if ( kk != round(kk) ) std::cerr << "Warning: number of correct trials is " << kk << " and not an integer\n";
		x[i] = xx;
		k[i] = kk;
		n[i] = nn;
		if ( nn < kk ) std::cerr << "Warning: number of correct trials is larger than the total number of trials!\n";
	}

	return new PsiData ( x, n, k, nafc );
}

PsiPsychometric * allocatePsychometric ( std::string core, std::string sigmoid, int nafc, PsiData* data, bool verbose=false ) {
	PsiSigmoid *psisigmoid;
	PsiCore *psicore;
	double alpha;
	PsiPsychometric * out;

	if ( verbose ) std::cerr << "setting up psychometric function model (";

	if ( !sigmoid.compare("cauchy") ) {
		if (verbose) std::cerr <<  PsiCauchy::getDescriptor();
		psisigmoid = new PsiCauchy();
	} else if ( !sigmoid.compare("exponential") ) {
		if (verbose) std::cerr <<  PsiExponential::getDescriptor();
		psisigmoid = new PsiExponential();
	} else if ( !sigmoid.compare("gauss") ) {
		if (verbose) std::cerr <<  PsiGauss::getDescriptor();
		psisigmoid = new PsiGauss();
	} else if ( !sigmoid.compare("gumbel_l") ) {
		if (verbose) std::cerr <<  PsiGumbelL::getDescriptor();
		psisigmoid = new PsiGumbelL();
	} else if ( !sigmoid.compare("gumbel_r") ) {
		if (verbose) std::cerr <<  PsiGumbelR::getDescriptor();
		psisigmoid = new PsiGumbelR();
	} else if ( !sigmoid.compare("logistic") ) {
		if (verbose) std::cerr <<  PsiLogistic::getDescriptor();
		psisigmoid = new PsiLogistic();
	} else {
		std::cerr << "Unknown sigmoid " << sigmoid << "\n";
		exit ( -1 );
	}

	if ( verbose ) std::cerr << ",";

	if ( !core.compare("ab") ) {
		if (verbose) std::cerr << abCore::getDescriptor();
		psicore = new abCore ();
	} else if ( !core.compare ( 0, 2, "mw" ) ) {
		alpha = atof ( core.substr( 2, 100 ).c_str() );
		if (verbose) std::cerr << mwCore::getDescriptor() << alpha;
		psicore = new mwCore ( data, psisigmoid->getcode(), alpha );
	} else if ( !core.compare ( "linear" ) ) {
		if (verbose) std::cerr << linearCore::getDescriptor();
		psicore = new linearCore ();
	} else if ( !core.compare ( "log" ) ) {
		if (verbose) std::cerr << logCore::getDescriptor();
		psicore = new logCore ( data );
	} else if ( !core.compare ( "poly" ) ) {
		if (verbose) std::cerr << polyCore::getDescriptor();
		psicore = new polyCore ( data );
	} else if ( !core.compare ( "weibull" ) ) {
		if (verbose) std::cerr << weibullCore::getDescriptor();
		psicore = new weibullCore ( data );
	} else {
		std::cerr << "Unknown core " << sigmoid << "\n";
		exit ( -1 );
	}

	if ( verbose ) std::cerr << ") ";

	out = new PsiPsychometric ( nafc, psicore, psisigmoid );

	delete psisigmoid;
	delete psicore;

	return out;
}

PsiPrior * allocatePrior ( std::string prior ) {
	PsiPrior *psiprior;
	double prm1,prm2;

	if ( !prior.compare ( 0, 4, "Beta" ) ) {
		sscanf ( prior.c_str(), "Beta(%lf,%lf)", &prm1, &prm2 );
		psiprior = new BetaPrior ( prm1, prm2 );
	} else if ( !prior.compare ( 0, 5, "Gamma" ) ) {
		sscanf ( prior.c_str(), "Gamma(%lf,%lf)", &prm1, &prm2 );
		psiprior = new GammaPrior ( prm1, prm2 );
	} else if ( !prior.compare ( 0, 6, "nGamma" ) ) {
		sscanf ( prior.c_str(), "nGamma(%lf,%lf)", &prm1, &prm2 );
		psiprior = new nGammaPrior ( prm1, prm2 );
	} else if ( !prior.compare ( 0, 5, "Gauss" ) ) {
		sscanf ( prior.c_str(), "Gauss(%lf,%lf)", &prm1, &prm2 );
		psiprior = new GaussPrior ( prm1, prm2 );
	} else if ( !prior.compare ( 0, 7, "Uniform" ) ) {
		sscanf ( prior.c_str(), "Uniform(%lf,%lf)", &prm1, &prm2 );
		psiprior = new UniformPrior ( prm1, prm2 );
	} else if ( !prior.compare ( "None" ) ) {
		psiprior = NULL;
	} else {
		std::cerr << "Unknown prior: " << prior << "\n";
		exit ( -1 );
	}

	return psiprior;
}

void setPriors ( PsiPsychometric * pmf, std::string prior1, std::string prior2, std::string prior3, std::string prior4 ) {
	PsiPrior * prior;

	prior = allocatePrior ( prior1 );
	if (prior!=NULL) pmf->setPrior ( 0, prior );
	delete prior;
	prior = allocatePrior ( prior2 );
	if (prior!=NULL) pmf->setPrior ( 1, prior );
	delete prior;
	prior = allocatePrior ( prior3 );
	if (prior!=NULL) pmf->setPrior ( 2, prior );
	delete prior;
	if ( pmf->getNparams()>3 ) {
		prior = allocatePrior ( prior4 );
		if (prior!=NULL) pmf->setPrior ( 3, prior );
		delete prior;
	}
}

std::vector<double> getCuts ( std::string cuts ) {
	size_t pos0(0),pos1(0);
	unsigned int nfound(1);

	pos0 = cuts.find_first_of ( "," );
	while ( pos0!=std::string::npos ) {
		nfound ++;
		pos0 = cuts.find_first_of ( ",", pos0+1 );
	}

	std::vector<double> out ( nfound );
	unsigned int i(0);

	pos0 = 0;
	pos1 = cuts.find_first_of ( "," );
	for (i=0; i<nfound; i++) {
		if (pos0==pos1) {
			pos1 = 1000;
		}
		out[i] = atof ( cuts.substr ( pos0, pos1-pos0 ).c_str() );
		pos0 = pos1+1;
		pos1 = cuts.find_first_of ( ",", pos0 );
	}

	return out;
}

void print ( std::vector<double> theta, bool matlabformat, std::string varname, FILE *ofile ) {
	unsigned int i;
	if ( matlabformat ) {
		fprintf ( ofile, "results.%s = [ %lf", varname.c_str(), theta[0] );
		for ( i=1; i<theta.size(); i++ ) {
			fprintf ( ofile, ", %lf", theta[i] );
		}
		fprintf ( ofile, "];\n" );
	} else {
		fprintf ( ofile, "\n# %s\n", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			fprintf ( ofile, " %lf ", theta[i] );
		}
		fprintf ( ofile, "\n" );
	}
	fprintf ( ofile, "\n" );
}

void print ( std::vector<int> theta, bool matlabformat, std::string varname, FILE *ofile ) {
	unsigned int i;

	if ( matlabformat ) {
		fprintf ( ofile, "results.%s = [ %d", varname.c_str(), theta[0] );
		for ( i=1; i<theta.size(); i++ ) {
			fprintf ( ofile, ", %d", theta[i] );
		}
		fprintf ( ofile, "];\n" );
	} else {
		fprintf ( ofile, "\n# %s\n", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			fprintf ( ofile, " %d ", theta[i] );
		}
		fprintf ( ofile, "\n" );
	}
}

void print ( double theta, bool matlabformat, std::string varname, FILE *ofile ) {
	if ( matlabformat )
		fprintf ( ofile, "results.%s = %lf;\n", varname.c_str(), theta );
	else
		fprintf ( ofile, "\n# %s\n %lf\n\n", varname.c_str(), theta );
}

void print_fisher ( PsiPsychometric *pmf, std::vector<double> theta, PsiData *data, FILE* ofile, bool matlabformat ) {
	unsigned int i,j;
	unsigned int nparameters ( pmf->getNparams() );

	Matrix * fisher = pmf->ddnegllikeli ( theta, data );

	if ( matlabformat ) {
		fprintf ( ofile, "results.fisher_info = [ " );
		for ( i=0; i<nparameters; i++ ) {
			for ( j=0; j<nparameters; j++ ) {
				fprintf ( ofile, " %lf ", (*fisher)(i,j) );
			}
			if ( i==nparameters-1 )
				fprintf ( ofile, "];\n" );
			else
				fprintf ( ofile, ";" );
		}
	} else {
		fprintf ( ofile, "# fisher_info\n" );
		for ( i=0; i<nparameters; i++ ) {
			fprintf ( ofile, " %lf", (*fisher)(i,0));
			for (j=1; j<nparameters; j++) {
				fprintf ( ofile, " %lf", (*fisher)(i,j) );
			}
			fprintf ( ofile, " \n" );
		}
		fprintf ( ofile, "\n" );
	}

	delete fisher;
}

void print ( std::vector< std::vector<int> >& theta, bool matlabformat, std::string varname, FILE *ofile ) {
	unsigned i, j;

	if ( matlabformat ) {
		fprintf ( ofile, "results.%s = [ ", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			for ( j=0; j<theta[i].size()-1; j++ ) {
				fprintf ( ofile, "%d, ", theta[i][j] );
			}
			if ( i<theta.size()-1 ) {
				fprintf ( ofile, "%d; ...\n", theta[i][j] );
			} else {
				fprintf ( ofile, "%d];\n", theta[i][j] );
			}
		}
	} else {
		fprintf ( ofile, "\n# %s\n", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			for ( j=0; j<theta[i].size(); j++ ) {
				fprintf ( ofile, " %d ", theta[i][j] );
			}
			fprintf ( ofile,"\n" );
		}
		fprintf ( ofile, "\n" );
	}
}

void print ( std::vector< std::vector<double> >& theta, bool matlabformat, std::string varname, FILE *ofile ) {
	unsigned i, j;

	if ( matlabformat ) {
		fprintf ( ofile, "results.%s = [ ", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			for ( j=0; j<theta[i].size()-1; j++ ) {
				fprintf ( ofile, "%lf, ", theta[i][j] );
			}
			if ( i<theta.size()-1 ) {
				fprintf ( ofile, "%lf; ...\n", theta[i][j] );
			} else {
				fprintf ( ofile, "%lf];\n", theta[i][j] );
			}
		}
	} else {
		fprintf ( ofile, "\n# %s\n", varname.c_str() );
		for ( i=0; i<theta.size(); i++ ) {
			for ( j=0; j<theta[i].size(); j++ ) {
				fprintf ( ofile, " %lf ", theta[i][j] );
			}
			fprintf ( ofile,"\n" );
		}
		fprintf ( ofile, "\n" );
	}
}
