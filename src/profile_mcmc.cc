#include "psipp.h"
#include <vector>
#include <iostream>

int main ( int argc, char ** argv ) {
	PsiData *data;
	PsiCore *core = new mwCore ();
	PsiSigmoid *sigmoid = new PsiLogistic();
	PsiPsychometric * pmf = new PsiPsychometric ( 2, core, sigmoid );

	setSeed(0);

	std::vector< double > prm ( 3 );
	prm[0] = 4; prm[1] = 1; prm[2] = .03;

	unsigned int nblocks ( 6 ), i, j;
	std::vector< double > x ( nblocks );
	std::vector< int    > n ( nblocks );
	std::vector< int    > k ( nblocks );
	double p;
	std::cout << "Data:\n";
	for ( i=0; i<nblocks; i++ ) {
		x[i] = 2+4.*double(i)/(nblocks-1);
		n[i] = 30;
		p = pmf->evaluate ( x[i], prm );
		k[i] = 0;
		for ( j=0; j<n[i]; j++ ) {
			k[i] += drand48() < p;
		}
		std::cout << x[i] << " " << k[i] << " " << n[i] << "\n";
	}
	data = new PsiData ( x, n, k, 2 );

	PsiRandom *proposal = new GaussRandom ();
	MetropolisHastings MH ( pmf, data, proposal );
	MH.setStepSize ( 0, 1. );
	MH.setStepSize ( 1, 1. );
	MH.setStepSize ( 2, .02 );

	std::cout << "Getting samples ... "; std::cout.flush();
	MCMCList mcmc = MH.sample( 1000000 );
	std::cout << "Done\n";

	delete data;
	delete core;
	delete sigmoid;
	delete pmf;
	delete proposal;
}
