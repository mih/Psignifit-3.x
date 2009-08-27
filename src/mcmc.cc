#include "mcmc.h"

PsiMetropolisHastings::PsiMetropolisHastings ( const PsiPsychometric * pmf, const PsiData * dat, const PsiRandom * proposal )
	: model(*pmf), data(*dat), currenttheta(pmf->getNparams(),0), newtheta(pmf->getNparams(),0), propose(proposal), stepwidths(pmf->getNparams,1),currentmomentum(pmf->getNparams()),newmomentum(pmf->getNparams()), qold(1)
{}


std::vector<double> PsiMetropolisHastings::sample ( void ) {
	double qnew, lratio, acc(log(drand48()));
	int prm, Nprm(model.getNparams);

	// propose a new point
	for (prm=0; prm<Nprm; prm++) {
		newtheta[prm] = currenttheta[prm] + stepwidths[prm] * propose->draw ( );
	}

	// negative log posterior of the point
	qnew = model.neglpost ( newtheta, &data );
	lratio = qold - qnew;

	if (acc<ratio) {
		// accept the new point
		qold = qnew;
		currenttheta = newtheta;
		currentdeviance = model.deviance ( currenttheta, &data );
	}

	return currenttheta;
}
