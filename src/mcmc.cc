#include "mcmc.h"

#ifdef DEBUG_MCMC
#include <iostream>
#include <iomanip>
#endif

/**********************************************************************
 *
 * MetropolisHastings sampling
 *
 */

MetropolisHastings::MetropolisHastings ( const PsiPsychometric * pmf, const PsiData * dat, PsiRandom * proposal )
	: PsiSampler ( pmf, dat ),
	currenttheta(pmf->getNparams(),0),
	newtheta(pmf->getNparams(),0),
	propose(proposal),
	stepwidths(pmf->getNparams(),.1),
	qold(1),
	accept(0)
{
	setTheta ( currenttheta );
}

std::vector<double> MetropolisHastings::draw ( void ) {
	double qnew, lratio, acc(drand48());
	const PsiPsychometric * model (getModel());
	const PsiData * data (getData());
	int prm, Nprm(model->getNparams());

	// propose a new point
	for (prm=0; prm<Nprm; prm++) {
		newtheta[prm] = currenttheta[prm] + stepwidths[prm] * propose->draw ( );
	}

	// negative log posterior of the point
	qnew = model->neglpost ( newtheta, data );
	lratio = exp(qold - qnew);
	if (lratio>1) lratio = 1;

	if (acc<lratio) {
		// accept the new point
		qold = qnew;
		currenttheta = newtheta;
		currentdeviance = model->deviance ( currenttheta, data );
		accept ++;
#ifdef DEBUG_MCMC
		std::cerr << " * ";
#endif
	}
#ifdef DEBUG_MCMC
	else
		std::cerr << "   ";

	std::cerr
		<< "Q_old: " << std::setiosflags ( std::ios::fixed ) << qold
		<< " Q_new: " << std::setiosflags ( std::ios::fixed ) << qnew
		<< " P(accept):" << std::setiosflags ( std::ios::fixed ) << exp(qold-qnew);

	std::cout << currenttheta[0] << "\n";
#endif

	return currenttheta;
}

void MetropolisHastings::setTheta ( const std::vector<double>& prm ) {
	if (prm.size()==currenttheta.size())
		currenttheta = prm;
	else
		throw BadArgumentError();
	qold = getModel()->neglpost( currenttheta, getData() );
}

void MetropolisHastings::setstepsize ( double size, unsigned int param ) {
	if ( param<getModel()->getNparams() )
		stepwidths[param] = size;
	else
		throw BadIndexError();
}

void MetropolisHastings::setstepsize ( const std::vector<double>& sizes ) {
	int i;
	for (i=0; i<stepwidths.size(); i++)
		stepwidths[i] = sizes[i];
}

PsiMClist MetropolisHastings::sample ( unsigned int N ) {
	const PsiPsychometric * model ( getModel() );
	PsiMClist out ( N, model->getNparams() );
	int i;

	for (i=0; i<N; i++) {
		out.setEst ( i, draw(), 0. );
		out.setdeviance ( i, getDeviance() );
#ifdef DEBUG_MCMC
		std::cerr << " accept: " << std::setiosflags ( std::ios::fixed ) << double(accept)/(i+1) << "\n";
#endif
	}

#ifdef DEBUG_MCMC
	std::cerr << "Acceptance rate: " << double(accept)/N << "\n";
#endif

	return out;
}

/**********************************************************************
 *
 * Hybird MCMC
 *
 */

HybridMCMC::HybridMCMC ( const PsiPsychometric* Model, const PsiData* Data, int Nleap )
	: PsiSampler ( Model, Data ),
	currenttheta ( Model->getStart( Data ) ),
	newtheta     ( Model->getNparams(), 0 ),
	momentum  ( Model->getNparams(), 0 ),
	gradient (     Model->getNparams(), 0 ),
	currentgradient ( Model->getNparams(), 0 ),
	stepsizes (    Model->getNparams(), .001),
	Naccepted (0),
	Nleapfrog(Nleap)
{
	proposal = new GaussRandom;

	setTheta ( currenttheta );

	stepsizes[0] = 0.001;
	stepsizes[1] = 0.001;
	stepsizes[2] = 0.0001;
}

std::vector<double> HybridMCMC::draw ( void ) {
	int i;
	const PsiPsychometric * model ( getModel() );
	const PsiData *         data  ( getData() );

	for (i=0; i<model->getNparams(); i++)
		momentum[i] = proposal->draw();

	currentH = 0;
	for (i=0; i<model->getNparams(); i++)
		currentH += momentum[i] * momentum[i];
	currentH *= 0.5;
	currentH += energy;

	leapfrog();

	newenergy = model->neglpost ( newtheta, data );
	newH = 0;
	for (i=0; i<model->getNparams(); i++)
		newH += momentum[i] * momentum[i];
	newH *= 0.5;
	newH += newenergy;

	if ( log(drand48()) < currentH-newH ) {
		// Accept
		for (i=0; i<model->getNparams(); i++) {
			currenttheta[i] = newtheta[i];
			currentgradient[i] = gradient[i];
		}
		energy = newenergy;
		Naccepted ++;
#ifdef DEBUG_MCMC
		std::cerr << " * ";
#endif
	}
#ifdef DEBUG_MCMC
	else
		std::cerr << "   ";
	std::cout << currenttheta[0] << "\n";
#endif
	return currenttheta;
}

void HybridMCMC::setTheta ( const std::vector<double>& theta ) {
	int i;
	currenttheta = theta;

	for (i=0; i<getModel()->getNparams(); i++) {
		gradient[i] = getModel()->dlposteri ( currenttheta, getData(), i );
	}
	energy = getModel()->neglpost ( currenttheta, getData() );
}

void HybridMCMC::setstepsize ( const std::vector<double>& sizes ) {
	if (sizes.size()==stepsizes.size())
		stepsizes = sizes;
	else
		throw BadArgumentError();
}

void HybridMCMC::setstepsize ( double size, unsigned int param ) {
	if ( param>=stepsizes.size() )
		throw BadIndexError();

	stepsizes[param] = size;
}

void HybridMCMC::leapfrog ( void ) {
	int i,n;
	int Nparams(getModel()->getNparams());
	const PsiPsychometric * model (getModel());

	gradient = currentgradient;
	newtheta = currenttheta;

	for (n=0; n<Nleapfrog; n++) {
		for (i=0; i<Nparams; i++)
			momentum[i] -= 0.5 * stepsizes[i] * gradient[i];

		for (i=0; i<Nparams; i++)
			newtheta[i] +=          stepsizes[i] * momentum[i];

		for (i=0; i<Nparams; i++)
			gradient[i] = model->dlposteri ( newtheta, getData(), i );

		for (i=0; i<Nparams; i++)
			momentum[i] -= 0.5 * stepsizes[i] * gradient[i];
	}
}

double HybridMCMC::getDeviance ( void ) {
	return getModel()->deviance ( currenttheta, getData() );
}

PsiMClist HybridMCMC::sample ( unsigned int N ) {
	PsiMClist out ( N, getModel()->getNparams() );
	int i;

	for (i=0; i<N; i++) {
		out.setEst ( i, draw(), 0. );
		out.setdeviance ( i, getDeviance() );
#ifdef DEBUG_MCMC
		std::cerr
			<< "H: "       << std::setiosflags(std::ios::fixed) << currentH
			<< " H_: "     << std::setiosflags(std::ios::fixed) << newH
			<< " E: "      << std::setiosflags(std::ios::fixed) << energy
			<< " accept: " << std::setiosflags(std::ios::fixed) << double(Naccepted)/(i+1) << "\n";
#endif
	}

#ifdef DEBUG_MCMC
	std::cerr << "Acceptance rate: " << double(Naccepted)/N << "\n";
#endif

	return out;
}
