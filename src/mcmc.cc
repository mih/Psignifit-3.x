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
	double qnew, lratio, acc(propose->rngcall());
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

MCMCList MetropolisHastings::sample ( unsigned int N ) {
	const PsiData * data ( getData() );
	const PsiPsychometric * model ( getModel() );
	MCMCList out ( N, model->getNparams(), data->getNblocks() );
	PsiData *localdata = new PsiData ( data->getIntensities(), data->getNtrials(), data->getNcorrect(), data->getNalternatives() );
	PsiData *reduceddata;
	std::vector<int> posterior_predictive ( data->getNblocks() );
	std::vector<double> probs ( data->getNblocks() );
	std::vector<double> est ( model->getNparams() );
	int i,j,k,l;

	std::vector<double> reducedx ( data->getNblocks()-1 );
	std::vector<int> reducedk ( data->getNblocks()-1 );
	std::vector<int> reducedn ( data->getNblocks()-1 );

	for (i=0; i<N; i++) {
		// Draw the nest sample
		est = draw();
		out.setEst ( i, est, 0. );
		out.setdeviance ( i, getDeviance() );

		// determine posterior predictives
		for ( k=0; k<data->getNblocks(); k++ )
			probs[k] = model->evaluate ( data->getIntensity(k), est );
		newsample ( localdata, probs, &posterior_predictive);
		localdata->setNcorrect ( posterior_predictive );
		out.setppData ( i, posterior_predictive, model->deviance ( est, localdata ) );
		probs = model->getDevianceResiduals ( est, localdata );
		out.setRpd ( i, model->getRpd ( probs, est, data ) );
		out.setRkd ( i, model->getRkd ( probs, data ) );
		out.setppRpd ( i, model->getRpd ( probs, est, localdata ) );
		out.setppRkd ( i, model->getRkd ( probs, localdata ) );

		// Store log posterior ratios for reduced data sets
		for ( k=0; k<data->getNblocks(); k++) {
			j=0;
			for ( l=0; l<data->getNblocks(); l++ ) {
				if ( l!=k ) {
					reducedx[j] = data->getIntensity(l);
					reducedk[j] = data->getNcorrect(l);
					reducedn[j] = data->getNtrials(l);
					j++;
				}
			}
			reduceddata = new PsiData ( reducedx, reducedn, reducedk, data->getNalternatives() );
			out.setlogratio ( i, k, model->neglpost(est,data)-model->neglpost(est,reduceddata) );
			delete reduceddata;
		}
#ifdef DEBUG_MCMC
		std::cerr << " accept: " << std::setiosflags ( std::ios::fixed ) << double(accept)/(i+1) << "\n";
#endif
	}

#ifdef DEBUG_MCMC
	std::cerr << "Acceptance rate: " << double(accept)/N << "\n";
#endif

	delete localdata;

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

	if ( log(proposal->rngcall()) < currentH-newH ) {
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

MCMCList HybridMCMC::sample ( unsigned int N ) {
	MCMCList out ( N, getModel()->getNparams(), getData()->getNblocks() );
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

/**********************************************************************
 *
 * Evidence
 *
 */

double ModelEvidence ( const PsiPsychometric* pmf, const PsiData* data )
{
	std::vector<double> prm ( pmf->getNparams() );
	double E(0);
	int i,k,n(50000);

	for ( i=0; i<n; i++ ) {
		for ( k=0; k<pmf->getNparams(); k++ )
			prm[k] = pmf->randPrior ( k );

		E += exp ( - pmf->negllikeli ( prm, data ) );
	}

	E /= n;

	return E;
}

std::vector<double> OutlierDetection ( const PsiPsychometric* pmf, OutlierModel* outl, const PsiData* data )
{
	int i;
	std::vector<double> out ( data->getNblocks() );
	double E ( ModelEvidence ( pmf, data ) );

	for ( i=0; i<data->getNblocks(); i++ ) {
		outl->setexclude ( i );
		out[i] = E/ModelEvidence ( outl, data );
	}

	return out;
}
