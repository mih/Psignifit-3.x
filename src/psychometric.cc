#include "psychometric.h"
#include "special.h"
#include "linalg.h"

// #ifdef DEBUG_PSYCHOMETRIC
#include <iostream>
// #endif

PsiPsychometric::PsiPsychometric (
	int nAFC,
	PsiCore * core,
	PsiSigmoid * sigmoid
	) : Nalternatives(nAFC), guessingrate(1./nAFC), priors( (nAFC==1 ? 4 : 3 ) )
{
	int k;
	Core = core;
	Sigmoid = sigmoid;
	for (k=0; k<priors.size(); k++)
		priors[k] = new PsiPrior;
}

PsiPsychometric::~PsiPsychometric ( void )
{
	int k;
	delete Core;
	delete Sigmoid;
	for (k=0; k<priors.size(); k++) {
		delete priors[k];
	}
}

double PsiPsychometric::evaluate ( double x, const std::vector<double>& prm ) const
{
	double gamma(guessingrate);
	if (Nalternatives==1) // Here we talk about a yes/no task
		gamma = prm[3];

#ifdef DEBUG_PSYCHOMETRIC
	std::cerr << "Evaluating psychometric function. Parameters: \n";
	std::cerr << "alpha = " << prm[0] << "\nbeta = " << prm[1] << "\nlambda = " << prm[2] << "\ngamma = " << gamma << "\n";
	std::cerr << "Value of sigmoid: " << Sigmoid->f(Core->g(x,prm)) << "\n";
#endif

	return gamma + (1-gamma-prm[2]) * Sigmoid->f(Core->g(x,prm));
}

double PsiPsychometric::negllikeli ( const std::vector<double>& prm, const PsiData* data ) const
{
	int i,n,k;
	double l(0);
	double x,p,lognoverk;

	for (i=0; i<data->getNblocks(); i++)
	{
		n = data->getNtrials(i);
		k = data->getNcorrect(i);
		x = data->getIntensity(i);
		lognoverk = data->getNoverK(i);
		p = evaluate(x, prm);
		l -= lognoverk;
		if (p>0)
			l -= k*log(p);
		else
			l += 1e10;
		if (p<1)
			l -= (n-k)*log(1-p);
		else
			l += 1e10;
	}

	return l;
}

double PsiPsychometric::leastfavourable ( const std::vector<double>& prm, const PsiData* data, double cut, bool threshold ) const
{
	if (!threshold) throw NotImplementedError();  // So far we only have this for the threshold

	std::vector<double> delta (prm.size(),0), du(prm.size(),0);
	Matrix * I;
	double ythres;
	double rz,nz,xz,pz,fac1;
	double l_LF(0);
	double c,s;
	int i,j,k,z;

	// Fill u
	ythres = Sigmoid->inv(cut);
	du[0] = Core->dinv(ythres,prm,0);
	du[1] = Core->dinv(ythres,prm,1);

	// Determine 2nd derivative
	I = ddnegllikeli ( prm, data );

	// Now we have to solve I*delta = du for delta
	delta = I->solve ( du );

	// I is not needed anymore
	delete I;

	// Normalize the result
	s = 0;
	for (i=0; i<prm.size(); i++)
		s += delta[i]*delta[i];
	s = sqrt(s);
	for (i=0; i<prm.size(); i++)
	    delta[i] /= s;

	// The result has to be multiplied by the gradient of the likelihood
	for (z=0; z<data->getNblocks(); z++) {
		rz = data->getNcorrect(z);
		nz = data->getNtrials(z);
		xz = data->getIntensity(z);
		pz = evaluate(xz,prm);
		fac1 = rz/pz - (nz-rz)/(1-pz);
		for (i=0; i<2; i++)
			l_LF += delta[i] * fac1 * Sigmoid->df(Core->g(xz,prm)) * Core->dg(xz,prm,i);
	
		for (i=2; i<prm.size(); i++)
			l_LF += delta[i] * fac1 * ( (i==2 ? 1 : 0) - Sigmoid->f(Core->g(xz,prm)) );
	}

	return l_LF;
}

Matrix * PsiPsychometric::ddnegllikeli ( const std::vector<double>& prm, const PsiData* data ) const
{
	Matrix * I = new Matrix ( prm.size(), prm.size() );

	double rz,nz,pz,xz,fac1,fac2;
	int z,i,j;

	// Fill I
	for (z=0; z<data->getNblocks(); z++) {
		rz = data->getNcorrect(z);
		nz = data->getNtrials(z);
		xz = data->getIntensity(z);
		pz = evaluate(xz,prm);
		fac1 = rz/pz - (nz-rz)/(1-pz);
		fac2 = rz/(pz*pz) + (nz-rz)/((1-pz)*(1-pz));

		// These parts must be determined
		for (i=0; i<2; i++) {
			for (j=i; j<2; j++) {
				(*I)(i,j) += fac1 * (1-guessingrate-prm[2]) * (Sigmoid->ddf(Core->g(xz,prm)) * Core->dg(xz,prm,i) * Core->dg(xz,prm,j) + Sigmoid->df(Core->g(xz,prm)) * Core->ddg(xz,prm,i,j));
				(*I)(i,j) -= fac2 * (1-guessingrate-prm[2]) * (1-guessingrate-prm[2]) * pow(Sigmoid->df(Core->g(xz,prm)),2) * Core->dg(xz,prm,i) * Core->dg(xz,prm,j);
			}
			for (j=2; j<prm.size(); j++) {
				(*I)(i,j) -= fac1 * Sigmoid->df(Core->g(xz,prm)) * Core->dg(xz,prm,i);
				(*I)(i,j) += fac2 * (1-guessingrate-prm[2]) * Sigmoid->df(Core->g(xz,prm)) * Core->dg(xz,prm,i) * ( (j==2 ? 1 : 0) - Sigmoid->f(Core->g(xz,prm)) );
			}
		}
	}

	// Average
	// TODO: Do we really need this step?
	for (i=0; i<prm.size(); i++)
		for (j=i; j<prm.size(); j++)
			(*I)(i,j) /= data->getNblocks();

	// The remaining parts of I can be copied
	for (i=1; i<prm.size(); i++)
		for (j=0; j<i; j++)
			(*I)(i,j) = (*I)(j,i);

	return I;
}

std::vector<double> PsiPsychometric::dnegllikeli ( const std::vector<double>& prm, const PsiData* data ) const
{
	std::vector<double> out (prm.size());
	double rz,xz,pz,nz,fac1;
	int z,i;

	for (z=0; z<data->getNblocks(); z++) {
		rz = data->getNcorrect(z);
		nz = data->getNtrials(z);
		xz = data->getIntensity(z);
		pz = evaluate(xz,prm);
		fac1 = rz/pz - (nz-rz)/(1-pz);
		for (i=0; i<2; i++)
			out[i] = fac1 * (1-guessingrate-prm[2]) * Sigmoid->df(Core->g(xz,prm)) * Core->dg(xz,prm,i);
	
		for (i=2; i<prm.size(); i++)
			out[i] = fac1 * ( (i==2 ? 1 : 0) - Sigmoid->f(Core->g(xz,prm)) );
	}

	return out;
}

double PsiPsychometric::deviance ( const std::vector<double>& prm, const PsiData* data ) const
{
	int i,n;
	double D(0);
	double x,y,p;

	for ( i=0; i<data->getNblocks(); i++ )
	{
		n = data->getNtrials(i);
		y = data->getPcorrect(i);
		x = data->getIntensity(i);
		p = evaluate( x, prm );
		if (y>0)
			D += n*y*log(y/p);
		if (y<1)
			D += n*(1-y)*log((1-y)/(1-p));
	}
	D *= 2;
	return D;
}

void PsiPsychometric::setPrior ( int index, PsiPrior* prior )
{
	delete priors[index];

	priors[index] = prior;
}

double PsiPsychometric::neglpost ( const std::vector<double>& prm, const PsiData* data ) const
{
	int i;
	double l;
	l = negllikeli( prm, data);

	for (i=0; i<prm.size(); i++) {
		l -= log( priors[i]->pdf(prm[i]) );
	}

	return l;
}

std::vector<double> PsiPsychometric::getStart ( const PsiData* data ) const
{
	int i;
	double a,b;
	std::vector<double> x (data->getIntensities());
	std::vector<double> p (data->getPcorrect());
	double minp(1000), maxp(-1000);
	double meanx(0), meanp(0);
	double varx(0),covxp(0);

	// Scale the data to the interval (0,1)
	for (i=0; i<x.size(); i++)
		if (minp>p[i]) minp = p[i];
#ifdef DEBUG_PSYCHOMETRIC
	std::cerr << "minp="<<minp << std::endl;
#endif

	for (i=0; i<x.size(); i++)
		p[i] -= 0.999*minp;

	for (i=0; i<x.size(); i++)
		if (maxp<p[i]) maxp = p[i];
#ifdef DEBUG_PSYCHOMETRIC
	std::cerr << "maxp="<<maxp << std::endl;
#endif

	for (i=0; i<x.size(); i++)
		p[i] /= 1.0001*maxp;

	// Apply logit
	for (i=0; i<x.size(); i++) {
#ifdef DEBUG_PSYCHOMETRIC
		std::cerr << "p["<<i<<"] = "<<p[i]<<"\t";
#endif
		p[i] = log(p[i]/(1-p[i]));
#ifdef DEBUG_PSYCHOMETRIC
		std::cerr << "lp["<<i<<"] = "<<p[i]<<"\n";
#endif
	}

	// Determine averages
	for (i=0; i<x.size(); i++) {
		meanx += x[i];
		meanp += p[i];
	}
	meanx /= x.size();
	meanp /= x.size();

	// Compute covariances
	for (i=0; i<x.size(); i++) {
		varx += (x[i]-meanx)*(x[i]-meanx);
		covxp += (x[i]-meanx)*(p[i]-meanp);
	}

	b = covxp/varx;
	a = meanp - meanx*b;

	std::vector<double> out (Core->transform( getNparams(), a,b));
	if (Nalternatives==1) {
		out[2] = 0.02;
		out[3] = .02;
	} else {
	    out[2] = 0.02;
	}

	return out;
}

std::vector<double> PsiPsychometric::getDevianceResiduals ( const std::vector<double>& prm, const PsiData* data ) const
{
	int i, n;
	double x,y,p;
	std::vector<double> out (data->getNblocks());

	for ( i=0; i<data->getNblocks(); i++ )
	{
		n = data->getNtrials(i);
		y = data->getPcorrect(i);
		x = data->getIntensity(i);
		p = evaluate(x,prm);
		out[i] = 0;
		if (y>0)
			out[i] += n*y*log(y/p);
		if (y<1)
			out[i] += n*(1-y)*log((1-y)/(1-p));
		out[i] = (y>p?1:-1) * sqrt(2*out[i]);
	}

	return out;
}

double PsiPsychometric::getRpd ( const std::vector<double>& devianceresiduals, const std::vector<double>& prm, const PsiData* data ) const {
	int k,N(data->getNblocks());
	double Ed(0),Ep(0),vard(0),varp(0),R(0);
	std::vector<double> p ( data->getPcorrect() );

	// Evaluate p values in advance
	for ( k=0; k<data->getNblocks(); k++ ) {
		p[k] = evaluate(data->getIntensity(k),prm);
	}

	// Calculate averages
	for ( k=0; k<data->getNblocks(); k++ ) {
		Ed += devianceresiduals[k];
		Ep += p[k];
	}
	Ed /= N;
	Ep /= N;

	// Calculate unnormalized variances and covariances
	for ( k=0; k<data->getNblocks(); k++ ) {
		vard += pow(devianceresiduals[k]-Ed,2);
		varp += pow(p[k]-Ep,2);
		R    += (devianceresiduals[k]-Ed)*(p[k]-Ep);
	}

	// Normalize and return
	R /= sqrt(vard);
	R /= sqrt(varp);

	return R;
}

double PsiPsychometric::getRkd ( const std::vector<double>& devianceresiduals ) const
{
	int k,N(devianceresiduals.size());
	double Ed(0), Ek(0), vard(0), vark(0), R(0);

	// Calculate averages
	for ( k=0; k<N; k++ ) {
		Ed += devianceresiduals[k];
		Ek += k;
	}
	Ed /= N;
	Ek /= N;

	// Calculate unnormalized variances and covariances
	for ( k=0; k<N; k++ ) {
		vard += pow(devianceresiduals[k]-Ed,2);
		vark += pow(k-Ek,2);
		R    += (devianceresiduals[k]-Ed)*(k-Ek);
	}

	// Normalize and return
	R /= sqrt(vard);
	R /= sqrt(vark);

	return R;
}

double PsiPsychometric::dllikeli ( std::vector<double> prm, const PsiData* data, unsigned int i ) const
{
	int k, Nblocks(data->getNblocks());
	double rz,pz,nz,xz,dl(0);
	double guess(1./Nalternatives);
	if (Nalternatives==1) // Here we talk about a yes/no task
		guess = prm[3];

	for (k=0; k<Nblocks; k++) {
		rz = data->getNcorrect ( k );
		nz = data->getNtrials  ( k );
		xz = data->getIntensity( k );
		pz = evaluate ( xz, prm );
		switch (i) {
			case 0: case 1:
				dl += (rz/pz - (nz-rz)/(1-pz)) * (1-guess-prm[2]) * Sigmoid->df ( Core->g ( xz, prm ) ) * Core->dg ( xz, prm, i );
				break;
			case 2:
				dl -= (rz/pz - (nz-rz)/(1-pz)) * Sigmoid->f ( Core->g ( xz, prm ) );
				break;
			case 3:
				if (Nalternatives==1) // gamma is a free parameter
					dl += (rz/pz - (nz-rz)/(1-pz)) * (1 - Sigmoid->f ( Core->g ( xz, prm ) ));
				break;
		}
	}

	return dl;
}

double PsiPsychometric::dlposteri ( std::vector<double> prm, const PsiData* data, unsigned int i ) const
{
	if ( i < getNparams() )
		return dllikeli ( prm, data, i ) + priors[i]->dpdf(prm[i]);
	else
		return 0;
}
