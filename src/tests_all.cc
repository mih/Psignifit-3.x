#include <iostream>
#include <cstdlib>
#include "psychometric.h"
#include "mclist.h"
#include "bootstrap.h"
#include "testing.h"
#include "mcmc.h"

int PsychometricValues ( TestSuite* T ) {
	int failures(0),i;
	char message[40];
	std::vector <double> x ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );
	std::vector <double> p ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 29;  k[1] = 31;  k[2] = 36;  k[3] = 42;  k[4] = 46;  k[5] = 49;
	p[0] = 0.5311852; p[1] = 0.60013209; p[2] = 0.74;
	p[3] = 0.87986791;p[4] = 0.9488148;  p[5] = 0.97136662;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	PsiPsychometric * pmf = new PsiPsychometric ( 2, new abCore(), new PsiLogistic() );
	std::vector<double> prm(3);
	prm[0] = 4; prm[1] = 1.5; prm[2] = 0.02;

	// Test forward probabilities
	for ( i=0; i<6; i++ ) {
		sprintf(message,"PsychometricValues at x=%g", x[i]);
		failures += T->isequal ( pmf->evaluate(x[i],prm), p[i], message);
	}

	// Test likelihood
	failures += T->isequal ( pmf->negllikeli(prm,data), 11.996474658154325, "PsychometricValues likelihood");

	return failures;
}

int OptimizerSolution ( TestSuite * T ) {
	int failures(0),i;
	char message[40];
	double deviance(0);

	/************************
	 * 2AFC
	 */
	std::vector <double> x ( 6 );
	std::vector <int>    n ( 6, 50 );
	std::vector <int>    k ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	PsiPsychometric * pmf = new PsiPsychometric ( 2, new abCore(), new PsiLogistic() );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, new UniformPrior(0.,0.1));
	std::vector<double> start (4);
	start = pmf->getStart( data );

	// Optimizer
	PsiOptimizer *opt = new PsiOptimizer (pmf, data);

	std::vector<double> solution (4);
	solution = opt->optimize(pmf,data);
	/* To check the solution, we can use the following python script:
#!/usr/bin/env python

from numpy import *
from scipy import stats
from scipy.optimize import fmin

x = arange(6)*2
k = array([24.,32,40,48,50,48])
n = ones(k.shape)*50

# Logistic regression for starting values
p = k/n
p -= 0.999*p.min()
p /= 1.0001*p.max()
lp = log(p/(1-p))
b,a = stats.linregress(x,lp)[:2]
bt = 1./b
al = -a*bt
prm0 = array([al,bt,0.02])

def model(prm,x,k,n):
    al,bt,lm = prm
    gm = 0.5
    psi = gm + (1-gm-lm)/(1.+exp(-(x-al)/bt))
    penalty = 0.
    if lm<0 or lm>.1:
        penalty = 1e10
    return -sum(log(stats.binom.pmf(k,n,psi))) + penalty

al,bt,lm = fmin(model,prm0,args=(x,k,n))
	*/
	std::vector<double> devianceresiduals (pmf->getDevianceResiduals ( solution, data ));
	for ( i=0; i<devianceresiduals.size(); i++ ) {
		deviance += devianceresiduals[i]*devianceresiduals[i];
	}

	failures += T->isequal(start[0],3.99317998,"Model->getStart 2AFC alpha",1e-4);
	failures += T->isequal(start[1],0.88126645,"Model->getStart 2AFC beta",1e-4);

	failures += T->isequal(solution[0],3.2963951658293693,"OptimizerSolution 2AFC alpha",1e-4);
	failures += T->isequal(solution[1],0.95916747050675411,"OptimizerSolution 2AFC beta",1e-4);
	failures += T->isequal(solution[2],0.019132769792808153,"OptimizerSolution 2AFC lambda",1e-4);

	failures += T->isequal(pmf->deviance(solution,data),3.98476,"OptimizerSolution 2AFC deviance",1e-4);
	failures += T->isequal(pmf->deviance(solution,data),deviance,"OptimizerSolution 2AFC deviance sum", 1e-7);

	failures += T->isequal(pmf->getRpd(devianceresiduals,solution,data),0.155395,"OptimizerSolution 2AFC Rpd",1e-2);
	// The following test fails for some reason.
	// Why is that? Python gives the same correlation between index and deviance residuals.
	// The number used here is from psignifit.
	failures += T->isequal(pmf->getRkd(devianceresiduals,data),-0.320889,"OptimizerSolution 2AFC Rkd",1e-2);

	delete pmf;
	delete data;
	delete opt;

	/************************
	 * Yes/No
	 */
	std::clog << "\n";

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 3;  k[1] = 10;  k[2] = 34;  k[3] = 45;  k[4] = 50;  k[5] = 50;
	data = new PsiData (x,n,k,2);

	// Set up psychometric function
	pmf = new PsiPsychometric ( 1, new abCore(), new PsiLogistic() );
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02; prm[3] = 0.1;
	pmf->setPrior( 2, new UniformPrior(0.,0.1));
	// pmf->setPrior( 3, new UniformPrior(0.,1.));
	start = pmf->getStart( data );

	// Optimizer
	opt = new PsiOptimizer (pmf, data);

	solution = opt->optimize(pmf,data);
	/* To check the solution, we can use the following python script:
#!/usr/bin/env python

from numpy import *
from scipy import stats
from scipy.optimize import fmin

x = arange(6)*2
k = array([ 3, 10, 34, 45, 50, 50])
n = ones(k.shape)*50

# Logistic regression for starting values
p = k/n
p -= 0.999*p.min()
p /= 1.0001*p.max()
lp = log(p/(1-p))
b,a = stats.linregress(x,lp)[:2]
bt = 1./b
al = -a*bt
prm0 = array([al,bt,0.02,0.5])

def model(prm,x,k,n):
    al,bt,lm,gm = prm
    psi = gm + (1-gm-lm)/(1.+exp(-(x-al)/bt))
    penalty = 0.
    if lm<0 or lm>.1:
        penalty = 1e10
    return -sum(log(stats.binom.pmf(k,n,psi))) + penalty

al,bt,lm,gm = fmin(model,prm0,args=(x,k,n))
	*/
	devianceresiduals = pmf->getDevianceResiduals ( solution, data );
	deviance = 0;
	for ( i=0; i<devianceresiduals.size(); i++ ) {
		deviance += devianceresiduals[i]*devianceresiduals[i];
	}

	failures += T->isequal(start[0],4.11079579,"Model->getStart Y/N alpha",1e-4);
	failures += T->isequal(start[1],0.54401123599182111,"Model->getStart Y/N beta",1e-4);

	failures += T->isequal(solution[0],3.4435404356684396,"OptimizerSolution Y/N alpha",2*1e-2);
	failures += T->isequal(solution[1],0.99142691236758229,"OptimizerSolution Y/N beta",2*1e-2);
	failures += T->isequal(solution[2],8.4692705817775732e-09,"OptimizerSolution Y/N lambda",1e-2);
	failures += T->isequal(solution[3],0.029071443401547103,"OptimizerSolution Y/N gamma",1e-2);

	failures += T->isequal(pmf->deviance(solution,data),2.08172,"OptimizerSolution Y/N deviance",2*1e-1);
	failures += T->isequal(pmf->deviance(solution,data),deviance,"OptimizerSolution Y/N deviance sum", 1e-7);

	failures += T->isequal(pmf->getRpd(devianceresiduals,solution,data),0.217146,"OptimizerSolution Y/N Rpd",1e-1);
	failures += T->isequal(pmf->getRkd(devianceresiduals,data),-0.477967,"OptimizerSolution Y/N Rkd",1e-2);

	delete pmf;
	delete data;
	delete opt;

	return failures;
}

int BootstrapTest ( TestSuite * T ) {
	srand48( 0 );
	int failures(0),i;
	std::vector<double> x ( 6 );
	std::vector<int>    n ( 6, 50 );
	std::vector<int>    k ( 6 );
	char testname[40];

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	PsiPsychometric * pmf = new PsiPsychometric ( 2, new abCore(), new PsiLogistic() );
	std::vector<double> prm(4);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, new UniformPrior(0.,0.1));

	std::vector<double> cuts (1, 0.5);
	BootstrapList boots = bootstrap ( 9999, data, pmf, cuts );

	// Check against psignifit results
	// These values are subject to statistical variation. "equality" is defined relatively coarse
	failures += T->isless(boots.getAcc(0),0.05,"Acceleration constant");
	failures += T->isequal(boots.getBias(0),-.08,"Bias",.1);
	failures += T->isequal(boots.getThres(.1,0),2.63745,"u(.1)",.1);
	failures += T->isequal(boots.getThres(.9,0),3.84851,"u(.9)",.1);
	failures += T->isequal(boots.getDeviancePercentile(0.975),9.2,"Deviance limits",.5);
	failures += T->isequal(boots.percRpd(.025), -0.534564, "Rpd( 2.5%)", .1); // Testing mean and standard error
	failures += T->isequal(boots.percRpd(.975), 0.587535, "Rpd(97.5%)",  .1);
	failures += T->isequal(boots.percRkd(.025), -0.925168, "Rkd( 2.5%)", .1);
	failures += T->isequal(boots.percRkd(.975), 0.612778, "Rkd(97.5%)",  .1);

	// Check for influential observations and outliers
	JackKnifeList jackknife = jackknifedata (data, pmf);

	std::vector<double> ci_lower ( pmf->getNparams() ), ci_upper ( pmf->getNparams() );
	for ( i=0; i<pmf->getNparams(); i++ ) {
		ci_lower[i] = boots.getPercentile(0.025,i);
		ci_upper[i] = boots.getPercentile(0.975,i);
	}

	for ( i=0; i<6; i++ ) {
		sprintf(testname,"influential %d",i);
		failures += T->conditional(jackknife.influential(i,ci_lower,ci_upper)<1,testname);
		sprintf(testname,"outliers %d",i);
		failures += T->conditional(!jackknife.outlier(i),testname);
	}

	return failures;
}

int MCMCTest ( TestSuite * T ) {
	int failures ( 0 );

	std::vector<double> x ( 6 );
	std::vector<int>    n ( 6, 50 );
	std::vector<int>    k ( 6 );

	// Set up data
	x[0] =  0.; x[1] =  2.; x[2] =  4.; x[3] =  6.; x[4] =  8.; x[5] = 10.;
	k[0] = 24;  k[1] = 32;  k[2] = 40;  k[3] = 48;  k[4] = 50;  k[5] = 48;
	PsiData * data = new PsiData (x,n,k,2);

	// Set up psychometric function
	PsiPsychometric * pmf = new PsiPsychometric ( 2, new abCore(), new PsiLogistic() );
	std::vector<double> prm(3);
	prm[0] = 4; prm[1] = 0.8; prm[2] = 0.02;
	pmf->setPrior( 2, new UniformPrior(0.,1));

	MetropolisHastings * mhS = new MetropolisHastings( pmf, data, new GaussRandom() );
	mhS->setTheta( prm );
	mhS->setstepsize(0.1,0);
	mhS->setstepsize(0.1,1);
	mhS->setstepsize(0.001,2);

	HybridMCMC * S = new HybridMCMC ( pmf, data, 20 );
	S->setTheta ( prm );
	S->setstepsize ( 0.001, 0 );
	S->setstepsize ( 0.001, 1 );
	S->setstepsize ( 0.0001, 2 );

	srand48(0);
	MCMCList post ( S->sample(10000) );
	srand48(0);
	MCMCList mhpost ( mhS->sample(10000) );

	failures += T->isequal ( post.getMean(0), 3.21657, "Hybrid MCMC alpha", .3 );
	failures += T->isequal ( post.getMean(1), 1.20476, "Hybrid MCMC beta", .2 );
	failures += T->isequal ( post.getMean(2), 0.0217217, "Hybrid MCMC lambda", .02 );
	failures += T->isequal ( mhpost.getMean(0), 3.22372, "Metropolis Hastings alpha", .2 );
	failures += T->isequal ( mhpost.getMean(1), 1.12734, "Metropolis Hastings beta", .2 );
	failures += T->isequal ( mhpost.getMean(2), 0.0199668, "Metropolis Hastings lambda", .02 );

	return failures;
}

int PriorTest ( TestSuite * T ) {
	int failures ( 0 );
	PsiPrior * prior;

	prior = new PsiPrior;
	failures += T->isequal ( prior->pdf ( 0 ), 1, "Flat prior at 0" );
	failures += T->isequal ( prior->dpdf ( 0 ) , 0, "Flat prior derivative at 0" );
	delete prior;

	prior = new UniformPrior ( 0, 1 );
	failures += T->isequal ( prior->pdf ( -.5 ) , 0, "Uniform prior at -0.5" );
	failures += T->isequal ( prior->pdf ( .5 ) , 1,  "Uniform prior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.5 ) , 0, "Uniform prior at 1.5" );
	failures += T->isequal ( prior->dpdf ( -.5 ) , 0, "Uniform prior derivative at -0.5" );
	failures += T->isequal ( prior->dpdf ( .5 ) , 0,  "Uniform prior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.5 ) , 0, "Uniform prior derivative at 1.5" );
	delete prior;

	prior = new GaussPrior ( 0, 1 );
	failures += T->isequal ( prior->pdf ( -1 ), 0.24197072, "Gaussian prior at -1" );
	failures += T->isequal ( prior->pdf ( 0 ), 0.39894228, "Gaussian prior at 0" );
	failures += T->isequal ( prior->pdf ( 1 ), 0.24197072, "Gaussian prior at 1" );
	failures += T->isequal ( prior->dpdf ( -1 ), 0.24197072, "Gaussian prior derivative at -1" );
	failures += T->isequal ( prior->dpdf ( 0 ), 0, "Gaussian prior derivative at 0" );
	failures += T->isequal ( prior->dpdf ( 1 ), -0.24197072, "Gaussian prior derivative at 1" );
	delete prior;

	prior = new BetaPrior ( 1.5, 3. );
	failures += T->isequal ( prior->pdf ( -.1 ), 0, "BetaPrior at 0" );
	failures += T->isequal ( prior->pdf ( .1 ), 1.68094822, "BetaPrior at 0.1" );
	failures += T->isequal ( prior->pdf ( .5 ), 1.16009706, "BetaPrior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.1 ), 0, "BetaPrior at 1.1" );
	failures += T->isequal ( prior->dpdf ( -.1 ), 0, "BetaPrior derivative at 0" );
	failures += T->isequal ( prior->dpdf ( .1 ), 12.14018158, "BetaPrior derivative at 0.1" );
	failures += T->isequal ( prior->dpdf ( .5 ), 5.80048531, "BetaPrior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.1 ), 0, "BetaPrior derivative at 1.1" );
	delete prior;

	prior = new GammaPrior ( 1.5, 3. );
	failures += T->isequal ( prior->pdf ( -0.5 ), 0., "GammaPrior at -0.5" );
	failures += T->isequal ( prior->pdf ( 0.5 ), 0.12997977, "GammaPrior at 0.5" );
	failures += T->isequal ( prior->pdf ( 1.0 ), 0.15559955, "GammaPrior at 1.0" );
	failures += T->isequal ( prior->pdf ( 1.5 ), 0.16131382, "GammaPrior at 1.5" );
	failures += T->isequal ( prior->dpdf ( -0.5 ), 0., "GammaPrior derivative at -0.5" );
	failures += T->isequal ( prior->dpdf ( 0.5 ), 0.08665318, "GammaPrior derivative at 0.5" );
	failures += T->isequal ( prior->dpdf ( 1.0 ), 0.02593326, "GammaPrior derivative at 1.0" );
	failures += T->isequal ( prior->dpdf ( 1.5 ), 0., "GammaPrior derivative at 1.5" );
	delete prior;

	return failures;
}

int SigmoidTests ( TestSuite * T ) {
	int failures(0);
	PsiSigmoid * sigmoid;

	// Check gaussian cdf themselves
	failures += T->isequal(Phi(0),.5,"Phi(0)",1e-5);
	failures += T->isequal(invPhi(0.5),0.,"invPhi(0.5)",1e-5);
	failures += T->isequal(invPhi(Phi(.3)),.3,"invPhi(Phi(0.3))",1e-5);
	failures += T->isequal(Phi(invPhi(0.3)),.3,"Phi(invPhi(0.3))",1e-5);

	sigmoid = new PsiLogistic ();
	// Check specific function values
	// f should be 0.5 at 0 and close to 0 resp. 1 at low resp. high values
	failures += T->isequal ( sigmoid->f ( 0), 0.5,  "PsiLogistic->f( 0)" );
	failures += T->isless  ( sigmoid->f (-3), 0.05, "PsiLogistic->f(-3)" );
	failures += T->ismore  ( sigmoid->f ( 3), 0.95, "PsiLogistic->f( 3)" );
	// Check symmetry
	failures += T->isequal ( sigmoid->f ( 3), 1-sigmoid->f(-3), "PsiLogistic( 3)-PsiLogistic(-3)" );
	// Check saturation
	// f should have low derivative at values far from 0
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiLogistic->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiLogistic->df(-3)");
	// Check monotonicity
	double mindf(1e20),x,df;
	for ( x=-5; x<5; x+=0.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiLogistic monotonically increasing" );
		
	// Check convexity
	// if x>0, ddf<0
	// if x<0, ddf>0
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiLogistic->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiLogistic->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGauss ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f ( 0), 0.5,  "PsiGauss->f( 0)" );
	failures += T->isless  ( sigmoid->f (-3), 0.01, "PsiGauss->f(-3)" );
	failures += T->ismore  ( sigmoid->f ( 3), 0.99, "PsiGauss->f( 3)" );
	// Check symmetry
	failures += T->isequal ( sigmoid->f ( 3), 1-sigmoid->f(-3), "PsiGauss( 3)-PsiGauss(-3)" );
	// Check monotonicity
	mindf = 1e20;
	for ( x=-5; x<5; x+=0.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGaussian monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.01, "PsiGauss->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.01, "PsiGauss->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGauss->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGauss->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGumbelL ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f (0), 0.63212055882855767, "PsiGumbelL->f(0)");
	failures += T->isequal ( sigmoid->f (-3), .048568007099546562, "PsiGumbelL->f(0)");
	failures += T->isequal ( sigmoid->f (3), .99999999810782125, "PsiGumbelL->f(0)");
	// Check asymmetry
	failures += T->ismore ( sigmoid->f ( 3 ), 1-sigmoid->f( -3 ), "PsiGumbelL( 3 )-PsiGumbelL (-3 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGumbelL monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.01, "PsiGumbelL->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiGumbelL->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGumbelL->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGumbelL->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiGumbelR ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f (0), 0.36787944117144233, "PsiGumbelR->f(0)");
	failures += T->isequal ( sigmoid->f (-3), 1.8921786948382924e-09, "PsiGumbelR->f(0)");
	failures += T->isequal ( sigmoid->f (3), .95143199290045344, "PsiGumbelR->f(0)");
	// Check asymmetry
	failures += T->ismore ( 1-sigmoid->f ( -3 ), sigmoid->f( 3 ), "PsiGumbelR( -3 )-PsiGumbelR (3 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 )
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	failures += T->ismore ( mindf, 0, "PsiGumbelR monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiGumbelR->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.01, "PsiGumbelR->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGumbelR->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGumbelR->ddf(-3)");
	delete sigmoid;

	sigmoid = new PsiCauchy ();
	// Check specific function values
	failures += T->isequal ( sigmoid->f ( 0 ), 0.5, "PsiCauchy->f(0)" );
	double zalpha = -2*tan(M_PI*(0.1-0.5));
	failures += T->isequal ( sigmoid->f(-0.5*zalpha), 0.1, "PsiCauchy->f(-z(0.1)*(-.5) )" );
	failures += T->isequal ( sigmoid->f(0.5*zalpha),  0.9, "PsiCauchy->f( z(0.9)*.5 )" );
	// Check monotonicity
	mindf = 1e20;
	for (x=-5; x<5; x+=.1 ) {
		if ( (df=sigmoid->df(x))<mindf )
			mindf = df;
	}
	failures += T->ismore ( mindf, 0, "PsiCauchy monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.05, "PsiCauchy->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.05, "PsiCauchy->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiCauchy->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiCauchy->ddf(-3)");
	delete sigmoid;


	return failures;
}

int CoreTests ( TestSuite * T ) {
	int failures(0);
	PsiCore * core;
	PsiData * data;
	std::vector<double> *x;
	std::vector<int> *k,*n;
	int i;
	double th;
	std::vector<double> prm(4,0);
	std::vector<double> prm2(4,0);

	core = new abCore;
	prm[0] = 3.;
	prm[1] = 2.;
	failures += T->isequal(core->g(3.,prm),0,            "abCore at threshold");
	failures += T->isequal(core->dg(3.,prm,0),-1./prm[1],"abCore derivative 0");
	failures += T->isequal(core->dg(3.,prm,1),0,         "abCore derivative 1");
	failures += T->isequal(core->ddg(3.,prm,0,0),0,      "abCore 2nd derivative 0,0");
	failures += T->isequal(core->ddg(3.,prm,0,1),1./(prm[1]*prm[1]),"abCore 2nd derivative 0,1");
	failures += T->isequal(core->ddg(3.,prm,1,1),0,      "abCore 2nd derivative 1,1");
	failures += T->isequal(core->g(core->inv(2.,prm),prm),2, "abCore inversion g(inv(2))");
	failures += T->isequal(core->inv(core->g(2.,prm),prm),2, "abCore inversion inv(g(2))");
	failures += T->isequal(core->dinv(2.,prm,0),1.,"abCore inversion dinv(2,0)");
	failures += T->isequal(core->dinv(2.,prm,1),2.,"abCore inversion dinv(2,1)");
	// TODO: Transform tests
	delete core;

	core = new mwCore (1,0.1);
	prm[0] = 3.;
	prm[1] = 2.;
	failures += T->isequal(core->g(3.,prm),0,            "mwCore at threshold");
	failures += T->isequal(core->dg(3.,prm,0),log(9.),"mwCore derivative 0");
	failures += T->isequal(core->dg(3.,prm,1),0,         "mwCore derivative 1");
	failures += T->isequal(core->ddg(3.,prm,0,0),0,      "mwCore 2nd derivative 0,0");
	failures += T->isequal(core->ddg(3.,prm,0,1),0.5*log(9.),"mwCore 2nd derivative 0,1");
	failures += T->isequal(core->ddg(3.,prm,1,1),0,      "mwCore 2nd derivative 1,1");
	failures += T->isequal(core->g(core->inv(2.,prm),prm),2, "mwCore inversion g(inv(2))");
	failures += T->isequal(core->inv(core->g(2.,prm),prm),2, "mwCore inversion inv(g(2))");
	failures += T->isequal(core->dinv(2.,prm,0),1.,"mwCore inversion dinv(2,0)");
	failures += T->isequal(core->dinv(2.,prm,1),1./log(9.),"mwCore inversion dinv(2,1)");
	// TODO: Transform Tests
	delete core;

	core = new linearCore;
	prm2 = prm;
	th = -2./3;
	failures += T->isequal ( core->g(th,prm), 0,                   "linearCore at threshold");
	failures += T->isequal ( core->dg(th,prm,0), th,               "linearCore derivative(0) at threshold");
	failures += T->isequal ( core->dg(th,prm,1), 1.,               "linearCore derivative(1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,0), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,0), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,1), 0,             "linearCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,1), core->ddg(th,prm,1,0), "linearCore 2nd derivative symmetry at threshold");
	failures += T->isequal ( core->inv(0,prm),th,                  "linearCore inverse threshold");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm), 2.,   "linearCore inverse g(inv(2))");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm), 2.,   "linearCore inverse inv(g(2))");
	prm2[0] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,0), (core->inv(2,prm2)-core->inv(2.,prm))/1e-8,              "linearCore inverse derivative(0)");
	prm2[0] = prm[0];
	prm2[1] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,1), (core->inv(2,prm2)-core->inv(2.,prm))/1e-8,              "linearCore inverse derivative(1)");
	prm2[1] = prm[1];
	delete core;

	x = new std::vector<double> (6,0);
	k = new std::vector<int> (6,0);
	n = new std::vector<int> (6,50);
	for (i=0; i<6; i++) (*x)[i] = 2*i+.1;
	(*k)[0] = 24; (*k)[1] = 32; (*k)[2] = 40; (*k)[3] = 48; (*k)[4] = 50; (*k)[5] = 48;
	data = new PsiData ( *x, *n, *k, 2 );
	core = new logCore ( data );
	th = exp(-2./3);

	failures += T->isequal ( core->g(th,prm), 0,                  "logCore at threshold");
	failures += T->isequal ( core->dg(th,prm,0), log(th),         "logCore derivative(0) at threshold");
	failures += T->isequal ( core->dg(th,prm,1), 1,               "logCore derivative(1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,0), 0,            "logCore derivative(0,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,0), 0,            "logCore derivative(1,0) at threshold");
	failures += T->isequal ( core->ddg(th,prm,1,1), 0,            "logCore derivative(1,1) at threshold");
	failures += T->isequal ( core->ddg(th,prm,0,1), core->ddg(th,prm,1,0), "logCore 2nd derivative symmetry at threshold");
	failures += T->isequal ( core->inv(0,prm), th,                "logCore inverse");
	failures += T->isequal ( core->inv(core->g(2.,prm),prm),2.,   "logCore inversion g(inv(2))");
	failures += T->isequal ( core->g(core->inv(2.,prm),prm),2.,   "logCore inversion inv(g(2))");
	prm2[0] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,0), (core->inv(2,prm2)-core->inv(2,prm))/1e-8,            "logCore inversion dinv(2,0)");
	prm2[0] = prm[0];
	prm2[1] += 1e-8;
	failures += T->isequal ( core->dinv(2.,prm,1), (core->inv(2,prm2)-core->inv(2,prm))/1e-8, "logCore inversion dinv(2,1)");
	prm2[1] = prm[1];

	delete data;
	delete core;
	delete x;
	delete k;
	delete n;

	return failures;
}

int LinalgTests ( TestSuite * T ) {
	// These tests compare the results with the respective numpy/scipy routines
	int failures (0);

	Matrix *M = new Matrix (3,3);
	std::vector<double> x(3),b(3);

	(*M)(0,0) = 0.75; (*M)(0,1) = 0.52; (*M)(0,2) = -.16;
	(*M)(1,0) = 0.52; (*M)(1,1) = 1.38; (*M)(1,2) = -.42;
	(*M)(2,0) = -.16; (*M)(2,1) = -.42; (*M)(2,2) = 0.53;

	Matrix *I = M->inverse();
	failures += T->isequal ( (*I)(0,0),  1.80488979, "Inverse (0,0)" );
	failures += T->isequal ( (*I)(1,0), -0.67772799, "Inverse (1,0)" );
	failures += T->isequal ( (*I)(2,0),  0.00780493, "Inverse (2,0)" );
	failures += T->isequal ( (*I)(0,1), -0.67772799, "Inverse (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.20943876, "Inverse (1,1)" );
	failures += T->isequal ( (*I)(2,1),  0.75382604, "Inverse (2,1)" );
	failures += T->isequal ( (*I)(0,2),  0.00780493, "Inverse (0,2)" );
	failures += T->isequal ( (*I)(1,2),  0.75382604, "Inverse (1,2)" );
	failures += T->isequal ( (*I)(2,2),  2.48652024, "Inverse (2,2)" );
	delete I;

	b[0] = 1; b[1] = 0.5; b[2] = 0;
	x = M->solve(b);
	failures += T->isequal ( x[0],  1.4660258,  "solving Ax=b, x[0]" );
	failures += T->isequal ( x[1], -0.0730086,  "solving Ax=b, x[1]" );
	failures += T->isequal ( x[2],  0.38471795, "solving Ax=b, x[2]" );

	I = M->cholesky_dec ();
	failures += T->isequal ( (*I)(0,0),  0.8660254,  "Cholesky (0,0)" );
	failures += T->isequal ( (*I)(1,0),  0.60044428, "Cholesky (1,0)" );
	failures += T->isequal ( (*I)(2,0), -0.18475209, "Cholesky (2,0)" );
	failures += T->isequal ( (*I)(0,1),  0.        , "Cholesky (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.00968642, "Cholesky (1,1)" );
	failures += T->isequal ( (*I)(2,1), -0.30610164, "Cholesky (2,1)" );
	failures += T->isequal ( (*I)(0,2),  0.        , "Cholesky (0,2)" );
	failures += T->isequal ( (*I)(1,2),  0.        , "Cholesky (1,2)" );
	failures += T->isequal ( (*I)(2,2),  0.63416753, "Cholesky (2,2)" );
	delete I;

	I = M->lu_dec ();
	failures += T->isequal ( (*I)(0,0),  0.75      , "LU (0,0)" );
	failures += T->isequal ( (*I)(1,0),  0.69333333, "LU (1,0)" );
	failures += T->isequal ( (*I)(2,0), -0.21333333, "LU (2,0)" );
	failures += T->isequal ( (*I)(0,1),  0.52      , "LU (0,1)" );
	failures += T->isequal ( (*I)(1,1),  1.01946667, "LU (1,1)" );
	failures += T->isequal ( (*I)(2,1), -0.30316505, "LU (2,1)" );
	failures += T->isequal ( (*I)(0,2), -0.16      , "LU (0,2)" );
	failures += T->isequal ( (*I)(1,2), -0.30906667, "LU (1,2)" );
	failures += T->isequal ( (*I)(2,2),  0.40216845, "LU (2,2)" );
	delete I;

	b = (*M)*x;
	failures += T->isequal ( b[0], 1., "Ax=b, b[0]" );
	failures += T->isequal ( b[1], .5, "Ax=b, b[1]" );
	failures += T->isequal ( b[2], 0., "Ax=b, b[2]" );

	T->isequal ( M->symmetric(), 1., "M should be symmetric" );

	int i,j;
	I = new Matrix (3,3);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		(*I)(i,j) = (*M)(i,j);
	I->scale(2);
	for (i=0; i<3; i++)
		for (j=0; j<3; j++)
			failures += T->isequal ( (*I)(i,j), 2*(*M)(i,j), "matrix scaling" );
	delete I;

	delete M;

	return failures;
}

int main ( int argc, char ** argv ) {
	TestSuite Tests ( "tests_all.log" );
	Tests.addTest(&PsychometricValues,"Values of the psychometric function");
	Tests.addTest(&OptimizerSolution, "Solutions of optimizer");
	Tests.addTest(&BootstrapTest,     "Bootstrap properties");
	Tests.addTest(&SigmoidTests,      "Properties of sigmoids");
	Tests.addTest(&CoreTests,         "Tests of core objects");
	Tests.addTest(&MCMCTest,          "MCMC");
	Tests.addTest(&PriorTest,         "Priors");
	Tests.addTest(&LinalgTests,       "Linear algebra routines");
	Tests.runTests();
}
