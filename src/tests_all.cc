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
	failures += T->isequal(pmf->getRkd(devianceresiduals),-0.320889,"OptimizerSolution 2AFC Rkd",1e-2);

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
	failures += T->isequal(pmf->getRkd(devianceresiduals),-0.477967,"OptimizerSolution Y/N Rkd",1e-2);

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
	BootstrapList boots = parametricbootstrap ( 9999, data, pmf, cuts );

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
		failures += T->conditional(!jackknife.influential(i,ci_lower,ci_upper),testname);
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
	PsiMClist post ( S->sample(10000) );
	std::cout << "\n\n";
	srand48(0);
	PsiMClist mhpost ( mhS->sample(10000) );

	failures += T->isequal ( post.getMean(0), 3.21657, "Hybrid MCMC alpha", 1e-5 );
	failures += T->isequal ( post.getMean(1), 1.20476, "Hybrid MCMC beta", 1e-5 );
	failures += T->isequal ( post.getMean(2), 0.0217217, "Hybrid MCMC lambda", 1e-5 );
	failures += T->isequal ( mhpost.getMean(0), 3.22372, "Metropolis Hastings alpha", 1e-5 );
	failures += T->isequal ( mhpost.getMean(1), 1.12734, "Metropolis Hastings beta", 1e-5 );
	failures += T->isequal ( mhpost.getMean(2), 0.0199668, "Metropolis Hastings lambda", 1e-5 );

	std::cerr << pmf->neglpost(prm,data);

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
	failures += T->ismore ( mindf, 0, "PsiLogistic monotonically increasing" );
	// Check saturation
	failures += T->isless  ( sigmoid->df( 3), 0.01, "PsiGauss->df(3)");
	failures += T->isless  ( sigmoid->df(-3), 0.01, "PsiGauss->df(-3)");
	// Check convexity
	failures += T->isless  ( sigmoid->ddf(3), 0,    "PsiGauss->ddf(3)");
	failures += T->ismore  ( sigmoid->ddf(-3),0,    "PsiGauss->ddf(-3)");
	delete sigmoid;


	return failures;
}

int CoreTests ( TestSuite * T ) {
	int failures(0);
	PsiCore * core;
	std::vector<double> prm(4);

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


	return failures;
}

int main ( int argc, char ** argv ) {
	TestSuite Tests ( "tests_all.log" );
	// Tests.addTest(&PsychometricValues,"Values of the psychometric function");
	// Tests.addTest(&OptimizerSolution, "Solutions of optimizer");
	// Tests.addTest(&BootstrapTest,     "Bootstrap properties");
	// Tests.addTest(&SigmoidTests,      "Properties of sigmoids");
	// Tests.addTest(&CoreTests,         "Tests of core objects");
	Tests.addTest(&MCMCTest,          "MCMC");
	Tests.runTests();
}
