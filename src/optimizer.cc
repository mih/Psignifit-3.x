#include "optimizer.h"

#ifdef DEBUG_OPTIMIZER
#include <iostream>
#endif

const double maxstep (1e-7);
const double maxfstep(1e-7);
const int    maxiter (80);

PsiOptimizer::PsiOptimizer ( const PsiPsychometric * model, const PsiData * data)
	: nparameters ( model->getNparams() ),
	simplex ( nparameters+1, std::vector<double> (nparameters) ),fx ( nparameters+1 ),
	x  ( nparameters ),
	xx ( nparameters ),
	start ( nparameters ),
	modified ( nparameters+1, true )
{}

PsiOptimizer::~PsiOptimizer ( void ) {}

double testfunction(const std::vector<double>& x) {
	double out(0);
	int k;
	for (k=0; k<x.size(); k++)
		out += x[k]*x[k];
	return out;
}

std::vector<double> PsiOptimizer::optimize ( const PsiPsychometric * model, const PsiData * data )
{
	start = model->getStart(data);

	int k, l, reoptimization;

	for ( k=0; k<nparameters+1; k++ ) {
		for ( l=0; l<nparameters; l++)
			simplex[k][l] = start[l];
		modified[k] = true;
	}

#ifdef DEBUG_OPTIMIZER
	std::cerr << "Starting values for optimization:\n";
	for (k=0; k<nparameters; k++)
		std::cerr << start[k] << "\n";
#endif


	double dl;

	double ffx;         // function value at potential new point
	int maxind(0);      // Index of simplex node with maximal function value
	int minind(0);      // Index of simplex node with minimal function value
	double stepsize(1); // Measure for the size of the step
	double fstepsize(1);// Measure of size of a step in function values
	int iter(0);        // Number of simplex iterations
	int run;            // the model should be rerun after convergence


	for (run=0; run<2; run++) {
		for (k=1; k<nparameters+1; k++) simplex[k][k-1] *= 1.5;
		// for (k=1; k<nparameters+1; k++) simplex[k][k-1] += .05;
		iter = 0;
		while (1) {
			// Evaluate model at every simplex node and determine maximum and minimum
			maxind = minind = 0;
			for (k=0; k<nparameters+1; k++) {
				if (modified[k]) {
					fx[k] = model->neglpost(simplex[k], data );
					modified[k] = false;
				}
				// fx[k] = testfunction(simplex[k]);
#ifdef DEBUG_OPTIMIZER
				// std::cerr << simplex[k][0] << " " << simplex[k][1] << " " << simplex[k][2] << " " << simplex[k][3] << " " << "\n";
#endif
				if (fx[k]<fx[minind]) minind = k;
				if (fx[k]>fx[maxind]) maxind = k;
			}

			// Check Stoping criteria based on simplex and function values
			stepsize = 0;
			for (k=0; k<nparameters; k++)
				stepsize += (simplex[maxind][k]-simplex[minind][k])*(simplex[maxind][k]-simplex[minind][k]);
			// Simplex size
			if (stepsize<maxstep) {
#ifdef DEBUG_OPTIMIZER
				std::cerr << "Terminating optimization due to small simplex size (" << stepsize << ") after " << iter << " iterations\n";
				std::cout << "\n";
#endif
				break;
			}
			// function value differences
			if ((fstepsize=(fx[maxind]-fx[minind])) < maxfstep ) {
#ifdef DEBUG_OPTIMIZER
				std::cerr << "Terminating optimization due to small function value variation (" << fstepsize << ") after " << iter << " iterations\n";
				std::cout << "\n";
#endif
				break;
			}

#ifdef DEBUG_OPTIMIZER
			std::cout << iter << " " << fx[minind] << " " << stepsize << "\n";
#endif

			// Calculate the average of the non maximal nodes
			for (k=0; k<nparameters; k++) x[k] = 0;
			for (k=0; k<nparameters+1; k++) {
				if (k!=maxind)
					for (l=0; l<nparameters; l++)
						x[l] += simplex[k][l];
			}
			for (k=0; k<nparameters; k++) x[k] /= nparameters;

			// Determine the reflection of the worst point
			for (k=0; k<nparameters; k++) xx[k] = x[k] - (simplex[maxind][k]-x[k]);

			// Now check what to do
			ffx = model->neglpost(xx,data);
			// ffx = testfunction(xx);
			if (ffx<fx[minind]) {
				// The reflected point is better than the previous worst point ~> Expand
				for (k=0; k<nparameters; k++) simplex[maxind][k] = x[k] - 2*(simplex[maxind][k] - x[k]);
				modified[maxind] = true;
			} else if (ffx>fx[maxind]) {
				// The reflected point is even worse than it was before ~> Shrink
				for (k=0; k<nparameters+1; k++) {
					for (l=0; l<nparameters; l++)
						simplex[k][l] = simplex[minind][l] + 0.5 * (simplex[k][l] - simplex[minind][l]);
					modified[k] = true;
				}
			} else {
				// The reflected point is somewhere in between
				for (k=0; k<nparameters; k++) simplex[maxind][k] = xx[k];
				fx[maxind] = ffx;
			}

			// Also cancel if the number of iterations gets to large
			if (iter++ > maxiter) {
#ifdef DEBUG_OPTIMIZER
				std::cerr << "Terminating optimization due to large number of iterations (" << iter << "). Final stepsize: " << stepsize << "\n";
				std::cout << "\n";
#endif
				break;
			}
		}
	}

	// Evaluate model at every simplex node and determine minimum
	minind = 0;
	for (k=0; k<nparameters+1; k++) {
		if (modified[k]) {
			fx[k] = model->neglpost(simplex[k], data );
			modified[k] = false;
		}
		// fx[k] = testfunction(simplex[k]);
		if (fx[k]<fx[minind]) minind = k;
	}

	for (k=0; k<nparameters; k++)
	    start[k] = simplex[minind][k];

	/*
	// Perform some Gradient descent steps
	for (k=0; k<40; k++) {
		x = model->dnegllikeli ( start, data );
		dl = 0;
		for (l=0; l<nparameters; l++) {
			start[l] -= .1*x[l];
			if (fabs(x[l])>dl)
				dl = fabs(x[l]);
		}
		if (dl<1e-6)
			break;
	}
	*/

	return start;
}
