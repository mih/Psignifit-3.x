#include "integrate.h"
#include "errors.h"

std::vector<double> raw_grid (
		const PsiData *data,
		const PsiPsychometric *pmf,
		unsigned int prmindex,
		unsigned int gridsize
		)
{
	unsigned int i;
	double xmin,xmax,dx;
	std::vector<double> x (gridsize);
	parameter_range ( data, pmf, prmindex, &xmin, &xmax );

	dx = (xmax-xmin)/(gridsize-1);
	for ( i=0; i<gridsize; i++ )
		x[i] = xmin + i*dx;

	return x;
}

std::vector<double> cdf_grid (
		const PsiPrior *fitted_dist,
		double pmin,
		double pmax,
		unsigned int gridsize
		) {
	if ( pmin<0.05 || pmin>0.45 ) throw PsiError ( "pmin should be between 0.05 and 0.45" );
	if ( pmax<0.55 || pmax>0.95 ) throw PsiError ( "pmax should be between 0.55 and 0.95" );

	std::vector<double> p(gridsize);
	std::vector<double> x(gridsize);
	double dp ((pmax-pmin)/(gridsize-1));
	unsigned int i;

	for ( i=0; i<gridsize; i++ )
		p[i] = pmin + i*dp;

	x[0] = fitted_dist->ppf ( p[0], fitted_dist->mean() );

	for ( i=1; i<gridsize; i++ ) {
		x[i] = fitted_dist->ppf ( p[i], x[i-1] );
#ifdef DEBUG_INTEGRATE
		std::cerr << "ppf(" << p[i] << ") = " << x[i] << ",";
#endif
	}
#ifdef DEBUG_INTEGRATE
	std::cerr << "\n";
#endif

	return x;
}

void integrate_3d (
		const PsiPsychometric *pmf,
		const PsiData *data,
		const std::vector<double>& grid1,
		const std::vector<double>& grid2,
		const std::vector<double>& grid3,
		std::vector<double> *margin1,
		std::vector<double> *margin2,
		std::vector<double> *margin3
		)
{
	unsigned int i,j,k;
	double f;
	std::vector<double> prm ( 3 );

	// Make sure the margins start with only zeros
	for ( i=0; i<grid1.size(); i++ )
		(*margin1)[i] = 0;
	for ( j=0; j<grid2.size(); j++ )
		(*margin2)[j] = 0;
	for ( k=0; k<grid3.size(); k++ )
		(*margin3)[k] = 0;

	for ( i=0; i<grid1.size(); i++ ) {
		prm[0] = grid1[i];
		for ( j=0; j<grid2.size(); j++ ) {
			prm[1] = grid2[j];
			for ( k=0; k<grid3.size(); k++ ) {
				prm[2] = grid3[k];
				f = exp ( - pmf->neglpost ( prm, data ) );
				(*margin1)[i] += f;
				(*margin2)[j] += f;
				(*margin3)[k] += f;
			}
		}
	}

	return;
}

void integrate_4d (
		const PsiPsychometric *pmf,
		const PsiData *data,
		const std::vector<double>& grid1,
		const std::vector<double>& grid2,
		const std::vector<double>& grid3,
		const std::vector<double>& grid4,
		std::vector<double> *margin1,
		std::vector<double> *margin2,
		std::vector<double> *margin3,
		std::vector<double> *margin4
		)
{
	unsigned int i,j,k,l;
	double f;
	std::vector<double> prm ( 4 );

	// Make sure the margins start with only zeros
	for ( i=0; i<grid1.size(); i++ )
		(*margin1)[i] = 0;
	for ( j=0; j<grid2.size(); j++ )
		(*margin2)[j] = 0;
	for ( k=0; k<grid3.size(); k++ )
		(*margin3)[k] = 0;
	for ( l=0; l<grid4.size(); l++ )
		(*margin4)[l] = 0;

	for ( i=0; i<grid1.size(); i++ ) {
		prm[0] = grid1[i];
		for ( j=0; j<grid2.size(); j++ ) {
			prm[1] = grid2[j];
			for ( k=0; k<grid3.size(); k++ ) {
				prm[2] = grid3[k];
				for ( l=0; l<grid4.size(); l++ ) {
					prm[3] = grid4[l];
					f = exp ( - pmf->neglpost ( prm, data ) );
					(*margin1)[i] += f;
					(*margin2)[j] += f;
					(*margin3)[k] += f;
					(*margin4)[k] += f;
				}
			}
		}
	}

	return;
}

void integrate_grid (
		const PsiPsychometric *pmf,
		const PsiData *data,
		const std::vector< std::vector<double> >& grid,
		std::vector<double> *margin1,
		std::vector<double> *margin2,
		std::vector<double> *margin3,
		std::vector<double> *margin4
		) throw (PsiError)
{
	unsigned int nprm ( pmf->getNparams () );

	if ( grid.size() != nprm )
		throw PsiError ( "integrate_grid: number of grid dimensions does not match with number of parameters" );
	if ( nprm==4 && margin4==NULL )
		throw PsiError ( "integrate_grid: marginalizations for four parameters requested but only three output vectors provided" );

	if ( grid[0].size() != margin1->size() )
		throw PsiError ( "integrate_grid: grid[0].size() != margin1->size()" );
	if ( grid[1].size() != margin2->size() )
		throw PsiError ( "integrate_grid: grid[0].size() != margin1->size()" );
	if ( grid[2].size() != margin3->size() )
		throw PsiError ( "integrate_grid: grid[0].size() != margin1->size()" );
	if ( nprm==4 )
		if ( grid[3].size() != margin4->size() )
			throw PsiError ( "integrate_grid: grid[0].size() != margin1->size()" );

	if ( nprm==3 )
		integrate_3d ( pmf, data, grid[0], grid[1], grid[2], margin1, margin2, margin3 );
	else if ( nprm==4 )
		integrate_4d ( pmf, data, grid[0], grid[1], grid[2], grid[3], margin1, margin2, margin3, margin4 );
	else
		throw PsiError ( "integrate_grid: neither 3 nor 4 parameters to integrate this should not happen" );
}

std::vector<double> fit_posterior (
		const std::vector<double>& x,
		const std::vector<double>& fx,
		const std::vector<double>& start,
		unsigned int index
		)
{
	double (*error) ( const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

	std::vector<double> e ( 4 );
	std::vector<double> centroid ( start ), suggested1 ( 3 ), suggested2(3);
	if ( index != 0 ) {
		centroid[0] = sqrt ( centroid[0] );
		centroid[1] = sqrt ( centroid[1] );
	}
	unsigned int i, j, iter;
	double mmax(0);
	std::vector<double> margin ( fx );
	for ( i=0; i<margin.size(); i++ )
		if ( margin[i] > mmax )
			mmax = margin[i];
	for ( i=0; i<margin.size(); i++ )
		margin[i] = log ( margin[i]/mmax );

	centroid[2] = 0;
	for ( i=0; i<margin.size(); i++ ) centroid[2] += margin[i];
	centroid[2] /= double ( margin.size() );
	std::vector< std::vector<double> > simplex (4, centroid );
	double emin(1e5),emax(0);
	unsigned int imin(0), imax(0);
	double esuggested1, esuggested2, minmax;
	double me,se;

	const double alpha ( 1.), gamma ( 2.), beta ( 0.5 );

	switch ( index ) {
		case 0:
			error = &error_gauss;
			break;
		case 1:
			error = &error_gamma;
			break;
		case 2: case 3:
			error = &error_beta;
			break;
		default:
			throw PsiError ( "fit_posterior: invalid parameter index" );
	}

	for ( i=0; i<3; i++ )
		simplex[i][i] = ( simplex[i][i]==0 ? 0.1 : 1.3*simplex[i][i] );

	// Evaluate simplex
	for ( i=0; i<4; i++ )
		e[i] = (*error)(simplex[i], x, margin );

	// Start simplex here
	for ( iter=0; iter<200; iter++ ) {
		emin = 1e5; emax = 0;
		for ( i=0; i<4; i++ ) {
			if ( e[i] < emin ) { emin = e[i]; imin = i; }
			if ( e[i] > emax ) { emax = e[i]; imax = i; }
		}

		// Check for convergence
		se = 0.; me = 0.;
		for ( i=0; i<4; i++ ) me += e[i];
		me /= 4.;
		for ( i=0; i<4; i++ ) se += (e[i]-me)*(e[i]-me);
		se = sqrt(se/4.);
		if ( se < 1e-7 ) {
#ifdef DEBUG_INTEGRATE
			std::cerr << "Posterior fit converged after " << iter << " iterations\n";
			std::cerr << "    mean residual error:    " << me << "\n";
			std::cerr << "    minimum residual error: " << emin << "\n";
#endif
			break;
		}

		// Determine the centroid
		for ( j=0; j<3; j++ ) centroid[j] = 0;
		for ( i=0; i<4; i++ )
			if ( i!=imax )
				for ( j=0; j<3; j++ )
					centroid[j] += simplex[i][j];
		for ( j=0; j<3; j++ ) centroid[j] /= 3.;

		// Reflect on the centroid
		for ( j=0; j<3; j++ ) suggested1[j] = (1+alpha) * centroid[j] - alpha * simplex[imax][j];
		esuggested1 = (*error)( suggested1, x, margin );

		if ( esuggested1 > emin && esuggested1 < emax ) {
			for ( j=0; j<3; j++ )
				simplex[imax][j] = suggested1[j];
			e[imax] = esuggested1;
		} else if ( esuggested1 < emin ) {
			for ( j=0; j<3; j++ )
				suggested2[j] = gamma * suggested1[j] + (1-gamma) * centroid[j];
			esuggested2 = (*error)( suggested2, x, margin );
			if ( esuggested2 < esuggested1 ) {
				for ( j=0; j<3; j++ ) simplex[imax][j] = suggested2[j];
				e[imax] = esuggested2;
			} else {
				for ( j=0; j<3; j++ ) simplex[imax][j] = suggested1[j];
				e[imax] = esuggested1;
			}
		} else {
			minmax = ( esuggested1 < emax ? esuggested1 : emax );
			if ( esuggested1 < emax ) {
				for ( j=0; j<3; j++ ) simplex[imax][j] = suggested1[j];
				emax = esuggested1;
			}
			for ( j=0; j<3; j++ ) suggested2[j] = beta * simplex[imax][j] + (1-beta) * centroid[j];
			esuggested2 = (*error)( suggested2, x, margin );
			if ( esuggested2 > minmax ) e[imax] = emax;
			else {
				for ( j=0; j<3; j++ ) simplex[imax][j] = suggested2[j];
				e[imax] = esuggested2;
			}
		}
	}

	// Transform back
	if ( index != 0 ) {
		simplex[imin][0] *= simplex[imin][0];
		simplex[imin][1] *= simplex[imin][1];
	}

	return simplex[imin];
}

double error_gauss ( const std::vector<double>& prm, const std::vector<double>& x, const std::vector<double>& fx ) {
	double e (0), ee;
	unsigned int i;
	double a,b,Z;
	a = prm[0];
	b = prm[1];
	Z = prm[2];

	for ( i=0; i<x.size(); i++ ) {
		ee = -(x[i]-a)*(x[i]-a)/(2*b*b) - Z - fx[i];
		e += ee*ee;
	}

	return e;
}

double error_gamma ( const std::vector<double>& prm, const std::vector<double>& x, const std::vector<double>& fx ) {
	double e(0), ee;
	unsigned int i;
	double k,th,Z;
	k = prm[0]*prm[0];
	th = prm[1]*prm[1];
	Z = prm[2];

	for ( i=0; i<x.size(); i++ ) {
		ee = (k-1)*log(x[i]) - x[i]/th - Z - fx[i];
		e += ee*ee;
	}

	return e;
}

double error_beta ( const std::vector<double>& prm, const std::vector<double>& x, const std::vector<double>& fx ) {
	double e(0), ee;
	unsigned int i;
	double al,bt,Z;

	al = prm[0]*prm[0];
	bt = prm[1]*prm[1];
	Z  = prm[2];

	for ( i=0; i<x.size(); i++ ) {
		ee = (al-1)*log(x[i]) + (bt-1)*log(1-x[i]) - Z - fx[i];
		e += ee*ee;
	}

	return e;
}

std::vector<double> start_gauss ( const std::vector<double>& x, const std::vector<double>& fx ) {
	unsigned int i;
	std::vector<double> start ( 3 );

	start[0] = 0;
	for ( i=0; i<fx.size(); i++ ) {
		start[0] += x[i];
	}
	start[0] /= fx.size();

	start[1] = 0;
	for ( i=0; i<fx.size(); i++ ) {
		start[1] += (x[i]-start[0])*(x[i]-start[0]);
	}
	start[1] /= fx.size();
	start[1] = sqrt ( start[1] );

	return start;
}

std::vector<double> start_gamma ( const std::vector<double>& x, const std::vector<double>& fx ) {
	std::vector<double> start ( start_gauss ( x, fx ) );
	double theta, k;
	double m ( start[0] ), v ( start[1]*start[1] );
	theta = v / m;
	k     = v / (theta*theta);
	start[0] = k;
	start[1] = theta;
	return start;
}

std::vector<double> start_beta ( const std::vector<double>& x, const std::vector<double>& fx ) {
	std::vector<double> start ( start_gauss ( x, fx ) );
	double alpha, beta;
	double m ( start[0] ), v ( 2*start[1]*start[1] );
	alpha = (1-m)*m*m - v*m;
	alpha /= v;
	beta = alpha * (1-m)/m;
	start[0] = alpha;
	start[1] = beta;
#ifdef DEBUG_INTEGRATE
	std::cerr << "start: Beta(" << alpha << ", " << beta << ")\n";
#endif
	return start;
}

PsiIndependentPosterior independent_marginals (
		const PsiPsychometric *pmf,
		const PsiData *data,
		unsigned int nrefinements,
		unsigned int gridsize
		)
{
	unsigned int nprm ( pmf->getNparams() ), i, j;
	std::vector< std::vector<double> > grids ( nprm );
	std::vector< std::vector<double> > margin ( nprm, std::vector<double>(gridsize) );
	std::vector< std::vector<double> > distparams (nprm, std::vector<double>(3) );
	std::vector<PsiPrior*> fitted_posteriors (nprm);

	for ( i=0; i<nprm; i++ )
		grids[i] = raw_grid ( data, pmf, i, gridsize );

	integrate_grid ( pmf, data, grids, &(margin[0]), &(margin[1]), &(margin[2]), (nprm==4 ? &(margin[3]) : NULL) );

	for ( i=0; i<nprm; i++ ) {
		switch (i) {
			case 0:
				distparams[i] = fit_posterior ( grids[i], margin[i], start_gauss(grids[i], margin[i]), i );
				fitted_posteriors[i] = new GaussPrior ( distparams[i][0], distparams[i][1] );
				break;
			case 1:
				distparams[i] = fit_posterior ( grids[i], margin[i], start_gamma(grids[i], margin[i]), i );
				fitted_posteriors[i] = new GammaPrior ( distparams[i][0], distparams[i][1] );
				break;
			case 2: case 3:
				distparams[i] = fit_posterior ( grids[i], margin[i], start_beta(grids[i], margin[i]), i );
				fitted_posteriors[i] = new BetaPrior ( distparams[i][0], distparams[i][1] );
				break;
		}
	}

	for ( j=0; j<nrefinements; j++ ) {
		for ( i=0; i<nprm; i++ )
			grids[i] = cdf_grid ( fitted_posteriors[i], 0.1, 0.9, gridsize );
		integrate_grid ( pmf, data, grids, &(margin[0]), &(margin[1]), &(margin[2]), (nprm==4 ? &(margin[3]) : NULL ) );
		for ( i=0; i<nprm; i++ ) {
			delete fitted_posteriors[i];
			switch (i) {
				case 0:
					distparams[i] = fit_posterior ( grids[i], margin[i], start_gauss(grids[i], margin[i]), i );
					fitted_posteriors[i] = new GaussPrior ( distparams[i][0], distparams[i][1] );
					break;
				case 1:
					distparams[i] = fit_posterior ( grids[i], margin[i], start_gamma(grids[i], margin[i]), i );
					fitted_posteriors[i] = new GammaPrior ( distparams[i][0], distparams[i][1] );
					break;
				case 2: case 3:
					distparams[i] = fit_posterior ( grids[i], margin[i], start_beta(grids[i], margin[i]), i );
					fitted_posteriors[i] = new BetaPrior ( distparams[i][0], distparams[i][1] );
					break;
			}
		}
	}

	return PsiIndependentPosterior ( nprm, fitted_posteriors, grids, margin );
}

MCMCList sample_posterior (
		const PsiPsychometric *pmf,
		const PsiData *data,
		PsiIndependentPosterior& post,
		unsigned int nsamples
		)
{
	unsigned int nprm ( pmf->getNparams() ), i, j, k;
	unsigned int nproposals ( nsamples*10 );
	MCMCList finalsamples ( nsamples, nprm, data->getNblocks() );
	double q,p;
	double nduplicate ( 0 );
	PsiRandom rng;
	std::vector < PsiPrior* > posteriors ( nprm );

	std::vector< std::vector<double> > proposed ( nproposals, std::vector<double> (nprm) );
	std::vector<double> weights ( nproposals );
	std::vector<double> cum_probs ( nproposals );
	std::vector<double> rnumbers ( nproposals );

	for ( j=0; j<nprm; j++ )
		posteriors[j] = post.get_posterior (j)->clone();

	for ( i=0; i<nproposals; i++ ) {
		// Propose
		for ( j=0; j<nprm; j++ )
			proposed[i][j] = posteriors[j]->rand();
		// determine weight
		q = 1.;
        for ( j=0; j<nprm; j++ )
            q *= post.get_posterior (j)->pdf ( proposed[i][j] );
        p = exp ( - pmf->neglpost ( proposed[i], data ) );
        weights[i] = p/q;

		// Sort make a cumulative distribution vector for the weights
		if (i>0)
			cum_probs[i] = cum_probs[i-1] + weights[i];
		else
			cum_probs[0] = weights[0];

		// And generate random numbers
		rnumbers[i] = rng.rngcall();
	}

	for ( i=0; i<nproposals; i++ )
		cum_probs[i] /= cum_probs[nproposals-1];

	sort ( rnumbers.begin(), rnumbers.end() );

	// resampling
	i = j = 0;
    while (i<nsamples) {
        k = 0;
        while (rnumbers[i] <= cum_probs[j]) {
			finalsamples.setEst ( i, proposed[j], pmf->deviance ( proposed[j], data ) );
			nduplicate += k;
			k=1;
            i++;
            if (i>=nsamples)
                break;
		}
        j++;
	}

	finalsamples.set_accept_rate ( double(nduplicate)/nsamples );

	for ( i=0; i<nprm; i++ )
		delete posteriors[i];

	return finalsamples;
}

void sample_diagnostics (
		const PsiPsychometric *pmf,
		const PsiData *data,
		MCMCList *samples
		)
{
	unsigned int i,j,k, nprm ( pmf->getNparams() ), nblocks ( data->getNblocks() );
	std::vector<double> probs ( nblocks );
	std::vector<double> est ( nprm );
	PsiData *localdata = new PsiData ( data->getIntensities(), data->getNtrials(), data->getNcorrect(), data->getNalternatives() );
	std::vector<int> posterior_predictive ( nblocks );

	std::vector<double> reducedx ( data->getNblocks()-1 );
	std::vector<int> reducedk ( data->getNblocks()-1 );
	std::vector<int> reducedn ( data->getNblocks()-1 );
	std::vector< PsiData* > reduceddata ( data->getNblocks() );

	for ( i=0; i<nblocks; i++ ) {
		j = 0;
		for (k=0; k<nblocks; k++ ) {
			if ( j!=k ) {
				reducedx[j] = data->getIntensity(k);
				reducedk[j] = data->getNcorrect(k);
				reducedn[j] = data->getNtrials(k);
				j++;
			}
		}
		reduceddata[i] = new PsiData ( reducedx, reducedn, reducedk, data->getNalternatives() );
	}

	for ( i=0; i<samples->getNsamples(); i++ ) {
		for ( j=0; j<nprm; j++ )
			est[j] = samples->getEst ( i, j );

		for ( j=0; j<nblocks; j++ )
			probs[j] = pmf->evaluate ( data->getIntensity(j), est );
		newsample ( localdata, probs, &posterior_predictive );
		localdata->setNcorrect ( posterior_predictive );
		samples->setppData ( i, posterior_predictive, pmf->deviance ( est, localdata ) );

		probs = pmf->getDevianceResiduals ( est, data );
		samples->setRpd ( i, pmf->getRpd ( probs, est, data ) );
		samples->setRkd ( i, pmf->getRkd ( probs, data ) );

		probs = pmf->getDevianceResiduals ( est, localdata );
		samples->setppRpd ( i, pmf->getRpd ( probs, est, localdata ) );
		samples->setppRkd ( i, pmf->getRkd ( probs, localdata ) );

		// Store log posterior ratios for reduced data sets
		for ( j=0; j<nblocks; j++ )
			samples->setlogratio ( i, k, pmf->neglpost(est,data) - pmf->neglpost(est,reduceddata[j]) );
	}
}
