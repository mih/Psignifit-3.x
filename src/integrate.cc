#include "integrate.h"
#include "errors.h"
#include "linalg.h"

PsiIndependentPosterior::PsiIndependentPosterior ( unsigned int nprm,
				std::vector<PsiPrior*> posteriors,
				std::vector< std::vector<double> > x,
				std::vector< std::vector<double> > fx
				) : nparams (nprm), fitted_posteriors ( posteriors ), grids ( x ), margins ( fx ) {
	unsigned int i,j;
	std::vector<double> w;
	Matrix M ( grids[0].size(), 2 );

	for ( i=0; i<nparams; i++ ) {
		for ( j=0; j<grids[i].size(); j++ ) {
			M(j,0) = margins[i][j];
			M(j,1) = posteriors[i]->pdf ( grids[i][j] );
		}
		w = leastsq ( &M );
#ifdef DEBUG_INTEGRATE
		std::cerr << "w = " << w[0] << "\n";
#endif
		for ( j=0; j<margins[i].size(); j++ )
			margins[i][j] *= w[0];
	}
}

std::vector<double> lingrid ( double xmin, double xmax, unsigned int gridsize ) {
	unsigned int i;
	double dx;
	std::vector<double> x (gridsize);

	if ( xmin>xmax ) {
		dx = xmin;
		xmin = xmax;
		xmax = dx;
	}

	// std::cerr << "xmax = " << xmax << " xmin = " << xmin << "\n";
	dx = (xmax-xmin)/(gridsize-1);
	// std::cerr << "dx = " << dx << "\n";
	for ( i=0; i<gridsize; i++ )
		x[i] = xmin + i*dx;

	return x;
}

std::vector<double> raw_grid (
		const PsiData *data,
		const PsiPsychometric *pmf,
		unsigned int prmindex,
		unsigned int gridsize
		)
{
	/*
	double xmin,xmax;
	parameter_range ( data, pmf, prmindex, &xmin, &xmax );
	return lingrid ( xmin, xmax, gridsize );
	*/
	std::vector<double> incr ( pmf->getNparams() );
	std::vector<double> mapest (getstart ( pmf, data, gridsize, 3, 3, &incr ));
	int m ( (gridsize-1)/2 );

	return lingrid ( mapest[prmindex] - m*incr[prmindex], mapest[prmindex] + m*incr[prmindex], gridsize );
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

void normalize_margin ( std::vector<double> *margin ) {
	double m(0);
	unsigned int i;

	/*
	// Mean normalization
	for ( i=0; i<margin->size(); i++ )
		m += (*margin)[i];
	m /= double ( margin->size() );
	*/
	// max normalization
	for ( i=0; i<margin->size(); i++ )
		if ( (*margin)[i] > m )
			m = (*margin)[i];

	for ( i=0; i<margin->size(); i++ )
		(*margin)[i] /= m;
}

double sd ( const std::vector<double>& margin ) {
	unsigned int i;
	double s(0),m(0);

	for ( i=0; i<margin.size(); i++ )
		m += margin[i];
	m /= double ( margin.size() );
	for ( i=0; i<margin.size(); i++ )
		s += (margin[i]-m)*(margin[i]-m);
	s /= double ( margin.size() );
	return sqrt ( s );
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
	unsigned int i, j, iter, maxiter(500);
	double mmax(0);
	std::vector<double> margin ( fx );
	for ( i=0; i<margin.size(); i++ )
		if ( margin[i] > mmax )
			mmax = margin[i];
	for ( i=0; i<margin.size(); i++ )
		margin[i] = margin[i]/mmax;

	if ( index!=0 ) {
		centroid[2] = 0;
		for ( i=0; i<margin.size(); i++ ) centroid[2] += margin[i];
		centroid[2] /= double ( margin.size() );
	} else {
		maxiter = 10000;
	}
	std::vector< std::vector<double> > simplex (4, centroid );
	double emin(1e5),emax(0);
	unsigned int imin(0), imax(0);
	double esuggested1, esuggested2, minmax;
	double me,se;

#ifdef DEBUG_INTEGRATE
	std::cerr << "x = [" << x[0];
	for ( i=1; i<x.size(); i++ ) std::cerr << ", " << x[i];
	std::cerr << "]\nfx = [" << fx[0];
	for ( i=1; i<x.size(); i++ ) std::cerr << ", " << fx[i];
	std::cerr << "]\n";
#endif

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
	for ( iter=0; iter<maxiter; iter++ ) {
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
			std::cerr << "Posterior fit for parameter " << index << " converged after " << iter << " iterations\n";
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
#ifdef DEBUG_INTEGRATE
	if ( iter>=maxiter )
		std::cerr << "Maximum number of iterations reached for optimization of parameter " << index << ":\n iter = " << iter << "\n";
#endif

	// Transform back
	if ( index == 0 ) {
		simplex[imin][1] = fabs ( simplex[imin][1] );
	} else if ( index==1 ) {
		simplex[imin][0] *= simplex[imin][0];
		simplex[imin][1] *= simplex[imin][1];
	} else {
		simplex[imin][0] *= simplex[imin][0];
		simplex[imin][0] += 1;
		simplex[imin][1] *= simplex[imin][1];
		simplex[imin][1] += simplex[imin][0];
	}

	return simplex[imin];
}

double error_gauss ( const std::vector<double>& prm, const std::vector<double>& x, const std::vector<double>& fx ) {
	double e (0), ee;
	unsigned int i;
	double a,b,Z,m(0),z;
	a = prm[0];
	b = prm[1];
	Z = prm[2];

	// This should be in non log coordinates
	for ( i=0; i<x.size(); i++ ) {
		z = (x[i]-a)/b;
		ee = exp(-0.5*z*z + Z) - fx[i];
		e += ee*ee;
		m += x[i];
	}
	m /= double ( x.size() );
	e *= exp(0.001*(a-m)*(a-m));
	e *= exp(0.001*b*b);

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
		if ( x[i]>0 ) {
			ee = (k-1)*log(x[i]) - x[i]/th + Z - log(fx[i]);
			e += ee*ee;
		}
	}
	// e *= log(1e-5+k)/(1e-5+k);
	e += 0.001*th*th;

	return e;
}

double error_beta ( const std::vector<double>& prm, const std::vector<double>& x, const std::vector<double>& fx ) {
	double e(0), ee;
	unsigned int i;
	double al,bt,Z;

	// This parameterization might seem strange at first, but it makes certain constraints explicit:
	// 1. alpha>1
	// 2. beta>alpha>1
	al = 1+prm[0]*prm[0];
	bt = al+prm[1]*prm[1];
	Z  = prm[2];

	for ( i=0; i<x.size(); i++ ) {
		if ( x[i]>=0 && x[i]<=1 ) {
			ee = (al-1)*log(x[i]) + (bt-1)*log(1-x[i]) + Z - log(fx[i]);
			// ee = pow(x[i],al-1) * pow(1-x[i],bt-1) * Z - fx[i];
			e += ee*ee;
		}
	}
	// e -= log(1e-5+al)/(1e-5+al);
	// e -= log(1e-5+bt)/(1e-5+bt);
	// e += exp ( log ( 1+10*fabs(al-bt) ) );

	return e;
}

std::vector<double> start_gauss ( const std::vector<double>& x, const std::vector<double>& fx ) {
	unsigned int i;
	std::vector<double> start ( 3 );
	Matrix X ( x.size(), 4 );
	double sgsq, mu, k;

	// Perform linear regression in log coordinates here

	for ( i=0; i<x.size(); i++ ) {
		X(i,0) = 1;
		X(i,1) = x[i];
		X(i,2) = x[i]*x[i];
		X(i,3) = log(fx[i]);
	}
	start = leastsq ( &X );

	sgsq = -1./(2*start[2]);
	mu   = start[1] * sgsq;
	k    = start[0] + mu*mu/(2*sgsq);

	start[0] = mu;
	start[1] = sqrt(sgsq);
	start[2] = k;

#ifdef DEBUG_INTEGRATE
	std::cerr << "start: Gauss( " << start[0] << ", " << start[1] << ")\n";
#endif

	return start;
}

std::vector<double> start_gamma ( const std::vector<double>& x, const std::vector<double>& fx ) {
	std::vector<double> start ( 3 );
	double theta, k;
	double m ( 0 ), v ( 0 );
	unsigned int i;
	double s(0);

	for ( i=0; i<fx.size(); i++ ) {
		m += x[i]*fx[i];
		s += fx[i];
	}
	m /= s;
	for (i=0; i<fx.size(); i++ ) {
		v += (x[i]-m)*(x[i]-m) * fx[i];
	}
	v /= s;

	theta = v / m;
	k     = v / (theta*theta);
	start[0] = k;
	start[1] = theta;
#ifdef DEBUG_INTEGRATE
	std::cerr << "start: Gamma( " << start[0] << ", " << start[1] << ")\n";
#endif
	return start;
}

std::vector<double> start_beta ( const std::vector<double>& x, const std::vector<double>& fx ) {
	std::vector<double> start ( 3 );
	double alpha, beta;
	double m ( 0 ), v ( 0 );
	unsigned int i;
	double s(0);

	for ( i=0; i<fx.size(); i++ ) {
		m += x[i]*fx[i];
		s += fx[i];
	}
	m /= s;
	for (i=0; i<fx.size(); i++ ) {
		v += (x[i]-m)*(x[i]-m) * fx[i];
	}
	v /= s;
	v *= 2;

	alpha = (1-m)*m*m - v*m;
	alpha /= v;
	beta = alpha * (1-m)/m;
	start[0] = sqrt(alpha);
	start[1] = sqrt(beta-start[0]);
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
	unsigned int maxntrials ( 0 );
	double minp,minm,maxm,maxw,s;
	std::vector< std::vector<double> > grids ( nprm );
	std::vector< std::vector<double> > margin ( nprm, std::vector<double>(gridsize) );
	std::vector< std::vector<double> > distparams (nprm, std::vector<double>(3) );
	std::vector<PsiPrior*> fitted_posteriors (nprm);
	std::vector<double> marginsd (nprm);
	std::vector<double> maxprm (nprm);
	std::vector<double> minprm (nprm);

	for ( j=0; j<data->getNblocks(); j++ ) {
		if ( data->getNtrials(j) > maxntrials )
			maxntrials = data->getNtrials(j);
	}

	for ( i=0; i<nprm; i++ )
		grids[i] = raw_grid ( data, pmf, i, gridsize );
	minm = grids[0][0];
	maxm = grids[0].back();
	maxw = grids[1].back();
	minp = 1./maxntrials;

	integrate_grid ( pmf, data, grids, &(margin[0]), &(margin[1]), &(margin[2]), (nprm==4 ? &(margin[3]) : NULL) );
	for ( i=0; i<nprm; i++ ) {
		normalize_margin ( &(margin[i]) );
		marginsd[i] = sd ( margin[i] );
		minprm[i] = grids[i][0];
		maxprm[i] = grids[i].back();
		for ( j=0; j<margin[i].size(); j++ ) {
			if ( margin[i][j]<.1 )
				minprm[i] = grids[i][j];
			else
				break;
		}
		for ( j=margin[i].size()-1;j>=0; j-- ) {
			if ( margin[i][j]<.1 )
				maxprm[i] = grids[i][j];
			else
				break;
		}
	}

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
			grids[i] = lingrid ( minprm[i], maxprm[i], gridsize );
		/*
		for ( i=0; i<nprm; i++ )
			grids[i] = cdf_grid ( fitted_posteriors[i], 0.1, 0.9, gridsize );
		*/

		// Sanity check for grids
		/*
		if ( grids[0][gridsize-1] < minm )
			grids[0] = lingrid(minm,maxm, gridsize );
		else if ( grids[0][0] > maxm )
			grids[0] = lingrid(minm,maxm, gridsize );
		if ( grids[0][0] < minm )
			grids[0] = lingrid(minm,grids[0][gridsize-1], gridsize);
		if ( grids[0][gridsize-1] > maxm )
			grids[0] = lingrid(grids[0][0], maxm, gridsize);
		*/

		if ( grids[1][0] < 0 )
			grids[1] = lingrid(0,grids[1][gridsize-1], gridsize );
		/*
		if ( grids[1][0]!=grids[1][0] )
			grids[1] = raw_grid ( data, pmf, i, gridsize );
		if ( grids[1][gridsize-1] > maxw )
			grids[1] = lingrid ( grids[1][0], maxw, gridsize );
		*/

		/*
		if ( grids[2][gridsize-1] < minp )
			grids[2] = lingrid(grids[2][0], minp, gridsize);
		if ( grids[2][gridsize-1] > 1./pmf->getNparams() )
			grids[2] = lingrid(grids[2][0], 1./pmf->getNparams(), gridsize );
		*/
		for ( i=0; i<gridsize; i++ )
		if ( grids[2][0] < 0 )
			grids[2] = lingrid(0, grids[2][gridsize-1], gridsize );
		if ( grids[2][gridsize-1] > 1 )
			grids[2] = lingrid(grids[2][0], 1, gridsize );

		integrate_grid ( pmf, data, grids, &(margin[0]), &(margin[1]), &(margin[2]), (nprm==4 ? &(margin[3]) : NULL ) );
		for ( i=0; i<nprm; i++ ) {
			s = sd ( margin[i] );
			// only refit, if the new grid has more variance
			if ( s < marginsd[i] )
				continue;
			else {
				marginsd[i] = s;
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

		for ( i=0; i<nprm; i++ ) {
			normalize_margin ( &(margin[i]) );
			marginsd[i] = sd ( margin[i] );
			minprm[i] = -1e5;
			maxprm[i] = 1e5;
			for ( j=0; j<margin[i].size(); j++ ) {
				if ( margin[i][j]<.1 )
					minprm[i] = grids[i][j];
				else
					break;
			}
			for ( j=margin[i].size()-1;j>=0; j++ ) {
				if ( margin[i][j]<.1 )
					maxprm[i] = grids[i][j];
				else
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
	unsigned int nproposals ( nsamples*25 );
	MCMCList finalsamples ( nsamples, nprm, data->getNblocks() );
	double q,p;
	double nduplicate ( 0 );
	PsiRandom rng;
	std::vector < PsiPrior* > posteriors ( nprm );
	double H(0);

	std::vector< std::vector<double> > proposed ( nproposals, std::vector<double> (nprm) );
	std::vector<double> weights ( nproposals );
	std::vector<double> cum_probs ( nproposals );
	std::vector<double> rnumbers ( nsamples );

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
	}

	for ( i=0; i<nsamples; i++ ) {
		// And generate random numbers
		rnumbers[i] = rng.rngcall();
	}

	for ( i=0; i<nproposals; i++ )
		cum_probs[i] /= cum_probs[nproposals-1];

	H = - cum_probs[0] * log(cum_probs[0]);
	// Avoid zeros
	for ( i=0; i<nproposals-1; i++ )
		H -= (cum_probs[i+1]-cum_probs[i]) * log ( cum_probs[i+1]-cum_probs[i] );
	H /= log(nproposals);
	std::cerr << "H = " << H << "\n";

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
		if (j>nproposals) {
#ifdef DEBUG_INTEGRATE
			std::cerr << "What's going on here? i=" <<  i << ", cum_probs.max() = " << cum_probs[nproposals-1] << "\n";
#endif
			break;
		}
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
			samples->setlogratio ( i, j, pmf->neglpost(est,data) - pmf->neglpost(est,reduceddata[j]) );
	}
}
