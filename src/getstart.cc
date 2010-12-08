#include "getstart.h"

std::vector<double> linspace ( double xmin, double xmax, unsigned int n ) {
	double dummy;
	if ( xmin>xmax ) { // make sure that xmin is really less than xmax
		dummy = xmin;
		xmin = xmax;
		xmax = dummy;
	}
	unsigned int i, n(out->size());
	double xstep ( (xmax-xmin)/(n-1) );
	std::vector<double> out (n)

	out[0] = xmin;
	for (i=1; i<n; i++) {
		out[i] = out[i-1] + xstep;
	}

	return out;
}

/************************************************** PsiGrid methods ********************************************/

PsiGrid::PsiGrid ( const std::vector<double>& xmin, const std::vector<double>& xmax, unsigned int gridsize ) :
	lower_bounds(xmin),
	upper_bounds(xmax)
{
	if ( lower_bounds.size() != upper_bounds.size() )
		throw Error ( "Upper and lower grid bounds are unequal" );
	unsigned int i;

	for ( i=0; i<lower_bounds.size(); i++ ) {
		grid1d.push_back ( linspace ( lower_bounds[i], upper_bounds[i], gridsize ) );
	}
}

PsiGrid& PsiGrid::shift ( const std::vector<double>& newposition ) const
{
	std::vector<double> xmin ( lower_bounds );
	std::vector<double> xmax ( upper_bounds );
	unsigned int i;
	double gridmid;

	for ( i=0; i<newposition.size (); i++ ) {
		gridmid = (xmax[i]-xmin[i])/2.;
		xmin[i] += newposition[i]-gridmid;
		xmax[i] += newposition[i]-gridmid;
	}

	return PsiGrid( xmin, xmax, get_gridsize() );
}

PsiGrid& PsiGrid::shrink ( const std::vector<double>& newposition ) const
{
	std::vector<double> xmin ( lower_bounds );
	std::vector<double> xmax ( upper_bounds );
	std::list< std::vector<double> >::iterator i_grid1d;
	unsigned int i;
	double xstep;

	for ( i=0, i_grid1d=grid1d.begin(); i<newposition.size(); i++,i_grid1d++() ) {
		xstep = i_grid1d->at(1)-i_grid1d->at(0);
		xmin[i] = newposition[i]-xstep;
		xmax[i] = newposition[i]+xstep;
	}

	return PsiGrid ( xmin, xmax, get_gridsize() );
}

PsiGrid& PsiGrid::subgrid ( void ) const
{
	std::vector<double> xmin ( lower_bounds.size()-1 );
	std::vector<double> xmax ( upper_bounds.size()-1 );
	unsigned int i;

	for ( i=0; i<xmin.size(); i++ ) {
		xmin[i] = lower_bounds[i+1];
		xmax[i] = upper_bounds[i+1];
	}

	return PsiGrid ( xmin, xmax, get_gridsize() );
}

/************************************************** PsiGrid functions ********************************************/

void makegridpoints (
		const PsiGrid& grid,
		std::vector<double> prm,
		unsigned int pos,
		std::list< std::vector<double> > *gridpoints
		)
{
	if ( grid.dimension() != prm.size()-pos ) {
		throw Error ( "grid and parameter vector don't match" );
	}
	std::vector<double> gridvector;

	if ( grid.empty() ) {
		// We are in the lowest level of iteration
		gridpoints->push_back ( prm );
		return;
	} else {
		// We have to loop over this level
		gridvector = grid.front();
		for ( i=0; i<gridvector->size(); i++ ) {
			prm[pos] = gridvector[i];
			makegridpoints ( grid.subgrid(), prm, pos+1, gridpoints );
		}
	}
}

void evalgridpoints (
		std::list< std::vector<double> > *grid,
		std::list< std::vector<double> > *bestprm,
		std::list< double > *L,
		const PsiData* data,
		const PsiPsychometric* pmf,
		unsigned int nbest,
		)
{
	std::list< std::vector<double> >::iterator griditer;
	double l,p;
	double a,b,lm,gm;
	std::vector<double> prm;
	PsiCore *core = pmf->getCore();
	unsigned int i,j;
	bool store;

	for ( griditer=grid->begin(); griditer!=grid->end(); griditer++ ) {
		// Transform parameters and get negative log posterior
		a = (*griditer)[0];
		b = (*griditer)[1];
		prm = core->transform ( pmf->getNparams(), a, b );
		prm[2] = (*griditer)[2];
		if ( pmf->getNparams() > 3 ) prm[3] = (*griditer)[3];
		l = pmf->neglpost ( prm, data );

		// Where does it belong?
		for ( iter_L=L->begin(), iter_prm=bestprm->begin() ; iter_L!=L->end(); iter_L++, iter_prm++ ) {
			if ( l==(*iter_L) ) {
				if ( (*iter_prm) == (*griditer) )
					store = false;
				else
					store = true;
				break;
			} else if ( l<(*iter_L) ) {
				store = true;
				break;
			} else {
				store = false;
			}
		}
		// insert the values if they are good enough
		if ( store ) {
			L->insert ( iter_L, l );
			bestprm->insert ( iter_prm, std::vector<double>(*griditer) );
		}

		// Reduce the size of the best parameters list
		while ( L->size() > nbest ) {
			L->pop_back();
			bestprm->pop_back();
		}
	}
}

void updategridpoints (
		const PsiGrid& grid,
		const std::list< std::vector<double> >& bestprm,
		std::list< std::vector<double> > *newgridpoints,
		std::list< PsiGrid > *newgrids
		)
{
	// modify grid size: If a point in bestprm is on the edge of the grid: make a new, larger grid
	// if a point is an interior point of the grid, shring the grid to the area around that grid

	std::list< std::vector<double> >::iterator iter_prm;
	std::vector<double> prm ( bestprm->front().size() );
	bool isedge (false);
	unsigned int i;
	PsiGrid newgrid;

	for ( iter_prm=bestprm->begin(); iter_prm!=bestprm->end(); iter_prm++ ) {
		// Check whether the current point is on the edge of the grid
		isedge = false;
		for ( i=0); i<iter_prm->size(); i++ ) {
			isedge += (*iter_prm)[i]==grid.get_lower(i);
			isedge += (*iter_prm)[i]==grid.get_upper(i);
		}

		if (isedge) {
			newgrid = grid.shift ( *iter_prm );
		} else {
			newgrid = grid.shrink ( *iter_prm)
		}
		makegridpoints ( newgrid, prm, 0, newgridpoints );
		newgrids->push_back ( newgrid );
	}
}

/*************************************** Range heuristics ****************************************/

void a_range ( const PsiData* data, double *xmin, double *xmax ) {
	double x;
	unsigned int i;
	*xmin = 1e5;
	*xmax = -1e5;

	// Heuristic:
	// a will be between lowest and highes stimulus level
	for ( i=0; i<data->getNbocks(); i++ ) {
		x = data->getIntensity ( i );
		if ( x<*xmin ) {
			*xmin = x;
		}
		if ( x>*xmax ) {
			*xmax = x;
		}
	}
}

void b_range ( const PsiData* data, double *xmin, double *xmax ) {
	double x,p(1),xx,pp(0),pc;
	std::vector<double> intensities ( data->getIntensities() );
	unsigned int i,j;
	*xmin = 1e5;
	*xmax = -1e5;

	// Heuristic:
	// b will be between rising over the whole range of stimulus levels and
	// rising between two adjacent stimulus levels

	// First determine step sizes
	for ( i=0; i<intensities.size(); i++ ) {
		for ( j=i; j<intensities.size(); j++ ) {
			d = fabs ( intensities[i] - intensities[j] );
			if (d==0) continue;

			if (d>*xmax) {
				*xmax = d;
			}

			if (d<*xmin) {
				*xmin = d;
			}
		}
	}

	// Is the psychometric function rising or falling overall
	for ( i=0; i<intensities.size(); i++ ) {
		pc = data->getPcorrect ( i );
		w
		if ( pc<p ) {
			p = pc;
			x = intensities[i];
		}
		if ( pc>pp ) {
			pp = pc;
			xx = intensities[i];
		}
	}
	if ( xx<x ) { // psychometric function is falling
		x = *xmin;
		*xmin = *xmax;
		*xmax = x;
	}

	// In any case, b should be scaled to be roughly equivalent to w
	x = 2*log(9.);
	*xmin /= x;
	*xmax /= x;
}

void lm_range ( const PsiData* data, double *xmin, double *xmax ) {
	double p,pmax(0);
	unsigned int i;

	// Heuristic:
	// lm goes from 0 to twice the distance between 1 and the highest response probability
	for ( i=0; i<data->getNblocks(); i++ ) {
		p = data->getPcorrect ( i );
		if ( p > pmax ) {
			pmax = p;
		}
	}
	*xmin = 0;
	*xmax = 2*(1-pmax);
	if (*xmax>1) *xmax=.99;
}

void gm_range ( const PsiData* data, double *xmin, double *xmax ) {
	double p, pmin(0);
	unsigned int i;

	// Heuristic:
	// gm goes from 0 to twice the lowsest response probability
	for ( i=0; i<data->getNblocks(); i++ ) {
		p = data->getPcorrect ( i );
		if ( p<pmin ) {
			pmin = p;
		}
	}
	*xmin = 0;
	*xmax = 2*pmin;
	if (*xmax>1) *xmax=.99;
}

void parameter_range ( const PsiData* data, unsigned int prmindex, double *xmin, double *xmax ) {
	// Call the correct initial range function for the parameter
	switch ( prmindex ) {
		case 0:
			a_range ( data, xmin, xmax );
			break;
		case 1:
			b_range ( data, xmin, xmax );
			break;
		case 2:
			lm_range ( data, xmin, xmax );
			break;
		case 3:
			gm_range ( data, xmin, xmax );
			break;
	}
}

std::vector<double> getstart (
		const PsiPsychometric* pmf,
		const PsiData* data,
		unsigned int gridsize,
		unsigned int nneighborhoods,
		unsigned int niterations )
{
	std::vector<double> xmin ( pmf->getNparams() );
	std::vector<double> xmax ( pmf->getNparams() );
	std::list< std::vector<double> > bestprm;
	std::list< double > L;

	// Set up the initial grid
	for ( i=0; i<pmf->getNparams(); i++ ) {
		parameter_range ( data, i, &(xmin[i]), &(xmax[i]) );
	}

	// Make an initial grid
	PsiGrid grid ( xmin, xmax, gridsize );
	std::list< PsiGrid > newgrids;
	std::list< PsiGrid >::iterator i_grid;
	newgrids.push_back ( grid );

	// Perform first evaluation on the grid
	std::list< std::vector<double> > gridpoints;
	makegridpoints ( grid, xmin, 0, &gridpoints );
	evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, nneighborhoods );

	// potentially more evaluations
	for ( i=0; i<niterations; i++ ) {
		gridpoints = std::list< std::vector<double> > ();
		while ( newgrids.size() > nneighborhoods ) newgrids.pop_front ();

		for ( i_grid=newgrids.begin(); i_grid!=newgrids.end(); i_grid++ )
			updategridpoints ( *i_grid, bestprm, &gridpoints, &newgrids );

		evalgridpoints ( gridpoints, &bestprm, &L, data, pmf, nneighborhoods );
	}


	// Now transform the best parameter to the suitable format
	PsiCore *core = pmf->getCore();
	std::vector<double> out = core->transform ( pmf->getNparams(), bestprm.front()[0], bestprm.front()[1] );
	out[2] = bestprm.front()[2];
	if ( pmf->getNparams() > 3 ) out[3] = bestprm.front()[3];

	return out;
}
