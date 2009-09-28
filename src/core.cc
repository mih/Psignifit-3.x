#include "core.h"

/************************************************************
 * abCore methods
 */

double abCore::dg ( double x, const std::vector<double>& prm, int i ) {
	switch (i) {
	case 0:
		return -1./prm[1];
		break;
	case 1:
		return -(x-prm[0])/(prm[1]*prm[1]);
		break;
	default:
		// If the parameter does not exist in the abCore the derivative with respect to it will always be 0
		return 0;
		break;
	}
}

double abCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) {
	if (i==j) {
		switch (i) {
		case 0:
			return 0;
			break;
		case 1:
			return (x-prm[0])/(prm[1]*prm[1]*prm[1]);
			break;
		default:
			// If the parameter does not exist in the abCore the derivative with respect to it will always be 0
			return 0;
			break;
		}
	} else if ((i==0 && j==1) || (i==1 && j==0)) {
	    return 1./(prm[1]*prm[1]);
	} else
		    // If the parameter does not exist in the abCore the derivative with respect to it will always be 0
		return 0;
}

double abCore::inv ( double y, const std::vector<double>& prm ) {
	return y*prm[1] + prm[0];
}

double abCore::dinv ( double y, const std::vector<double>& prm, int i ) {
	switch (i) {
	case 0:
		return 1;
		break;
	case 1:
		return y;
		break;
	default:
		return 0;
		break;
	}
}

std::vector<double> abCore::transform ( int nprm, double a, double b ) {
	std::vector<double> out ( nprm, 0 );
	out[1] = 1./b;
	out[0] = -a/b;
	return out;
}

/************************************************************
 * mwCore methods
 */
mwCore::mwCore ( int sigmoid, double al )
	: sigmtype(sigmoid), alpha(al), zshift(0) {
	switch (sigmoid) {
		case 1:
			// logistic
			zalpha = 2*log(1./alpha-1.);
			break;
		case 2:
			// gaussian
			zalpha = invPhi(1-alpha)-invPhi(alpha);
			break;
		case 3:
			// gumbel
			zalpha = log(-log(alpha))-log(-log(1.-alpha));
			zshift = log(-log(0.5));
			break;
		case 4:
			// cauchy
			zalpha = -2*tan(M_PI*(alpha-0.5));
			zshift = 0;
			break;
		default:
			throw NotImplementedError();
	}
}

double mwCore::g ( double x, const std::vector<double>& prm ) {
	return zalpha*(x-prm[0])/prm[1] + zshift;
}

double mwCore::dg ( double x, const std::vector<double>& prm, int i ) {
	switch (i) {
		case 0:
			return zalpha/prm[1];
			break;
		case 1:
			return -zalpha*(x-prm[0])/(prm[1]*prm[1]);
			break;
		default:
			// Irrelevant parameter ~> derivative is 0.
			return 0;
			break;
	}
}

double mwCore::ddg ( double x, const std::vector<double>& prm, int i, int j ) {
	if (i==j) {
		if (i==0)
			return 0;
		else if (i==1)
			return 2*zalpha*(x-prm[0])/(prm[1]*prm[1]*prm[1]);
		else
			return 0;
	} else if ( (i==0 && j==1) || (i==1 && j==0) )
		return zalpha/(prm[1]*prm[1]);
	else
		return 0;
}

double mwCore::inv ( double y, const std::vector<double>& prm ) {
	return prm[0] + prm[1]*(y-zshift)/zalpha;
}

double mwCore::dinv ( double p, const std::vector<double>& prm, int i ) {
	switch (i) {
		case 0:
			return 1;
			break;
		case 1:
			return (p-zshift)/zalpha;
			break;
		default:
			return 0;
			break;
	}
}

std::vector<double> mwCore::transform ( int nprm, double a, double b ) {
	std::vector<double> out ( nprm, 0 );
	out[1] = zalpha/b;
	out[0] = out[1]*(zshift-a)/zalpha;
	return out;
}

/************************************************************
 * logarithmicCore
 */

logCore::logCore ( const PsiData* data ) : scale(0) {
	int i;
	// we need this to scale starting values obtained from logistic regression so that they are correct "on average"
	for (i=0; i<data->getNblocks(); i++)
		scale += data->getIntensity(i)/log(data->getIntensity(i));
	scale /= data->getNblocks();
}

double logCore::dg ( double x, const std::vector<double>& prm, int i ) {
	switch (i) {
		case 0:
			return log(x);
			break;
		case 1:
			return 1;
			break;
		default:
			return 0;
			break;
	}
}

double logCore::dinv ( double y, const std::vector<double>& prm, int i ) {
	switch (i) {
		case 0:
			return exp((y-prm[1])/prm[0]) * (prm[1]-y)/(prm[0]*prm[0]);
			break;
		case 1:
			return -exp((y-prm[1])/prm[0])/prm[0];
			break;
		default:
			return 0;
			break;
	}
}

std::vector<double> logCore::transform ( int nprm, double a, double b ) {
	std::vector<double> prm ( nprm, 0 );
	prm[0] = b*scale;  // we scale the intercept so that it is correct "on average"
	prm[1] = a;
	return prm;
}
