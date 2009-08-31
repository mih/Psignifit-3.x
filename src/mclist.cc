#include "mclist.h"


/************************************************************
 * PsiMClist methods
 */

std::vector<double> PsiMClist::getEst ( int i ) const
{
	// Check that the call does not ask for something we don't have
	if ( i>=getNsamples() )
		throw BadIndexError();

	int k;
	std::vector<double> out ( getNparams() );

	for (k=0; k<getNparams(); k++)
		out[k] = mcestimates[k][i];

	return out;
}

double PsiMClist::getEst ( int i, int prm ) const
{
	// check that the call does not ask for something we don't have
	if ( i>=getNsamples() )
		throw BadIndexError();
	if ( prm>=getNparams() )
		throw BadIndexError();

	return mcestimates[prm][i];
}

void PsiMClist::setEst ( int i, const std::vector<double> est, double deviance )
{
	// Check that the call does not ask for something we can't do
	if ( i>=getNsamples() )
		throw BadIndexError();

	int k;
	for ( k=0; k<getNparams(); k++ )
		mcestimates[k][i] = est[k];
	deviances[i] = deviance;
}

double PsiMClist::getPercentile ( double p, int prm ) {
	if ( prm>=getNparams() )
		throw BadIndexError();
	if ( p>1 || p<0 )
		throw BadArgumentError();

	int position;
	sort( mcestimates[prm].begin(),mcestimates[prm].end() );

	position = getNsamples()*p;

	return mcestimates[prm][position];
}


void PsiMClist::setdeviance ( int i, double deviance ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	deviances[i] = deviance;
}

double PsiMClist::getdeviance ( int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return deviances[i];
}

double PsiMClist::getDeviancePercentile ( double p ) {
	if ( p<=0 || p>= 1 )
		throw BadArgumentError();

	int ind ( p*deviances.size() );

	sort( deviances.begin(), deviances.end() );

	return deviances[ind];
}

double PsiMClist::getMean ( unsigned int prm ) const {
	double m(0);
	int i,Nsamples(getNsamples());
	if ( prm>=getNparams() )
		throw BadIndexError();

	for (i=0; i<Nsamples; i++)
		m += getEst ( i, prm );

	m /= Nsamples;
	return m;
}

/************************************************************
 * BootstrapList methods
 */

void BootstrapList::setData ( unsigned int i, const std::vector<int> newdata )
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError();

	int k;
	for ( k=0; k<getNblocks(); k++ )
		data[i][k] = newdata[k];
}

std::vector<int> BootstrapList::getData ( unsigned int i ) const
{
	if ( i>=getNsamples() || i<0 )
		throw BadIndexError();

	return data[i];
}

double BootstrapList::getThres ( double p, unsigned int cut ) {
	if ( cut>=cuts.size() )
		throw BadIndexError();
	if ( p>1 || p<0 )
		throw BadArgumentError();

	int position;
	sort( thresholds[cut].begin(), thresholds[cut].end() );

	// Bias correction of p
	if (BCa)
		p = Phi(bias[cut] + (invPhi(p) + bias[cut])/(1-acceleration[cut]*(invPhi(p) + bias[cut])));

	position = int(getNsamples()*p);

	return thresholds[cut][position];
}

void BootstrapList::setThres ( double thres, unsigned int i, unsigned int cut )
{
	if ( i>=getNsamples() )
		throw BadIndexError();
	if (cut>=cuts.size() )
		throw BadIndexError();

	thresholds[cut][i] = thres;
}

double BootstrapList::getCut ( unsigned int i ) const
{
	if ( i>=cuts.size() || i<0 )
		throw BadIndexError();

	return cuts[i];
}

void BootstrapList::setRpd ( unsigned int i, double r_pd ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	Rpd[i] = r_pd;
}

double BootstrapList::getRpd ( unsigned int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return Rpd[i];
}

double BootstrapList::percRpd ( double p ) {
	if ( p<0 || p>1 )
		throw BadArgumentError();

	int index ( p*getNsamples() );

	sort ( Rpd.begin(), Rpd.end() );

	return Rpd[index];
}

void BootstrapList::setRkd ( unsigned int i, double r_kd ) {
	if ( i>=getNsamples() )
		throw BadIndexError();

	Rkd[i] = r_kd;
}

double BootstrapList::getRkd ( unsigned int i ) const {
	if ( i>=getNsamples() )
		throw BadIndexError();

	return Rkd[i];
}

double BootstrapList::percRkd ( double p ) {
	if ( p<0 || p>1 )
		throw BadIndexError();

	int index ( p*getNsamples() );

	sort ( Rkd.begin(), Rkd.end() );

	return Rkd[index];
}

/************************************************************
 * JackKnifeList methods
 */

bool JackKnifeList::influential ( unsigned int block, const std::vector<double>& ci_lower, const std::vector<double>& ci_upper ) const {
	int prm;
	double est;

	for ( prm=0; prm<getNparams(); prm++ ) {
		est = getEst(block,prm);
		if ( est<ci_lower[prm] || est>ci_upper[prm] )
			return true;   // This block is an influential observation for parameter prm
	}
	return false;
}

bool JackKnifeList::outlier ( unsigned int block ) const {
	if ( block>=getNblocks() )
		throw BadIndexError();

	if ( maxdeviance-getdeviance(block) > 6.63 )
		return true;
	return false;
}
