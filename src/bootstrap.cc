#include "bootstrap.h"

#ifdef DEBUG_BOOTSTRAP
#include <iostream>
#endif


BootstrapList parametricbootstrap ( int B, const PsiData * data, const PsiPsychometric* model, std::vector<double> cuts, bool BCa )
{
#ifdef DEBUG_BOOTSTRAP
	std::cerr << "Starting bootstrap\n Cuts size=" << cuts.size() << " "; std::cerr.flush();
#endif
	BootstrapList bootstrapsamples ( B, model->getNparams(), data->getNblocks(), cuts );
	int b,k,l,cut;
	std::vector< std::vector<double> > l_LF (B, std::vector<double>(cuts.size()));
	std::vector< std::vector<double> > u_i (B, std::vector<double>(cuts.size()));
	PsiOptimizer opt ( model, data );
	PsiData * localdataset = new PsiData ( data->getIntensities(),
			data->getNtrials(),
			data->getNcorrect(),
			data->getNalternatives() );

	std::vector<double> initialfit ( opt.optimize( model, data ) );
	std::vector<double> localfit   ( model->getNparams() );
	std::vector<int>    sample     ( data->getNblocks() );
	std::vector<double> initialthresholds ( cuts.size() );
	std::vector<double> devianceresiduals ( data->getNblocks() );
	double p,deviance;

	for (cut=0; cut<cuts.size(); cut++)
		initialthresholds[cut] = model->getThres(initialfit,cuts[cut]);

	for ( b=0; b<B; b++ ) {
		// Resampling
		for ( k=0; k<data->getNblocks(); k++ ) {
			p = model->evaluate( data->getIntensity(k), initialfit );
			sample[k] = 0;
			// TODO: use a better random number generator here
			for ( l=0; l<data->getNtrials(k); l++ )
				sample[k] += int( drand48() < p );
		}
		localdataset->setNcorrect ( sample );
		bootstrapsamples.setData ( b, sample );

		// Fit
		localfit = opt.optimize (model, localdataset );
#ifdef DEBUG_BOOTSTRAP
		if (false)
		std::cerr << localfit[0] << " " << localfit[1] << " " << localfit[2] << "\n";
#endif
		deviance = model->deviance ( localfit, localdataset );
		devianceresiduals = model->getDevianceResiduals ( localfit, localdataset );
		bootstrapsamples.setEst ( b, localfit, deviance );
		bootstrapsamples.setRpd ( b, model->getRpd( devianceresiduals, localfit, localdataset ) );
		bootstrapsamples.setRkd ( b, model->getRkd( devianceresiduals ) );

		// Store what we need for the BCa stuff
		for (cut=0; cut<cuts.size(); cut++) {
			l_LF[b][cut] = model->leastfavourable ( localfit, localdataset, cuts[cut] );
#ifdef DEBUG_BOOTSTRAP
			if (l_LF[b][cut] != l_LF[b][cut])
			std::cerr << "l_LF["<<b<<"]["<<cut<<"] = " << l_LF[b][cut] << "\n";
#endif
			u_i[b][cut]  = model->getThres(localfit,cuts[cut]);
#ifdef DEBUG_BOOTSTRAP
			if (l_LF[b][cut]!= l_LF[b][cut])
			std::cerr << "u_i["<<b<<"]["<<cut<<"] = " << u_i[b][cut] << "\n";
#endif
			bootstrapsamples.setThres(u_i[b][cut], b, cut);

			if (l_LF[b][cut] != l_LF[b][cut]) {
				// TODO: if l_LF is nan we don't take this sample
				// TODO: This is not the best solution but it works
				b--;
				continue;
			}
		}
	}

	// Calculate BCa constants
	double E_l,E_l3, var_l,w;
	for (cut=0; cut<cuts.size(); cut++) {
		// Get expectations
		E_l = 0; E_l3 = 0; w=0;
		for ( b=0; b<B; b++ ) {
			E_l += l_LF[b][cut];
			E_l3 += l_LF[b][cut]*l_LF[b][cut]*l_LF[b][cut];
			w += int( u_i[b][cut] < initialthresholds[cut] );
		}
		var_l = 0;
		for ( b=0; b<B; b++ ) {
			var_l += (l_LF[b][cut]-E_l)*(l_LF[b][cut]-E_l);
		}
		var_l /= B-1;
		var_l = sqrt(var_l); // note that this is not variance anymore but standard deviation

		// Store BCa stuff
#ifdef DEBUG_BOOTSTRAP
	std::cerr << " E_l3=" << E_l3 << "  var_l="<< var_l << "\n";
	std::cerr << " w = " << w/(B+1) << " real_w = " << invPhi(w/(B+1)) << "\n";
#endif
		bootstrapsamples.setBCa(cut, invPhi(w/(B+1)), E_l3/(6*var_l*var_l*var_l));
	}

	delete localdataset;

	return bootstrapsamples;
}

JackKnifeList jackknifedata ( const PsiData * data, const PsiPsychometric* model )
{
	PsiOptimizer *opt = new PsiOptimizer( model, data );
	std::vector<double> mlestimate ( opt->optimize( model, data ) );
	delete opt;
	JackKnifeList jackknife ( data->getNblocks(), model->getNparams(), model->deviance(mlestimate, data) );
	PsiData * localdata;

	std::vector<double> x ( data->getNblocks()-1 );
	std::vector<int> k ( data->getNblocks()-1 );
	std::vector<int> n ( data->getNblocks()-1 );
	int i,j,l,exclude(0);

	for ( i=0; i<data->getNblocks(); i++ ) {
		j=0;
		for (l=0; l<data->getNblocks(); l++) {
			if (l!=exclude) {
				x[j] = data->getIntensity(l);
				k[j] = data->getNcorrect(l);
				n[j] = data->getNtrials(l);
				j++;
			}
		}

		localdata = new PsiData ( x,n,k,data->getNalternatives() );
		opt       = new PsiOptimizer ( model, localdata );

		mlestimate = opt->optimize( model, localdata );
		jackknife.setEst ( i, mlestimate, model->deviance(mlestimate,localdata) );

		delete localdata;
		delete opt;
	}

	return jackknife;
}

