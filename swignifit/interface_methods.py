#/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################
import numpy as np
import swignifit_raw as sfr
import swignifit.utility as sfu
import operator as op

def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True, gammaislambda=False ):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    cuts = sfu.get_cuts(cuts)
    ncuts = len(cuts)
    if start is not None:
        start = sfu.get_start(start, nparams)

    bs_list = sfr.bootstrap(nsamples, dataset, pmf, cuts, start, True, parametric)
    jk_list = sfr.jackknifedata(dataset, pmf)

    nblocks = dataset.getNblocks()

    # construct the massive tuple of return values
    samples = np.zeros((nsamples, nblocks), dtype=np.int32)
    estimates = np.zeros((nsamples, nparams))
    deviance = np.zeros((nsamples))
    thres = np.zeros((nsamples, ncuts))
    slope = np.zeros((nsamples, ncuts))
    Rpd = np.zeros((nsamples))
    Rkd = np.zeros((nsamples))
    for row_index in xrange(nsamples):
        samples[row_index] = bs_list.getData(row_index)
        estimates[row_index] = bs_list.getEst(row_index)
        deviance[row_index] = bs_list.getdeviance(row_index)
        thres[row_index] = [bs_list.getThres_byPos(row_index, j) for j in xrange(ncuts)]
        slope[row_index] = [bs_list.getSlope_byPos(row_index, j) for j in xrange(ncuts)]
        Rpd[row_index] = bs_list.getRpd(row_index)
        Rkd[row_index] = bs_list.getRkd(row_index)

    thacc = np.zeros((ncuts))
    thbias = np.zeros((ncuts))
    slacc = np.zeros((ncuts))
    slbias = np.zeros((ncuts))
    for cut in xrange(ncuts):
        thacc[cut] = bs_list.getAcc_t(cut)
        thbias[cut] = bs_list.getBias_t(cut)
        slacc[cut] = bs_list.getAcc_t(cut)
        slbias[cut] = bs_list.getBias_t(cut)

    ci_lower = sfr.vector_double(nparams)
    ci_upper = sfr.vector_double(nparams)

    for param in xrange(nparams):
        ci_lower[param] = bs_list.getPercentile(0.025, param)
        ci_upper[param] = bs_list.getPercentile(0.975, param)

    outliers = np.zeros((nblocks), dtype=np.bool)
    influential = np.zeros((nblocks))

    for block in xrange(nblocks):
        outliers[block] = jk_list.outlier(block)
        influential[block] = jk_list.influential(block, ci_lower, ci_upper)

    return samples, estimates, deviance, thres, thbias, thacc, slope, slbias, slacc, Rpd, Rkd, outliers, influential

def mcmc( data, start=None, nsamples=10000, nafc=2, sigmoid='logistic',
        core='mw0.1', priors=None, stepwidths=None, sampler="MetropolisHastings", gammaislambda=False):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    if start is not None:
        start = sfu.get_start(start, nparams)
    else:
        # use mapestimate
        opt = sfr.PsiOptimizer(pmf, dataset)
        start = opt.optimize(pmf, dataset)

    proposal = sfr.GaussRandom()
    if sampler not in sfu.sampler_dict.keys():
        raise sfu.PsignifitException("The sampler: " + sampler + " is not available.")
    else:
        sampler  = sfu.sampler_dict[sampler](pmf, dataset, proposal)
    sampler.setTheta(start)

    if stepwidths != None:
        stepwidths = np.array(stepwidths)
        if len(stepwidths.shape)==2:
            if isinstance ( sampler, sfr.GenericMetropolis ):
                sampler.findOptimalStepwidth ( sfu.make_pilotsample ( stepwidths ) )
            elif isinstance ( sampler, sfr.MetropolisHastings ):
                sampler.setStepSize ( sfr.vector_double( stepwidths.std(0) ) )
            else:
                raise sfu.PsignifitException("You provided a pilot sample but the selected sampler does not support pilot samples")
        elif len(stepwidths) != nparams:
            raise sfu.PsignifitException("You specified \'"+str(len(stepwidths))+\
                    "\' stepwidth(s), but there are \'"+str(nparams)+ "\' parameters.")
        else:
            if isinstance ( sampler, sfr.DefaultMCMC ):
                for i,p in enumerate(stepwidths):
                    p = sfu.get_prior(p)
                    sampler.set_proposal(i, p)
            else:
                sampler.setStepSize(sfr.vector_double(stepwidths))

    post = sampler.sample(nsamples)

    nblocks = dataset.getNblocks()

    estimates = np.zeros((nsamples, nparams))
    deviance = np.zeros(nsamples)
    posterior_predictive_data = np.zeros((nsamples, nblocks))
    posterior_predictive_deviances = np.zeros(nsamples)
    posterior_predictive_Rpd = np.zeros(nsamples)
    posterior_predictive_Rkd = np.zeros(nsamples)
    logposterior_ratios = np.zeros((nsamples, nblocks))

    for i in xrange(nsamples):
        for j in xrange(nparams):
            estimates[i, j] = post.getEst(i, j)
        deviance[i] = post.getdeviance(i)
        for j in xrange(nblocks):
            posterior_predictive_data[i, j] = post.getppData(i, j)
            logposterior_ratios[i,j] = post.getlogratio(i,j)
        posterior_predictive_deviances[i] = post.getppDeviance(i)
        posterior_predictive_Rpd[i] = post.getppRpd(i)
        posterior_predictive_Rkd[i] = post.getppRkd(i)

    accept_rate = post.get_accept_rate()

    return (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios, accept_rate)

def mapestimate ( data, nafc=2, sigmoid='logistic', core='ab', priors=None,
        cuts = None, start=None, gammaislambda=False):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda)

    cuts = sfu.get_cuts(cuts)

    opt = sfr.PsiOptimizer(pmf, dataset)
    estimate = opt.optimize(pmf, dataset, sfu.get_start(start, nparams) if start is not
            None else None)
    H = pmf.ddnegllikeli(estimate, dataset)
    thres = [pmf.getThres(estimate, c) for c in cuts]
    slope = [pmf.getSlope(estimate, th) for th in thres]
    deviance = pmf.deviance(estimate, dataset)

    # convert to numpy stuff
    estimate = np.array(estimate)
    fisher = np.zeros((nparams,nparams))
    for (i,j) in ((i,j) for i in xrange(nparams) for j in xrange(nparams)):
        fisher[i,j] = sfr.doublep_value(H(i,j))
    thres = np.array(thres)
    slope = np.array(slope)
    deviance = np.array(deviance)

    return estimate, fisher, thres, slope, deviance

def diagnostics(data, params, nafc=2, sigmoid='logistic', core='ab', cuts=None, gammaislambda=False):
    # here we need to hack stuff, since data can be either 'real' data, or just
    # a list of intensities, or just an empty sequence

    # in order to remain compatible with psipy we must check for an empty
    # sequence here, and return a specially crafted return value in that case.
    # sorry..
    if op.isSequenceType(data) and len(data) == 0:
        pmf, nparams =  sfu.make_pmf(sfr.PsiData([0],[0],[0],1), nafc, sigmoid, core, None, gammaislambda=gammaislambda )
        thres = np.array([pmf.getThres(params, cut) for cut in sfu.get_cuts(cuts)])
        slope = np.array([pmf.getSlope(params, th ) for th in thres])
        return np.array([]), np.array([]), 0.0, thres, np.nan, np.nan

    shape = np.shape(np.array(data))
    intensities_only = False
    if len(shape) == 1:
        # just intensities, make a dataset with k and n all zero
        k = n = [0] * shape[0]
        data  = [[xx,kk,nn] for xx,kk,nn in zip(data,k,n)]
        intensities_only = True
    else:
        # data is 'real', just do nothing
        pass

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, None, gammaislambda=gammaislambda)
    cuts = sfu.get_cuts(cuts)
    params = sfu.get_params(params, nparams)
    predicted = np.array([pmf.evaluate(intensity, params) for intensity in
            dataset.getIntensities()])

    if intensities_only:
        return predicted
    else:
        deviance_residuals = pmf.getDevianceResiduals(params, dataset)
        deviance = pmf.deviance(params, dataset)
        thres = np.array([pmf.getThres(params, cut) for cut in cuts])
        slope = np.array([pmf.getSlope(params, th ) for th in thres])
        rpd = pmf.getRpd(deviance_residuals, params, dataset)
        rkd = pmf.getRkd(deviance_residuals, dataset)
        return predicted, deviance_residuals, deviance, thres, slope, rpd, rkd

def asir ( data, nsamples=2000, nafc=2, sigmoid="logistic",
        core="mw0.1", priors=None, gammaislambda=False ):
    dataset, pmf, nparams = sfu.make_dataset_and_pmf ( data, nafc, sigmoid, core, priors, gammaislambda=gammaislambda )

    posterior = sfr.independent_marginals ( pmf, data, 1, 7 )
    samples   = sfr.sample_posterior ( pmf, data, posterior, nsamples )
    sfr.sample_diagnostics ( pmf, data, samples )

    out = {'mcestimates': np.array( [ [samples.getEst ( i, par ) for par in xrange ( nparams ) ] for i in xrange ( nsamples )]),
            'mcdeviance': np.array( [ samples.getdeviance ( i ) for i in xrange ( nparams ) ] ),
            'mcRpd':                    np.array ( [ samples.getRpd ( i ) for i in xrange ( nsamples ) ] ),
            'mcRkd':                    np.array ( [ samples.getRkd ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_data': np.array ( [ samples.getppData ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_deviance': np.array ( [ samples.getppDeviance ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_Rpd': np.array ( [ samples.getppRpd ( i ) for i in xrange ( nsamples ) ] ),
            'posterior_predictive_Rkd': np.array ( [ samples.getppRkd ( i ) for i in xrange ( nsamples ) ] ),
            'logposterior_ratios':      np.array ( [ samples.getlogratio ( i ) for i in xrange ( nsamples ) ] ),
            'duplicates':               samples.get_accept_rate (),
            'posterior_approximations_py': [posterior.get_posterior(i) for i in xrange ( nparams ) ],
            'posterior_approximations_str': None,
            'posterior_grids':          [ posterior.get_grid ( i ) for i in xrange ( nparams ) ],
            'posterior_margin':         [ posterior.get_margin ( i ) for i in xrange ( nparams ) ]
            }

    return out
