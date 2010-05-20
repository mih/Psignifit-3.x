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

def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True ):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors)

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
    Rpd = np.zeros((nsamples))
    Rkd = np.zeros((nsamples))
    for row_index in xrange(nsamples):
        samples[row_index] = bs_list.getData(row_index)
        estimates[row_index] = bs_list.getEst(row_index)
        deviance[row_index] = bs_list.getdeviance(row_index)
        thres[row_index] = [bs_list.getThres_byPos(row_index, j) for j in xrange(ncuts)]
        Rpd[row_index] = bs_list.getRpd(row_index)
        Rkd[row_index] = bs_list.getRkd(row_index)

    acc = np.zeros((ncuts))
    bias = np.zeros((ncuts))
    for cut in xrange(ncuts):
        acc[cut] = bs_list.getAcc(cut)
        bias[cut] = bs_list.getBias(cut)

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

    return samples, estimates, deviance, thres, bias, acc, Rpd, Rkd, outliers, influential

def mcmc( data, start=None, nsamples=10000, nafc=2, sigmoid='logistic',
        core='mw0.1', priors=None, stepwidths=None ):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors)

    if start is not None:
        start = sfu.get_start(start, nparams)
    else:
        # use mapestimate
        opt = sfr.PsiOptimizer(pmf, dataset)
        start = opt.optimize(pmf, dataset)

    proposal = sfr.GaussRandom()
    sampler  = sfr.MetropolisHastings(pmf, dataset, proposal)
    sampler.setTheta(start)

    if stepwidths != None:
        if len(stepwidths) != nparams:
            raise sfu.PsignifitException("You specified \'"+str(len(start))+\
                    "\' stepwidth(s), but there are \'"+str(nparams)+ "\' parameters.")
        else:
            sampler.setstepsize(sfr.vector_double(stepwidths))

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

    return (estimates, deviance, posterior_predictive_data,
        posterior_predictive_deviances, posterior_predictive_Rpd,
        posterior_predictive_Rkd, logposterior_ratios)

def mapestimate ( data, nafc=2, sigmoid='logistic', core='ab', priors=None,
        cuts = None, start=None):

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, priors)

    cuts = sfu.get_cuts(cuts)

    opt = sfr.PsiOptimizer(pmf, dataset)
    estimate = opt.optimize(pmf, dataset, sfu.get_start(start, nparams) if start is not
            None else None)
    H = pmf.ddnegllikeli(estimate, dataset)
    thres = [pmf.getThres(estimate, c) for c in cuts]
    deviance = pmf.deviance(estimate, dataset)

    # convert to numpy stuff
    estimate = np.array(estimate)
    fisher = np.zeros((nparams,nparams))
    for (i,j) in ((i,j) for i in xrange(nparams) for j in xrange(nparams)):
        fisher[i,j] = sfr.doublep_value(H(i,j))
    thres = np.array(thres)
    deviance = np.array(deviance)

    return estimate, fisher, thres, deviance

def diagnostics(data, params, nafc=2, sigmoid='logistic', core='ab', cuts=None):
    # here we need to hack stuff, since data can be either 'real' data, or just
    # a list of intensities.
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

    dataset, pmf, nparams = sfu.make_dataset_and_pmf(data, nafc, sigmoid, core, None)
    cuts = sfu.get_cuts(cuts)
    # TODO length check params
    params = sfr.vector_double(params)
    predicted = np.array([pmf.evaluate(intensity, params) for intensity in
            dataset.getIntensities()])

    if intensities_only:
        return predicted
    else:
        deviance_residuals = pmf.getDevianceResiduals(params, dataset)
        deviance = pmf.deviance(params, dataset)
        thres = [pmf.getThres(params, cut) for cut in cuts]
        rpd = pmf.getRpd(deviance_residuals, params, dataset)
        rkd = pmf.getRkd(deviance_residuals, dataset)
        return predicted, deviance_residuals, deviance, thres, rpd, rkd
