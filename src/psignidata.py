#!/usr/bin/env python

import sys,os
import numpy as N
import pylab as p
from scipy import stats
import _psipy

import psigniplot as pp
import pygibbsit

from psignierrors import NosamplesError

__all__ = ["BootstrapInference","BayesInference"]
__doc__ = """
This module contains data objects to store psychophysical data and perform inference on them. Two general approaches
have been suggested to fit psychometric functions

1. *Constrained maximum likelihood (Wichmann & Hill, 2001a,b)* This procedure starts by fitting a psychometric function
   to the data and then performs parametric bootstrap to obtain confidence intervals for parameters, associated
   thresholds,... This approach is implemented by the :BootstrapInference: class.
2. *Baysian Inference (Kuss et al., 2005)* This procedure starts with a number of prior distributions for each of
   the models parameters and then uses Bayes rule to derive the posterior distribution from the data. As the
   posterior distribution can only be determined up to a normalization constant, inference on the posterior distribution
   has to be based on samples. These samples are typically obtained using Markoc Chain Monte Carlo (MCMC). This
   approach is implemented in the :BayesInference: class.

The module also defines a :PsiInference: base class.
"""
warnred = [.7,0,0]

# Helper function to create properties with one function
def Property(function):
    keys = 'fget', 'fset', 'fdel'
    func_locals = {'doc':function.__doc__}
    def probeFunc(frame, event, arg):
        if event == 'return':
            locals = frame.f_locals
            func_locals.update(dict((k,locals.get(k)) for k in keys))
            sys.settrace(None)
        return probeFunc
    sys.settrace(probeFunc)
    function()
    return property(**func_locals)

##############################################################################################################################
class PsiInference ( object ):
    def __init__ ( self ):
        """This is just a dummy function"""
        self.data = None
        self.model = {
                "sigmoid":  "logistic",
                "core":     "ab",
                "priors":   None,
                "nafc":     2
                }
        self.estimate          = None
        self.deviance          = None
        self.devianceresiduals = None
        self.Rpd               = None
        self.Rkd               = None
        self.__outl            = None
        self.__infl            = None

    def evaluate ( self, x, prm=None ):
        """Evaluate the psychometric function model at positions given by x"""
        if prm==None:
            prm = self.estimate

        return N.array( _psipy.diagnostics ( x, prm, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

    def getThres ( self, cut=0.5 ):
        """Get thresholds at cut"""
        if self.data == None:
            raise NotImplementedError
        return float(_psipy.diagnostics ( self.data, self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3])

    def __repr__ ( self ):
        return "< PsiInference object >"

    desc = property ( fget=lambda self: "sigmoid: %(sigmoid)s\ncore: %(core)s\nnAFC: %(nafc)d" % self.model,
            doc="A short description of the employed model")
    outl = property ( fget=lambda self: self.__outl, doc="A boolean array indicating whether or not a block was an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean array indicating whether or not a block was an influential observation" )

##############################################################################################################################
class BootstrapInference ( PsiInference ):
    def __init__ ( self, data, sample=False, cuts=(.25,.5,.75), conf=(.025,.975), **kwargs ):
        """Set up an object of bootstrapped data

        :Parameters:
            *data* :
                an array or a list of lists containing stimulus intensities in the
                first column, number of correct responses (nAFC) or number of YES-
                responses in the second column, and number of trials in the third
                column. Each row should correspond to one experimental block. In
                addition, the sequence of the rows is taken as the sequence of
                data aquisition.
            *sample* :
                if sample is True, bootstrap samples are drawn. If sample is an
                integer, it gives the number of samples that are drawn
            *sigmoid* :
                shape of the sigmoid function. Valid choices are
                    - 'logistic'   [Default]
                    - 'gauss'
                    - 'gumbel_l'
                    - 'gumbel_r'
            *core* :
                term inside the sigmoid function. Valid choices are
                    - 'ab'         (x-a)/b        [Default]
                    - 'mw%g'       midpoint and width
                    - 'linear'     a+b*x
                    - 'log'        a+b*log(x)
            *priors* :
                a list of prior names. Valid choices are
                    - 'Uniform(%g,%g)'   Uniform distribution on an interval
                    - 'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                    - 'Beta(%g,%g)'      Beta distribution
                    - 'Gamma(%g,%g)'     Gamma distribution
                If no valid prior is selected, the parameter remains unconstrained.
                Alternatively, priors can be given as a dictionary that only specifies
                priors for those parameters you want to set in that case you can use
                'a','b','m','w','guess','gamma','lapse','lambda' as keys.
            *nafc* :
                number of response alternatives. If nafc==1, this indicates a Yes/No
                task
            *cuts* :
                performance values that should be considered 'thresholds'. This means that a
                'cut' of 0.5 corresponds to an expected performance of roughly 75%% correct in
                a 2AFC task.
            *conf* :
                limits of confidence intervals. The default gives 95%% confidence intervals.
                Any other sequence can be used alternatively. In addition, conf can be 'v1.0'
                to give the default values of the classical psignifit version (i.e. .023,.159,.841,.977,
                corresponding to -2,-1,1,2 standard deviations for a gaussian).
        """
        # Call the base constructor
        PsiInference.__init__(self)
        self.__nsamples = 0

        # Store basic data
        self.data = N.array(data)
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "ab"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }

        if self.model["core"][:2] == "mw":
            self.parnames = ["m","w"]
        else:
            self.parnames = ["a","b"]
        self.parnames.append("lambda")
        if self.model["nafc"]<2:
            self.parnames.append("guess")

        self.cuts = cuts
        if conf=="v1.0":
            self.conf = (0.023, 0.159, 0.841, 0.977)
        else:
            self.conf = conf

        # Store point estimates
        self.estimate,self.thres,self.deviance = _psipy.mapestimate(self.data,cuts=self.cuts,**self.model)
        self.predicted,self.devianceresiduals,self.deviance,thres,self.Rpd,self.Rkd = _psipy.diagnostics(self.data,self.estimate)

        # The _psipy arrays are not numpy arrays
        self.estimate          = N.array(self.estimate)
        self.predicted         = N.array(self.predicted)
        self.devianceresiduals = N.array(self.devianceresiduals)

        # Bootstrap parameters are empty first
        self.__bdata     = None
        self.__bestimate = None
        self.__bdeviance = None
        self.__bRpd      = None
        self.__bRkd      = None
        self.__outl      = None
        self.__infl      = None
        self.__bthres    = None
        self.__th_bias   = None
        self.__th_acc    = None

        # If we want direct sampling this is done here
        if sample:
            if isinstance(sample,bool):
                self.sample ()
            elif isinstance(sample,int):
                if sample>0:
                    self.sample (sample)
                else:
                    raise ValueError, "Negative number of bootstrap samples selected"
        else:
            self.__nsamples = 0

    def sample ( self, Nsamples=2000 ):
        """Draw bootstrap samples

        :Parameters:
            *Nsamples* :
                number of bootstrapsamples to be drawn
        """
        self.__nsamples = Nsamples
        self.__bdata,self.__bestimate,self.__bdeviance,self.__bthres,self.__th_bias,self.__th_acc,\
                self.__bRkd,self.__bRpd,self.__outl,self.__infl = _psipy.bootstrap(self.data,self.estimate,Nsamples,cuts=self.cuts,**self.model)

        # Cast sampled data to numpy arrays
        self.__bdata = N.array(self.__bdata)
        self.__bestimate = N.array(self.__bestimate)
        self.__bdeviance = N.array(self.__bdeviance)
        self.__bthres    = N.array(self.__bthres)
        self.__th_bias   = N.array(self.__th_bias)
        self.__th_acc    = N.array(self.__th_acc)
        self.__bRkd      = N.array(self.__bRkd)
        self.__bRpd      = N.array(self.__bRpd)
        self.__outl      = N.array(self.__outl,dtype=bool)
        self.__infl      = N.array(self.__infl,dtype=bool)

    def getCI ( self, cut, conf=None ):
        """Determine the confidence interval of a cut

        :Parameters:
            *cut* :
                index(!) of the cut of interest
            *conf* :
                is currently ignored
        """
        bias = self.__th_bias[cut]
        acc  = self.__th_bias[cut]

        vals = []
        for pp in self.conf:
            vals.append(stats.norm.cdf( bias + ( stats.norm.ppf(pp) + bias ) / (1-acc*(stats.norm.ppf(pp) + bias )) ))

        return p.prctile ( self.__bthres[:,cut], 100*N.array(vals) )

    def __repr__ ( self ):
        return "< BootstrapInference object with %d blocks and %d samples >" % ( self.data.shape[0], self.nsamples )

    outl = property ( fget=lambda self: self.__outl, doc="A boolean vector indicating whether a block should be considered an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean vector indicating whether a block should be considered an influential observation" )
    mcestimates = property ( fget=lambda self: self.__bestimate, doc="An array of bootstrap estimates of the fitted paramters" )
    mcdeviance = property ( fget=lambda self: self.__bdeviance, doc="A vector of bootstrapped deviances" )
    mcRpd = property ( fget=lambda self: self.__bRpd, doc="A vector of correlations between model prections and deviance residuals in all bootstrap samples" )
    mcRkd = property ( fget=lambda self: self.__bRkd, doc="A vector of correlations between block index and deviance residuals in all bootstrap samples" )
    @Property
    def nsamples ():
        """number of bootstrap samples (setting this attribute results in resampling!!!)"""
        def fget ( self ):
            return self.__nsamples
        def fset ( self, n ):
            self.__nsamples = n
            self.sample ( self.__nsamples )

##############################################################################################################################
class BayesInference ( PsiInference ):
    def __init__ ( self, data, sample=True, cuts=(.25,.5,.75), conf=(.025,.975), automatic=True, resample=False, **kwargs ):
        """Bayesian Inference for psychometric functions using MCMC

        :Parameters:
            *data* :
                an array or a list of lists containing stimulus intensities in the
                first column, number of correct responses (nAFC) or number of YES-
                responses in the second column, and number of trials in the third
                column. Each row should correspond to one experimental block. In
                addition, the sequence of the rows is taken as the sequence of
                data aquisition.
            *sample* :
                if sample is True, bootstrap samples are drawn. If sample is an
                integer, it gives the number of samples that are drawn
            *sigmoid* :
                shape of the sigmoid function. Valid choices are
                    - 'logistic'   [Default]
                    - 'gauss'
                    - 'gumbel_l'
                    - 'gumbel_r'
            *core* :
                term inside the sigmoid function. Valid choices are
                    - 'ab'         (x-a)/b        [Default]
                    - 'mw%g'       midpoint and width
                    - 'linear'     a+b*x
                    - 'log'        a+b*log(x)
            *priors* :
                a list of prior names. Valid choices are
                    - 'Uniform(%g,%g)'   Uniform distribution on an interval
                    - 'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                    - 'Beta(%g,%g)'      Beta distribution
                    - 'Gamma(%g,%g)'     Gamma distribution
                If no valid prior is selected, the parameter remains unconstrained.
                Alternatively, priors can be given as a dictionary that only specifies
                priors for those parameters you want to set in that case you can use
                'a','b','m','w','guess','gamma','lapse','lambda' as keys.
            *nafc* :
                number of response alternatives. If nafc==1, this indicates a Yes/No
                task
            *cuts* :
                performance values that should be considered 'thresholds'. This means that a
                'cut' of 0.5 corresponds to an expected performance of roughly 75%% correct in
                a 2AFC task.
            *conf* :
                limits of confidence intervals. The default gives 95%% confidence intervals.
                Any other sequence can be used alternatively. In addition, conf can be 'v1.0'
                to give the default values of the classical psignifit version (i.e. .023,.159,.841,.977,
                corresponding to -2,-1,1,2 standard deviations for a gaussian).
            *automatic* :
                do everything automatically
            *resample* :
                if a chain is considered "bad" in terms of convergence should it
                automatically be resampled?
        """
        PsiInference.__init__(self)

        # Store basic data
        self.data = N.array(data)
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "mw0.1"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }
        self.resample = resample

        if self.model["core"][:2] == "mw":
            self.parnames = ["m","w"]
        else:
            self.parnames = ["a","b"]
        self.parnames.append("lambda")
        if self.model["nafc"]<2:
            self.parnames.append("guess")

        self.mapestimate,thres,self.mapdeviance = _psipy.mapestimate(self.data,**self.model)

        if cuts is None:
            self.cuts = (.25,.5,.75)
        else:
            self.cuts = cuts
        if isinstance(cuts,float):
            self.Ncuts = 1
        else:
            self.Ncuts = len(self.cuts)

        self.Rpd,self.Rkd = _psipy.diagnostics ( self.data, self.mapestimate, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[4:]

        self.__meanestimate = None
        self.__meandeviance = None

        self.__mcmc_chains    = []
        self.__mcmc_deviances = []

        self.__pRpd   = None
        self.__pRkd   = None
        self.__pthres = None

        self.conf = conf

        self.burnin = 0
        self.thin   = 1
        self.nsamples = None

        self._steps = (.4,4,.01)

        if automatic:
            self.__determineoptimalsampling ()
            sample = True

        if sample:
            self.sample()

    def sample ( self, Nsamples=None, start=None ):
        """Draw samples from the posterior distribution using MCMC

        :Parameters:
            *Nsamples* :
                number of samples that should be drawn from the posterior. If Nsamples is
                None, an optimal number of samples is tried to obtain.
            *start* :
                starting value of the chain. If this is None, the chain starts
                at the MAP estimate. However, if you sample multiple chains, you
                might want to start from different (overdispersed) starting values
                to diagnose convergence
        """
        if isinstance (Nsamples,int):
            self.nsamples = Nsamples
        elif Nsamples is None:
            Nsamples = self.nsamples+self.burnin
        else:
            Nsamples = 10000

        if start is None:
            start = self.mapestimate
        # TODO: stepwidths are not good for all situations
        # TODO: it would be better to have the sampler itself set the stepwidths based on asymptotic properties of the mapestimate
        stepwidths = (0.4,4,1e-2)
        chain,deviance = _psipy.mcmc ( self.data, start, Nsamples, stepwidths=self._steps, **self.model )
        # print N.cov(N.array(chain[self.burnin::self.thin]).T)
        self.__mcmc_chains.append(N.array(chain))
        self.__mcmc_deviances.append(N.array(deviance))

    def __repr__ ( self ):
        return "< BayesInference object with %d blocks and %d mcmc chains of %d samples each >" % (self.data.shape[0],len(self.__mcmc_chains), self.nsamples)

    ############################################
    # Setters and getters
    def getsamples ( self, chain=None, raw=False ):
        """Get sampes from the posterior

        :Parameters:
            *chain* :
                if chain is None, samples are aggregated over all chains
                sampled so far. If chain is an integer only data from the
                chain indexed by this number are returned

        :Output:
            an array of nsamplesXnparams samples from the posterior
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_chains:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance (chain,int):
            # Get a single chain
            if raw:
                return self.__mcmc_chains[chain]
            else:
                return self.__mcmc_chains[chain][self.burnin::self.thin]
        else:
            raise IndexError, "chain should be either None or an integer"

    def getmcdeviance ( self, chain=None ):
        """Get samples from the posterior distribution of deviances

        :Parameters:
            *chain* :
                if chain is None, the samples are combined across all chains
                sampled so far. If chain is an integer, it is interpreted as
                the index of the chain to be returned

        :Output:
            an array of samples from the posterior distribution of deviances. This
            array respects the burnin and thin settings.
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_deviances:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance ( chain, int ):
            return self.__mcmc_deviances[chain][self.burnin::self.thin]
        else:
            raise ValueError, "chain should be either None or an integer"

    def getCI ( self, cut=None, conf=(.025,0.5,.975), param="thres" ):
        """Get a posterior credibility interval for a particular parameter

        :Parameters:
            *conf* :
                percentiles that should be returned
            *param* :
                parameter of interest. Currently, only thres/threshold
                and Rkd,Rpd,deviance are defined.
        """
        if param[:5]=="thres":
            # We have to handle thresholds separately because there could be multiple cuts.
            mcthres = self.mcthres
            out = []
            if cut==None:
                for k in xrange(self.Ncuts):
                    out.append(p.prctile ( mcthres[:,k], 100*N.array(conf) ))
                return N.array(out)
            else:
                return p.prctile ( mcthres[:,cut], 100*N.array(conf) )
        else:
            if param=="Rkd":
                vals = self.mcRkd
            elif param=="Rpd":
                vals = self.mcRpd
            elif param=="deviance":
                vals = self.mcdeviance
            else:
                raise NotImplementedError
            return p.prctile ( vals, 100*N.array(conf) )

    ############################################
    # Plotting routines
    def drawposteriorexamples ( self, ax=None, Nsamples=20 ):
        """plots the mean estimate of the psychometric function and a number of samples from the posterior

        :Parameters:
            *ax* :
                axes object in which to draw the plot. If this is None,
                a new axes object is created.
            *Nsamples* :
                number of psychometric functions that should be drawn
                from the posterior
        """
        if ax is None:
            ax = p.axes()

        # Plot the psychometric function
        xmin = self.data[:,0].min()
        xmax = self.data[:,0].max()
        x = N.mgrid[xmin:xmax:100j]

        # Now we sample Nsamples psychometric functions from all the chains we have
        samples = self.getsamples()
        deviances = self.getmcdeviance()
        # Scale deviance to 0,1
        deviances -= deviances.min()
        deviances /= deviances.max()
        deviances = N.clip(.4+4*deviances,0,1)
        for k in xrange(Nsamples):
            i = N.random.randint(samples.shape[0])
            psi = N.array(_psipy.diagnostics ( x, samples[i,:], sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
            # Lines are colored according to their deviance
            ax.plot(x,psi,color=[deviances[i]]*2+[1])

    ############################################
    # Convergence diagnostics
    def geweke ( self, parameter=0, nsegments=10 ):
        """Geweke test for stationarity of a chain.

        The Geweke test first transforms all samples to mean 0 and standard deviation 1.
        In a second step, it calculates the sample average in a number of segments and
        checks whether these subaverages differ significantly from 0.

        :Parameters:
            *parameter* :
                parameter of interest
            *nsegments* :
                number of subaverages to be calculated

        :Output:
            a boolean value indicating whether the chain is "good" or "bad"
        """
        z = N.zeros ( (nsegments, self.nchains), 'd' )
        for k in xrange ( self.nchains ):
            samples = self.getsamples ( k ) [:,parameter]
            w = len(samples)/nsegments
            m = samples.mean()
            s = samples.std()
            for l in xrange ( nsegments ):
                z[l,k] = (samples[l*w:(l+1)*w].mean()-m)/s

        # warn about bad points
        if abs(z).max() > 2:
            bad = []
            for k in xrange(self.nchains):
                if abs(z[:,k]).max() > 2:
                    bad.append(k)
            return False,z
        else:
            return True,z

    def Rhat ( self, parameter=0 ):
        """Gelman Rhat statistic for convergence using multiple chains

        This is also called the 'estimated potential scale reduction'.
        A value Rhat > 1.1 is indicative of poor convergence.
        """
        # See p.137 in Gilks, Richardson, Spiegelhalter (1996)
        m = self.nchains
        n = self.getsamples(0).shape[0]

        psi_i = N.zeros(m,'d')    # within chain averages
        si2 = N.zeros(m,'d')      # within chain variances

        for chain in xrange(m):
            psi_i[chain] = self.getsamples(chain)[:,parameter].mean()
            si2[chain]   = sum ( (self.getsamples(chain)[:,parameter]-psi_i[chain])**2 )/(n-1)

        psi_mean = psi_i.mean()
        B = n * sum( (psi_i-psi_mean)**2 ) / (m-1)
        W = si2.mean()

        return (float(n-1)/n * W + B/n)/W;

    ############################################
    # Properties
    nchains = property ( fget=lambda self: len(self.__mcmc_chains), doc="Number of chains that have been sampled" )
    @Property
    def estimate ():
        """Estimate of the parameters.

        If sampling has already occurred, this will be the mean estimate, otherwise it will be the mapestimate.
        """
        def fget (self):
            if self.__meanestimate is None:
                # We don't have a mean estimate
                if len(self.__mcmc_chains) > 0:
                    # But we have samples!
                    self.__meanestimate = self.getsamples().mean(0)
                    self.devianceresiduals,self.__meandeviance,self.thres,self.Rpd,self.Rkd = _psipy.diagnostics ( \
                            self.data, self.__meanestimate, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[1:]
                else:
                    # We have no samples ~> return mapestimate
                    return self.mapestimate
            # In this case, we seem to have a meanestimate, so we return it
            return self.__meanestimate
        def fset (self,v):
            pass

    @Property
    def deviance ():
        """Deviance of the estimate.

        If sampling has already occurred, this will be the deviance of the mean estimate. Otherwise it will be
        the deviance of the mapestimate.
        """
        def fget (self):
            if self.__meandeviance is None:
                return self.mapdeviance
            else:
                return self.__meandeviance
        def fset (self,v):
            pass

    @Property
    def burnin ():
        "Burnin: Number of samples to be discarded at the beginning of each chain"
        def fget (self):
            return self.__burnin
        def fset (self,b):
            """Set the burnin

            :Parameters:
                b   new burnin value, i.e. number of samples that are discarded at
                    the beginning of each chain
            """
            self.__burnin = b
            # Set all values that depend on burnin to None. This way, they are
            # recomputed on access
            self.__meanestimate = None
            self.__meandeviance = None
            self.__pRpd = None
            self.__pRkd = None
            self.__pthres = None

    @Property
    def thin ():
        "Thinning: Subsample chains to reduce autocorrelation"
        def fget (self):
            return self.__thin
        def fset (self,t):
            self.__thin = t
            # Set all values that depend on thin to None. This way, they are recomputed
            # on access
            self.__meanestimate = None
            self.__meandeviance = None
            self.__pRpd = None
            self.__pRkd = None
            self.__pthres = None

    @Property
    def mcRpd ():
        "Monte Carlo samples of posterior correlation between model predictions and data"
        def fget (self):
            """Get samples from the posterior distribution of correlation between model prediction and deviance residuals"""
            if self.__pRpd is None:
                # pRpd is currently undefined
                if len(self.__mcmc_chains) > 0:
                    # We have samples ~> recompute the correlations
                    self.__recomputeCorrelationsAndThresholds()
                else:
                    raise NosamplesError, "Samples from the posterior have not yet been drawn"
            return self.__pRpd
        def fset (self, v):
            pass

    @Property
    def mcRkd ():
        "Monte Carlo samples of posterior correlation between bock index and data"
        def fget (self):
            """Get samples from the posterior distribution of correlation between block index and deviance residuals"""
            if self.__pRkd is None:
                # pRkd is currently undefined
                if len(self.__mcmc_chains) > 0:
                    # We have samples ~> recompute the correlations
                    self.__recomputeCorrelationsAndThresholds()
                else:
                    raise NosamplesError, "Samples from the posterior have not yet been drawn"
            return self.__pRkd
        def fset (self, v):
            pass

    mcestimates = property ( fget=getsamples, doc="Monte Carlo samples from the posterior distribution of parameters" )
    mcdeviance = property ( fget=getmcdeviance , doc="Deviances of monte carlo samples from the posterior" )

    @Property
    def mcthres ():
        "Monte Carlo Samples from the posterior distribution of thresholds"
        def fget (self):
            """Get samples of the posterior distribution of thresholds"""
            if self.__pthres is None:
                # pthres is currently undefined
                if len(self.__mcmc_chains) > 0:
                    # We have samples ~> recompute the thresholds
                    self.__recomputeCorrelationsAndThresholds()
                else:
                    raise NosamplesError, "Samples from the posterior have not yet been drawn"
            return self.__pthres
        def fset (self, t):
            pass

    @Property
    def evidence ():
        """model evidence or marginal likelihood

        Model evidence is typically given as the integral of the likelihood over the parameter space.
        We replace the integral by a discrete sum over samples, such that we have

        E = 1/N sum P(D|theta)

        Model evidence is typically used in Bayesian model selection: If E1 is the evidence for model
        1 and E2 is the evidence for model 2, then the ratio E1/E2 can be interpreted as "how much more
        evidence is there for model 2 than for model 1".
        """
        def fget (self):
            dev = self.mcdeviance
            return N.exp(-0.5*dev).mean()

    @Property
    def pD ():
        """effective number of parameters"""
        def fget ( self ):
            return self.mcdeviance.mean()-self.deviance

    @Property
    def DIC ():
        """Deviance information criterion

        This is an information criterion based on the posterior distribution of deviance.
        In contrast, to other information criteria, the deviance information criterion
        determines the effective number of free parameters from the posterior distribution.
        """
        def fget ( self ):
            meandev = self.mcdeviance.mean()
            return 2*meandev-self.deviance

    ############################################
    # Private methods
    def __recomputeCorrelationsAndThresholds ( self ):
        """This method is called whenever the sample basis from the
        posterior changes. This can have three reasons:
            - burnin: the burnin is changed resulting in samples being
                added or removed at the beginning of each chain
            - thin: the thinning is changed resulting in samples being
                discarded from within the chains
            - sample: an additional chain is acquired. In this case,
                a large number of samples is added.
        """
        samples = self.getsamples()

        self.__pRpd = N.zeros(samples.shape[0],'d')
        self.__pRkd = N.zeros(samples.shape[0],'d')
        self.__pthres = N.zeros((samples.shape[0],self.Ncuts),'d')

        for k,theta in enumerate(samples):
            self.__pthres[k,:],self.__pRpd[k],self.__pRkd[k] = _psipy.diagnostics (\
                    self.data, theta, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3:]

    def __determineoptimalsampling ( self, noptimizations=10, verbose=False ):
        """Determine optimal sampling parameters using the Raftery&Lewis (1995) procedure

        Automatically set burnin,thin,nsamples.
        In addition, an object, that contains more detailed information about the sampling
        is stored in self.mcmcpars

        :Parameters:
            *noptimizations* :
                maximum number of optimization iterations. If the same
                sampling parameters are obtained before, the method
                terminates earlier
            *verbose* :
                display status messages
        """
        if noptimizations==0:
            return
        mcmcpars = {}

        if len(self.__mcmc_chains)>0:
            mcmc_chains    = self.__mcmc_chains.copy()
            mcmc_deviances = self.__mcmc_deviances.copy()
            self.__mcmc_chains    = []
            self.__mcmc_deviances = []
        else:
            mcmc_chains = []
            mcmc_deviances = []

        # Determine size of initial test run
        if self.nsamples is None:
            NN = 0
            for q in self.conf:
                Nmin = pygibbsit.gibbsit ( q=q )["Nmin"]
                NN = max(NN,Nmin)
            self.nsamples = NN

        oldburnin = 0
        oldthin   = 1
        oldnsamples = NN
        for n in xrange ( noptimizations ):
            self.sample ()           # Test run
            testrun = self.mcthres    # Thresholds from testrun
            samples = self.__mcmc_chains.pop()      # throw the samples away, don't use them for "real" inference
            deviances = self.__mcmc_deviances.pop()

            # Check all desired thresholds
            for q in self.conf:
                for k in xrange ( self.Ncuts ):
                    try:
                        mcmcpars = pygibbsit.gibbsit ( testrun[:,k], q=q )
                    except IndexError:
                        continue
                    self.burnin = max ( self.burnin, mcmcpars.burnin )
                    self.thin   = max ( self.thin,   mcmcpars.thin )
                    self.nsamples = max ( self.nsamples, mcmcpars.Nsamples )
            self._steps = N.sqrt(N.diag(N.cov ( samples[self.burnin::self.thin].T )))

            if verbose:
                print "Burnin:",self.burnin,"Thinning:",self.thin,"Nsamples:",self.nsamples
                print "Steps:",self._steps
            if oldburnin==self.burnin and oldthin==self.thin and oldnsamples==self.nsamples:
                break
            else:
                oldburnin,oldthin,oldnsamples = self.burnin,self.thin,self.nsamples
        self.mcmcpars = mcmcpars

        if len(mcmc_chains)>0:
            self.__mcmc_chains = mcmc_chains
            self.__mcmc_deviances = mcmc_deviances


