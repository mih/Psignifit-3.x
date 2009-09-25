#!/usr/bin/env python

import sys
import numpy as N
import pylab as p
from scipy import stats
import _psipy

import psigniplot as pp
import pygibbsit

from psignierrors import NosamplesError

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
        self.__outl              = None
        self.__infl              = None

    def pmfanddata ( self, ax=None, xlabel_text="Stimulus intensity", ylabel_text=None ):
        """Show the psychometric function and data in an axes system

        This function plots the best fitting psychometric function and with the
        corresponding data points. If data points are labelled influential, they
        are plotted as red squares, if data points are labelled as outliers, they
        are plotted as red triangles.
        The function uses its internal knowledge about the task (nAFC or Yes/No)
        to put the correct labels to the y-axis.

        :Parameters:
            ax          axes object in which the plot should go
            xlabel_text label for the x-axis
            ylabel_text label for the y-axis, if this is None, the functions
                        determines the correct label from its internal knowledge
                        about the task
        """
        if ax==None:
            ax = p.axes()

        # Plot the psychometric function
        xmin = self.data[:,0].min()
        xmax = self.data[:,0].max()
        x = N.mgrid[xmin:xmax:100j]
        psi = N.array(_psipy.diagnostics ( x, self.estimate, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
        ax.plot(x,psi,'b')

        # Plot the data
        xd = self.data[:,0]
        pd = self.data[:,1].astype("d")/self.data[:,2]
        goodpoints = N.ones(pd.shape,bool)
        if not self.outl==None:
            ax.plot(xd[self.outl],pd[self.outl],'r^')
            goodpoints = N.logical_and(goodpoints,N.logical_not(self.outl))
        if not self.infl==None:
            ax.plot(xd[self.infl],pd[self.infl],"rs")
            goodpoints = N.logical_and(goodpoints,N.logical_not(self.infl))
        ax.plot(xd[goodpoints],pd[goodpoints],'bo')

        # Check axes limits
        if self.model["nafc"]>1:
            ymin,ymax = 1./self.model["nafc"]-.05,1.05
            if ylabel_text is None:
                ylabel_text = "P(correct)"
        else:
            ymin,ymax = -.05,1.05
            if ylabel_text is None:
                ylabel_text = "P(Yes)"

        # Determine tics
        p.setp(ax,frame_on=False,ylim=(ymin,ymax))
        xtics = p.getp(ax,'xticks')
        ytics = p.getp(ax,'yticks').tolist()
        # Clean up ytics
        if self.model["nafc"]==1:
            for k,yt in enumerate(ytics):
                if yt<0 or yt>1:
                    ytics.pop(k)
        else:
            for k,yt in enumerate(ytics):
                if yt<(1./self.model["nafc"]) or yt>1:
                    ytics.pop(k)
        ytics = N.array(ytics)

        pp.drawaxes ( ax, xtics, "%g", ytics, "%g", xlabel_text, ylabel_text )

        # Write some model information
        if not self.deviance is None:
            ax.text(0.5*(xmin+xmax),ymin+.05,"D=%g" % ( self.deviance, ) )
        ax.text ( 0.5*(xmin+xmax),ymin+.1,self.desc )

    def plotRd ( self, ax=None, regressor="p" ):
        """plot deviance residuals against a regressor

        Deviance residuals are used plotted agains either predicted performance or
        block index to check for systematic deviations of the data from the fitted
        function.

        :Parameters:
            ax          an axes object where the plot should go
            regressor   plot deviance residuals against model prediction (p) or
                        against block index (k)
        """
        if ax==None:
            ax = p.axes()

        # Plot the data points
        if regressor=="p":
            psi = N.array(_psipy.diagnostics ( self.data[:,0], self.estimate, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
        elif regressor=="k":
            psi = N.arange(len(self.data[:,0]))
        else:
            raise ValueError,"regressor %s is unknown" % regressor
        psilims = N.array([psi.min(),psi.max()])
        ax.plot ( psi, self.devianceresiduals, "bo" )

        # Linear regression
        A = N.ones((len(psi),2),'d')
        A[:,1] = psi
        a,b = N.linalg.lstsq(A,self.devianceresiduals)[0]
        ax.plot(psilims,a+b*psilims,'b:')

        if regressor=="p":
            if self.model["nafc"]==1:
                p.setp(ax,xlim=(0,1))
            else:
                p.setp(ax,xlim=(1./self.model["nafc"],1))
        xtics = p.getp(ax,"xticks").tolist()
        if regressor=="p":
            # In this case predictions larger than 1 and less than 0 are impossible
            for k,xt in enumerate(xtics):
                if xtics[k]>1. or xtics[k]<0.:
                    xtics.pop(k)
        xtics = N.array(xtics)
        ytics = p.getp(ax,"yticks")

        # Generate the respective labels
        if regressor=="p":
            ax.text(psilims.mean(),ytics[-2],"Rpd=%.3f" % ( self.Rpd, ) )
            xname = "model prediction"
        elif regressor=="k":
            ax.text(psilims.mean(),ytics[-2],"Rkd=%.3f" % ( self.Rkd, ) )
            xname = "block index"

        pp.drawaxes ( ax, xtics, "%g", ytics, "%g", xname, "deviance residuals" )

    def plothistogram ( self, simdata, observed, xname, shortname=None, ax=None, hideobserved=False ):
        """plot a histogram and compare observed data to it

        :Parameters:
            simdata     an array of monte-carlo samples of the parameter of interest
            observed    observed value of the parameter of interest (for MCMC samples, it is often
                        reasonable to use this as the value of 'no effect' or something)
            xname       name of the paramter of interest
            shortname   short name of the parameter of interest
            ax          axes object defining the area where the plot should go
            hideobserved if this is True, the observed value is not plotted

        :Output:
            returns a boolean value indicating whether or not the Null-Hypothesis that
                observed was drawn from the same distribution as simdata is true
        """
        if ax is None:
            ax = p.axes()

        # Make sure we have a useful shortname
        if shortname is None:
            shortname = xname

        # Correlations plots should be treated differently
        if shortname[0] == "R":
            ax.hist ( simdata, bins=N.arange(-1,1,.1) )
            p.setp(ax,xlim=(-1,1))
        else:
            ax.hist ( simdata, bins=20 )

        # Get the tics and ranges
        xtics = p.getp(ax,"xticks")
        ytics = p.getp(ax,"yticks")
        xr = xtics.max()-xtics.min()
        yy = [ytics.min(),ytics.max()+0.02*xr]

        # Plot percentile bars
        p25,p975 = p.prctile ( simdata, (2.5,97.5) )
        if not hideobserved:
            ax.plot ( [observed]*2, yy, 'r', linewidth=2 )
        ax.plot ( [p25]*2, yy, 'r:', [p975]*2, yy, 'r:' )

        # Draw the full plot
        pp.drawaxes ( ax, xtics, "%g", ytics, "%d", xname, "number per bin" )

        # Write diagnostics
        yt = ytics.max()
        ax.text ( xtics.min(), yt+.1, "%s=%.3f, c(2.5%%)=%.3f, c(97.5%%)=%.3f" % (shortname,observed,p25,p975),\
                horizontalalignment="left",verticalalignment="bottom", fontsize=8 )

        if observed>p25 and observed<p975:
            return True
        else:
            return False

    desc = property ( fget=lambda self: "sigmoid: %(sigmoid)s\ncore: %(core)s\nnAFC: %(nafc)d" % self.model,
            doc="A short description of the employed model")
    outl = property ( fget=lambda self: self.__outl, doc="A boolean array indicating whether or not a block was an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean array indicating whether or not a block was an influential observation" )

##############################################################################################################################
class BootstrapInference ( PsiInference ):
    def __init__ ( self, data, sample=False, cuts=(.25,.5,.75), conf=(.025,.975), **kwargs ):
        """Set up an object of bootstrapped data

        :Parameters:
            data    an array or a list of lists containing stimulus intensities in the
                    first column, number of correct responses (nAFC) or number of YES-
                    responses in the second column, and number of trials in the third
                    column. Each row should correspond to one experimental block. In
                    addition, the sequence of the rows is taken as the sequence of
                    data aquisition.
            sample  if sample is True, bootstrap samples are drawn. If sample is an
                    integer, it gives the number of samples that are drawn
            sigmoid shape of the sigmoid function. Valid choices are
                        'logistic'   [Default]
                        'gauss'
                        'gumbel_l'
                        'gumbel_r'
            core    term inside the sigmoid function. Valid choices are
                        'ab'         (x-a)/b        [Default]
                        'mw%g'       midpoint and width
                        'linear'     a+b*x
                        'log'        a+b*log(x)
            priors  a list of prior names. Valid choices are
                        'Uniform(%g,%g)'   Uniform distribution on an interval
                        'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                        'Beta(%g,%g)'      Beta distribution
                        'Gamma(%g,%g)'     Gamma distribution
                    If no valid prior is selected, the parameter remains unconstrained.
                    Alternatively, priors can be given as a dictionary that only specifies
                    priors for those parameters you want to set in that case you can use
                    'a','b','m','w','guess','gamma','lapse','lambda' as keys.
            nafc    number of response alternatives. If nafc==1, this indicates a Yes/No
                    task
            cuts    performance values that should be considered 'thresholds'. This means that a
                    'cut' of 0.5 corresponds to an expected performance of roughly 75%% correct in
                    a 2AFC task.
            conf    limits of confidence intervals. The default gives 95%% confidence intervals.
                    Any other sequence can be used alternatively. In addition, conf can be 'v1.0'
                    to give the default values of the classical psignifit version (i.e. .023,.159,.841,.977,
                    corresponding to -2,-1,1,2 standard deviations for a gaussian).
        """
        # Call the base constructor
        PsiInference.__init__(self)

        # Store basic data
        self.data = data
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "ab"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }

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

    def sample ( self, Nsamples=2000 ):
        """Draw bootstrap samples

        :Parameters:
            Nsamples    number of bootstrapsamples to be drawn
        """
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
        self.__outl  = N.array(self.__outl,dtype=bool)
        self.__infl  = N.array(self.__infl,dtype=bool)

    def drawthreshold ( self, ax, cut ):
        """draw the threshold into an axes system

        :Parameters:
            ax      axes system in which to draw
            cut     index(!) of the cut of interest
        """
        cc = self.getCI(cut)

        ylev =_psipy.diagnostics ( [self.thres[cut]], self.estimate, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] )

        ax.plot ( cc, [ylev]*len(cc), 'b-|' )
        ax.plot ( [self.thres[cut]], [ylev], 'b|' )

    def getCI ( self, cut ):
        """Determine the confidence interval of a cut

        :Parameters:
            cut     index(!) of the cut of interest
        """
        bias = self.__th_bias[cut]
        acc  = self.__th_bias[cut]

        vals = []
        for pp in self.conf:
            vals.append(stats.norm.cdf( bias + ( stats.norm.ppf(pp) + bias ) / (1-acc*(stats.norm.ppf(pp) + bias )) ))

        return p.prctile ( self.__bthres[:,cut], 100*N.array(vals) )

    def gof ( self, warn=True ):
        """Make a diagnostic plot

        This plots a figure that resembles the diagnostic plot of the old psignifit matlab interface.
        The figure has the following form:

            +-----+-----+-----+
            |  1  |  2  |  3  |
            +-----+-----+-----+
            |  4  |  5  |  6  |
            +-----+-----+-----+

        At position 1, the psychometric function is shown with the fitted data. Influential observations
            are marked as red squares, outliers as red triangles.
        At position 2, the deviance residuals are plotted against model predictions. The best fitting
            straigt line is shown in addition. This plot should not show any obvious trends.
        At position 3, the deviance residuals are plotted against block index. The best fitting
            straigt line is shown in addition. This plot should not show any obvious trends. Trends
            in this plot might indicate perceptual learning (which is interesting in itself but
            makes the statistics ambiguous).
        At position 4, a histogram of deviances for an observer that perfectly agrees with the model fit
            in 1 is shown. In addition, the observed deviance is shown as a solid red line and the
            2.5% and 97.5% percentile are drawn as dotted red lines. If the observed deviance is outside
            the interval marked by the two dotted red lines, this indicates a bad fit.
        At position 5 and 6, histograms of the correlations between deviance residuals and model
            prediction (5) or block index (6) are shown for the assumption of an observer that perfectly
            agrees with the model fit in 1. The observed correlation (from positions 2 or 3) is drawn
            as a solid red line and the 2.5% and 97.5% percentile are marked by dotted red lines.
            If the observed correlation is outside the interval marked by the two dotted red lines, this
            indicates that the data systematically deviate from the model. In this case, the model
            does not capture all the structure in the data.
        """
        p.figure(figsize=(10,8))
        self.pmfanddata ( p.axes([0,.5,.33,.5] ) )
        for k in xrange(len(self.cuts)):
            self.drawthreshold ( p.gca(), k )
        self.plothistogram ( self.bdeviance, self.deviance, "deviance", "D", p.axes([0,0,.33,.5]) )
        self.plotRd ( p.axes([.33,.5,.33,.5]), "p" )
        self.plothistogram ( self.bRpd, self.Rpd, "Rpd", "Rpd", p.axes([.33,0,.33,.5]) )
        self.plotRd ( p.axes([.66,.5,.33,.5]), "k" )
        self.plothistogram ( self.bRkd, self.Rkd, "Rkd", "Rkd", p.axes([.66,0,.33,.5]) )

    outl = property ( fget=lambda self: self.__outl, doc="A boolean vector indicating whether a block should be considered an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean vector indicating whether a block should be considered an influential observation" )
    bdeviance = property ( fget=lambda self: self.__bdeviance, doc="A vector of bootstrapped deviances" )
    bRpd = property ( fget=lambda self: self.__bRpd, doc="A vector of correlations between model prections and deviance residuals in all bootstrap samples" )
    bRkd = property ( fget=lambda self: self.__bRkd, doc="A vector of correlations between block index and deviance residuals in all bootstrap samples" )

##############################################################################################################################
class BayesInference ( PsiInference ):
    def __init__ ( self, data, sample=True, cuts=(.25,.5,.75), conf=(.025,.975), automatic=True, resample=False, **kwargs ):
        """Bayesian Inference for psychometric functions using MCMC

        :Parameters:
            data    an array or a list of lists containing stimulus intensities in the
                    first column, number of correct responses (nAFC) or number of YES-
                    responses in the second column, and number of trials in the third
                    column. Each row should correspond to one experimental block. In
                    addition, the sequence of the rows is taken as the sequence of
                    data aquisition.
            sample  if sample is True, bootstrap samples are drawn. If sample is an
                    integer, it gives the number of samples that are drawn
            sigmoid shape of the sigmoid function. Valid choices are
                        'logistic'   [Default]
                        'gauss'
                        'gumbel_l'
                        'gumbel_r'
            core    term inside the sigmoid function. Valid choices are
                        'ab'         (x-a)/b        [Default]
                        'mw%g'       midpoint and width
                        'linear'     a+b*x
                        'log'        a+b*log(x)
            priors  a list of prior names. Valid choices are
                        'Uniform(%g,%g)'   Uniform distribution on an interval
                        'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                        'Beta(%g,%g)'      Beta distribution
                        'Gamma(%g,%g)'     Gamma distribution
                    If no valid prior is selected, the parameter remains unconstrained.
                    Alternatively, priors can be given as a dictionary that only specifies
                    priors for those parameters you want to set in that case you can use
                    'a','b','m','w','guess','gamma','lapse','lambda' as keys.
            nafc    number of response alternatives. If nafc==1, this indicates a Yes/No
                    task
            cuts    performance values that should be considered 'thresholds'. This means that a
                    'cut' of 0.5 corresponds to an expected performance of roughly 75%% correct in
                    a 2AFC task.
            conf    limits of confidence intervals. The default gives 95%% confidence intervals.
                    Any other sequence can be used alternatively. In addition, conf can be 'v1.0'
                    to give the default values of the classical psignifit version (i.e. .023,.159,.841,.977,
                    corresponding to -2,-1,1,2 standard deviations for a gaussian).
            automatic   do everything automatically
            resample if a chain is considered "bad" in terms of convergence should it
                    automatically be resampled?
        """
        PsiInference.__init__(self)

        # Store basic data
        self.data = data
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
            self.determineoptimalsampling ()
            sample = True

        if sample:
            self.sample()

    def determineoptimalsampling ( self, noptimizations=10 ):
        """Determine optimal sampling parameters using the Raftery&Lewis (1995) procedure

        Automatically set burnin,thin,nsamples.
        In addition, an object, that contains more detailed information about the sampling
        is stored in self.mcmcpars

        :Parameters:
            noptimizations  maximum number of optimization iterations. If the same
                            sampling parameters are obtained before, the method
                            terminates earlier
        """
        if noptimizations==0:
            return
        mcmcpars = {}

        if len(self.__mcmc_chains)>0:
            mcmc_chains    = self.__mcmc_chains
            mcmc_deviances = self.__mcmc_deviances
            self.__mcmc_chains    = []
            self.__mcmc_deviances = []
        else:
            mcmc_chains = []

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
            testrun = self.pthres    # Thresholds from testrun
            samples = self.__mcmc_chains.pop() # throw the samples away, don't use them for "real" inference

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
            print "Steps:",self._steps

            print "Burnin:",self.burnin,"Thinning:",self.thin,"Nsamples:",self.nsamples
            if oldburnin==self.burnin and oldthin==self.thin and oldnsamples==self.nsamples:
                break
            else:
                oldburnin,oldthin,oldnsamples = self.burnin,self.thin,self.nsamples
        self.mcmcpars = mcmcpars

        if len(mcmc_chains)>0:
            self.__mcmc_chains = mcmc_chains
            self.__mcmc_deviances = mcmc_deviances

    def sample ( self, Nsamples=None, start=None ):
        """Draw samples from the posterior distribution using MCMC

        :Parameters:
            Nsamples    number of samples that should be drawn from the posterior. If Nsamples is
                        None, an optimal number of samples is tried to obtain.
            start       starting value of the chain. If this is None, the chain starts
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

    ############################################
    # Setters and getters
    def getsamples ( self, chain=None ):
        """Get sampes from the posterior

        :Parameters:
            chain   if chain is None, samples are aggregated over all chains
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
            return self.__mcmc_chains[chain][self.burnin::self.thin]
        else:
            raise IndexError, "chain should be either None or an integer"

    def getpdeviance ( self, chain=None ):
        """Get samples from the posterior distribution of deviances

        :Parameters:
            chain   if chain is None, the samples are combined across all chains
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

    def getPI ( self, conf=(.025,0.5,.975), param="thres" ):
        """Get a posterior interval for a particular parameter

        :Parameters:
            conf    percentiles that should be returned
            param   parameter of interest. Currently, only thres/threshold
                    and Rkd,Rpd,deviance are defined.
        """
        if param[:5]=="thres":
            # We have to handle thresholds separately because there could be multiple cuts.
            pthres = self.pthres
            out = []
            for k in xrange(self.Ncuts):
                out.append(p.prctile ( pthres[:,k], 100*N.array(conf) ))
            return N.array(out)
        else:
            if param=="Rkd":
                vals = self.pRkd
            elif param=="Rpd":
                vals = self.pRpd
            elif param=="deviance":
                vals = self.pdeviance
            else:
                raise NotImplementedError
            return p.prctile ( vals, 100*N.array(conf) )

    ############################################
    # Plotting routines
    def drawposteriorexamples ( self, ax=None, Nsamples=20 ):
        """plots the mean estimate of the psychometric function and a number of samples from the posterior

        :Parameters:
            ax      axes object in which to draw the plot. If this is None,
                    a new axes object is created.
            Nsamples number of psychometric functions that should be drawn
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
        for k in xrange(Nsamples):
            psi = N.array(_psipy.diagnostics ( x, samples[N.random.randint(samples.shape[0]),:], sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
            ax.plot(x,psi,color=[.6,.6,1])

        # We further visualize the posterior distribution py plotting some posterior intervals
        for k,cut in enumerate(self.cuts):
            c25,c975 = self.getPI ( conf=(.025,.975) )[k]
            thres = float(_psipy.diagnostics ( self.data, self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3])
            ylev  = _psipy.diagnostics ( [thres],   self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )
            ax.plot ( [c25,thres,c975],[ylev]*3, 'b-|' )

        # This plots the 'real' psychometric function and the axes
        self.pmfanddata ( ax=ax )

    def gof ( self, warn=True ):
        """Draw a diagnostic figure to help assessing goodness of fit

        This graphic is intended to help the user determine how well the fitted function describes
        the data. The plot has 6 fields:

        +-----+-----+-----+
        |  1  |  3  |  5  |
        +-----+-----+-----+
        |  2  |  4  |  6  |
        +-----+-----+-----+

        The fields provide the following information:
        1.  The data and the fitted psychometric function. "fitted" here means the parameters are
            the mean of the posterior. To get an idea of the posterior distribution, posterior
            intervals are plotted at some positions (the location and width of the posterior
            intervals is given in the constructor). To make the posterior distribution really
            "plastic", a number of samples from the posterior distribution over psychometric
            functions are also drawn in light blue
        2.  A histogram to approximate the posterior distribution of deviances.
        3.  A plot of model predictions (of the mean estimate) against deviance residuals. If
            there is no obvious interrelation between model prediction and deviance residuals,
            this indicates that the model describes the data reasonably well. To get an idea
            of the interrelation between model prediction and deviance residuals, the best
            fitting line is plotted as a dotted line.
        4.  A histogram of samples from the posterior distribution of correlations between
            model prediction and deviance residuals. If this distribution is clearly shifted
            away from 0, this is strong evidence, that something is wrong with your model or
            your data.
        5,6 Similar to 3 and 4 but form correlations between block index and deviance residuals.
            Correlations between block index and deviance residuals indicate nonstationary
            data as should be found during e.g. perceptual learning.

        :Parameters:
            warn    if warn is set to True, red warning messages are displayed
                    whenever the fit does not seem to describe the data well.
        """
        p.figure(figsize=(10,8))

        # First part: Data and fitted function, bottom deviance
        self.drawposteriorexamples ( p.axes([0,.5,.33,.5] ) )
        self.plothistogram ( self.pdeviance, self.deviance, "posterior deviance", "D", p.axes ( [0,0,.33,.5] ) )

        # Second part: Correlations between model prediction and residuals
        self.plotRd ( p.axes([.33,.5,.33,.5]), "p" )
        good = self.plothistogram ( self.pRpd, self.Rpd, "posterior Rpd", "Rpd", p.axes([.33,0,.33,.5]) )
        if not good and warn==True:
            p.text ( 0, p.getp(p.gca(),'ylim').mean() , "Rpd is different from 0!\nModel deviates systematically from data", \
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

        # Third part: Correlations between model prediction and block index
        self.plotRd ( p.axes([.66,.5,.33,.5]), "k" )
        good = self.plothistogram ( self.pRkd, self.Rkd, "posterior Rkd", "Rkd", p.axes([.66,0,.33,.5]) )
        if not good and warn==True:
            p.text ( 0, p.getp(p.gca(),'ylim').mean(), "Rkd is different from 0!\nData are nonstationary!",\
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    def convergence ( self, parameter=0, warn=True ):
        """A simple convergence plot"""
        p.figure()
        ax = p.axes([0,.5,.5,.5])
        self.chainplot ( parameter,ax=ax, warn=warn )
        ax = p.axes([.5,.5,.5,.5])
        self.geweke ( parameter, ax=ax, warn=warn )

    ############################################
    # Convergence diagnostics
    def geweke ( self, parameter=0, nsegments=10, ax=None, warn=True ):
        """Geweke test for stationarity of a chain.

        The Geweke test first transforms all samples to mean 0 and standard deviation 1.
        In a second step, it calculates the sample average in a number of segments and
        checks whether these subaverages differ significantly from 0.

        :Parameters:
            parameter       parameter of interest
            nsegments       number of subaverages to be calculated
            ax              axes to plot in (if this is None, no plot is created)
            warn            should a warning message be plotted if the chain did not
                            converge?

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
            if not ax is None:
                p.plot(z[:,k],'o-')
        if not ax is None:
            xtics = N.array(p.getp(ax,"xticks"))
            p.setp(ax,"xticks",xtics,"yticks",N.array((-3,-2,-1,0,1,2,3)))
            p.plot ( [xtics.min(),xtics.max()],[-2]*2,'k:')
            p.plot ( [xtics.min(),xtics.max()],[ 2]*2,'k:')
            pp.drawaxes ( ax, xtics, "%g", N.array((-3,-2,-1,0,1,2,3)), "%g", "chain segment", "z-score" )

        # warn about bad points
        if abs(z).max() > 2:
            bad = []
            for k in xrange(self.nchains):
                if abs(z[:,k]).max() > 2:
                    bad.append(k)
            if warn:
                p.text(0.5*nsegments,0,"chains did not converge: %s" % (bad,), color="red", fontsize=16, rotation=45, verticalalignment="center", horizontalalignment="center" )
            return False
        else:
            return True

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

    def chainplot ( self, parameter=0, raw=False, ax=None, warn=True ):
        """Simply plot all chains for a single parameter

        :Parameters:
            parameter   index of the model parameter to plot
            raw         plot raw samples instead of thinned samples after burnin
            ax          axes in which to print
            warn        if True, warnings are written into the plot
        """
        # Do we have an appropriate axis?
        if ax==None:
            ax = p.axes()

        # Plot the chains
        for c in xrange(self.nchains):
            if raw:
                samples = samples.__mcmc_chains[c]
            else:
                samples = self.getsamples(c)
            p.plot ( samples[:,parameter] )

        # Learn something about the axes
        xtics = N.array(ax.get_xticks())
        x0    = xtics.min()
        xr    = xtics.max()-xtics.min()
        ytics = ax.get_yticks()
        y0    = ytics.min()
        yr    = N.array(ytics.max()-ytics.min())

        p.text(x0+0.6*xr,y0+0.95*yr,"R^ = %.4f" % (self.Rhat ( parameter ) ) )

        if warn and self.Rhat(parameter)<1.1:
            p.text(x0+0.5*xr,y0+0.5*yr,"Chains do not seem to sample\nfrom the same distribution!",
                    horizontalalignment="center",verticalalignment="center",fontsize=16,rotation=45,color=warnred)

        pp.drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", "sample #", self.parnames[parameter] )

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
    def pRpd ():
        "Determine posterior correlation of model predictions and data"
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
    def pRkd ():
        "Determine posterior correlation of model predictions and data"
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

    pdeviance = property ( fget=getpdeviance , doc="Deviances of the posterior samples" )

    @Property
    def pthres ():
        "Samples from the posterior distribution of thresholds"
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
            dev = self.pdeviance
            return N.exp(-0.5*dev).mean()

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

def main ( ):
    "If we call the file directly, we perform a test run"

    bootstrap = False

    x = [float(2*k) for k in xrange(6)]
    k = [34,32,40,48,50,48]
    n = [50]*6
    d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    d = N.array(zip(x,k,n))
    priors = ("flat","flat","Uniform(0,0.1)")

    if bootstrap:
        b = BootstrapInference ( d, sample=2000, priors=priors )
        b.gof()
    else:
        priors = ("Gauss(0,100)","Gamma(1,3)","Beta(2,100)")
        mcmc = BayesInference ( d, priors=priors )
        # mcmc.thin = 10
        # mcmc.burnin = 200
        mcmc.sample(start=(6,4,.3))
        mcmc.sample(start=(1,1,.1))
        mcmc.gof()
        print "Posterior Intervals",mcmc.getPI()
        print "Model Evidence", mcmc.evidence
        print "Rhat (m):",mcmc.Rhat ()
        print "Nsamples:",mcmc.nsamples

        mcmc.convergence(0)

    p.show()

if __name__ == "__main__":
    main()
