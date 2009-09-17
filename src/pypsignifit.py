#!/usr/bin/env python

import numpy as N
import pylab as p
from scipy import stats
import _psipy
import psigniplot as pp

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

    desc = property ( fget=lambda self: "sigmoid: %(sigmoid)s\ncore: %(core)s\nnAFC: %(nafc)d" % self.model,
            doc="A short description of the employed model")
    outl = property ( fget=lambda self: self.__outl, doc="A boolean array indicating whether or not a block was an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean array indicating whether or not a block was an influential observation" )

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

    def plothistogram ( self, simdata, observed, xname, ax=None ):
        """plot a histogram and compare observed data to it

        :Parameters:
            simdata     an array of monte-carlo samples of the parameter of interest
            observed    observed value of the parameter of interest
            xname       name of the paramter of interest
            ax          axes object defining the area where the plot should go

        :Output:
            returns a boolean value indicating whether or not the Null-Hypothesis that
                observed was drawn from the same distribution as simdata is true
        """
        if ax is None:
            ax = p.axes()

        if xname[0] == "R":
            ax.hist ( simdata, bins=N.arange(-1,1,.1) )
            p.setp(ax,xlim=(-1,1))
        else:
            ax.hist ( simdata, bins=20 )

        xtics = p.getp(ax,"xticks")
        ytics = p.getp(ax,"yticks")
        p25,p975 = p.prctile ( simdata, (2.5,97.5) )

        xr = xtics.max()-xtics.min()
        yy = [ytics.min(),ytics.max()+0.02*xr]
        ax.plot ( [observed]*2, yy, 'r', linewidth=2 )
        ax.plot ( [p25]*2, yy, 'r:', [p975]*2, yy, 'r:' )

        pp.drawaxes ( ax, xtics, "%g", ytics, "%d", xname, "number per bin" )

        # Write diagnostics
        yt = ytics.max()
        ax.text ( xtics.min(), yt+.1, "%s=%.3f, c(2.5%%)=%.3f, c(97.5%%)=%.3f" % (xname,observed,p25,p975), horizontalalignment="left",verticalalignment="bottom", fontsize=8 )

        if observed>p25 and observed<p975:
            return True
        else:
            return False

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

    def diagnostics ( self ):
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
        self.plothistogram ( self.bdeviance, self.deviance, "deviance", p.axes([0,0,.33,.5]) )
        self.plotRd ( p.axes([.33,.5,.33,.5]), "p" )
        self.plothistogram ( self.bRpd, self.Rpd, "Rpd", p.axes([.33,0,.33,.5]) )
        self.plotRd ( p.axes([.66,.5,.33,.5]), "k" )
        self.plothistogram ( self.bRkd, self.Rkd, "Rkd", p.axes([.66,0,.33,.5]) )

    outl = property ( fget=lambda self: self.__outl, doc="A boolean vector indicating whether a block should be considered an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean vector indicating whether a block should be considered an influential observation" )
    bdeviance = property ( fget=lambda self: self.__bdeviance, doc="A vector of bootstrapped deviances" )
    bRpd = property ( fget=lambda self: self.__bRpd, doc="A vector of correlations between model prections and deviance residuals in all bootstrap samples" )
    bRkd = property ( fget=lambda self: self.__bRkd, doc="A vector of correlations between block index and deviance residuals in all bootstrap samples" )

class BayesInference ( PsiInference ):
    def __init__ ( self, data, sample=True, cuts=(.25,.5,.75), conf=(.025,.975), **kwargs ):
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

        self.burnin = 0
        self.thin   = 1

        if sample:
            self.sample ()

    def sample ( self, Nsamples=10000, start=None ):
        """Draw samples from the posterior distribution using MCMC"""
        if start is None:
            start = self.mapestimate
        # TODO: stepwidths are not good so far
        stepwidths = (0.4,4,1e-2)
        chain,deviance = _psipy.mcmc ( self.data, start, Nsamples, stepwidths=stepwidths, **self.model )
        # print N.cov(N.array(chain[self.burnin::self.thin]).T)
        self.__mcmc_chains.append(N.array(chain))
        self.__mcmc_deviances.append(N.array(deviance))

    def getsamples ( self, chain=None ):
        """Get the chain at index"""
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_chains:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        else:
            return self.__mcmc_chains[chain][self.burnin::self.thin]

    def getestimate ( self ):
        """Get the mean estimate"""
        if self.__meanestimate is None:
            if len(self.__mcmc_chains) > 0:
                self.__meanestimate = self.getsamples().mean(0)
                self.devianceresiduals,self.__meandeviance,self.thres,self.Rpd,self.Rkd = _psipy.diagnostics ( \
                        self.data, self.__meanestimate, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[1:]
            else:
                return self.mapestimate
        return self.__meanestimate

    def getdeviance ( self ):
        """Get the deviance"""
        if self.__meandeviance is None:
            return self.mapdeviance
        else:
            return self.__meandeviance

    def __recomputeCorrelationsAndThresholds ( self ):
        samples = self.getsamples()
        self.__pRpd = N.zeros(samples.shape[0],'d')
        self.__pRkd = N.zeros(samples.shape[0],'d')
        self.__pthres = N.zeros((samples.shape[0],self.Ncuts),'d')
        for k,theta in enumerate(samples):
            self.__pthres[k,:],self.__pRpd[k],self.__pRkd[k] = _psipy.diagnostics (\
                    self.data, theta, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3:]

    def getpRpd ( self ):
        """Get posterior distribution of correlation between model prediction and deviance residuals"""
        if self.__pRpd is None:
            if len(self.__mcmc_chains) > 0:
                self.__recomputeCorrelationsAndThresholds()
            else:
                raise Exception
        return self.__pRpd

    def getpRkd ( self ):
        """Get posterior distribution of correlation between block index and deviance residuals"""
        if self.__pRkd is None:
            if len(self.__mcmc_chains) > 0:
                self.__recomputeCorrelationsAndThresholds()
            else:
                raise Exception
        return self.__pRkd

    def getpthres ( self ):
        """Get posterior distribution of thresholds"""
        if self.__pthres is None:
            if len(self.__mcmc_chains) > 0:
                self.__recomputeCorrelationsAndThresholds()
            else:
                raise Exception
        return self.__pthres

    def getpdeviance ( self, chain=None ):
        """Get the chain at index"""
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_deviances:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        else:
            return self.__mcmc_deviances[chain][self.burnin::self.thin]

    def set_burnin ( self, b ):
        """Set the burnin"""
        self.__burnin = b
        self.__meanestimate = None
        self.__meandeviance = None
        self.__pRpd = None
        self.__pRkd = None

    def set_thin ( self, t ):
        self.__thin = t
        self.__meanestimate = None
        self.__meandeviance = None

    def getPI ( self, param="thres", conf=(.025,0.5,.975) ):
        if param=="thres":
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
                raise Exception
            return p.prctile ( vals, 100*N.array(conf) )

    def drawposteriorexamples ( self, ax=None, Nsamples=20 ):
        if ax is None:
            ax = p.axes()

        # Plot the psychometric function
        xmin = self.data[:,0].min()
        xmax = self.data[:,0].max()
        x = N.mgrid[xmin:xmax:100j]

        for k in xrange(Nsamples):
            chain = N.random.randint ( 0, len(self.__mcmc_chains) )
            sample = N.random.randint ( 0, self.__mcmc_chains[chain].shape[0] )
            psi = N.array(_psipy.diagnostics ( x, self.__mcmc_chains[chain][sample,:], sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
            ax.plot(x,psi,color=[.6,.6,1])

        for k,cut in enumerate(self.cuts):
            c25,c975 = self.getPI ( conf=(.025,.975))[k]
            thres = float(_psipy.diagnostics ( self.data, self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3])
            ylev  = _psipy.diagnostics ( [thres],   self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )
            ax.plot ( [c25,thres,c975],[ylev]*3, 'b-|' )

        self.pmfanddata ( ax=ax )

    def drawposteriorRd ( self, ax=None, regressor="p", Nsamples=50 ):
        if ax is None:
            ax = p.axes()

        # Make shure we have everything
        dummy = self.estimate

        for k in xrange(Nsamples):
            chain = N.random.randint ( 0, len(self.__mcmc_chains) )
            sample = N.random.randint ( 0, self.__mcmc_chains[chain].shape[0] )

            psi,devianceresiduals = _psipy.diagnostics ( self.data, self.__mcmc_chains[chain][sample,:], sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] )[:2]
            psi = N.array(psi)
            devianceresiduals = N.array(devianceresiduals)
            if regressor=="p":
                pass
            elif regressor=="k":
                psi = N.arange(len(self.data[:,0]))
            else:
                raise ValueError,"regressor %s is unknown" % regressor
            psilims = N.array([psi.min(),psi.max()])
            ax.plot ( psi, devianceresiduals, ".", color=(.6,.6,1) )

            # Linear regression
            A = N.ones((len(psi),2),'d')
            A[:,1] = psi
            a,b = N.linalg.lstsq(A,devianceresiduals)[0]
            ax.plot(psilims,a+b*psilims,':',color=(.6,.6,1))

        self.plotRd ( ax=ax, regressor=regressor )

    def plotposterior ( self, ax=None, param="threshold", cut=0 ):
        if ax is None:
            ax = p.axes()

        if param=='Rkd':
            val = self.pRkd
            ax.hist(val,N.arange(-1,1,.1))
        elif param=="Rpd":
            val = self.pRpd
            ax.hist(val,N.arange(-1,1,.1))
        else:
            if param=="threshold":
                val = self.pthres[:,cut]
            elif param=="deviance":
                val = self.pdeviance
            ax.hist(val,20)

        if param[0]=="R":
            p.setp(ax,xlim=(-1,1))

        xtics = p.getp(ax,"xticks")
        ytics = p.getp(ax,"yticks")

        c25,c975 = self.getPI ( param=param, conf = (.025,.975) )
        ydat = [ytics.min(),ytics.max()]

        p.plot ( [c25]*2, ydat, 'r:',[c975]*2, ydat,'r:' )
        p.plot ( [0]*2, ydat, 'r-', linewidth=2 )

        pp.drawaxes (ax,xtics,"%g", ytics, "%d", param, "number per bin" )

    def gof ( self ):
        p.figure(figsize=(10,8))
        self.drawposteriorexamples ( p.axes([0,.5,.33,.5] ) )
        self.plotposterior ( p.axes([0,0,.33,.5]), "deviance" )
        self.drawposteriorRd ( p.axes([.33,.5,.33,.5]), "p", Nsamples=0 )
        self.plotposterior ( p.axes([.33,0,.33,.5]), "Rpd" )
        self.drawposteriorRd ( p.axes([.66,.5,.33,.5]), "k", Nsamples=0 )
        self.plotposterior ( p.axes([.66,0,.33,.5]), "Rkd" )


    nchains = property ( fget=lambda self: len(self.__mcmc_chains), doc="Number of chains that have been sampled" )
    estimate = property ( fget=getestimate, fset=lambda self,v: v, doc="Determine the estimate of the parameters" )
    deviance = property ( fget=getdeviance, fset=lambda self,v: v, doc="Deviance of the estimate" )
    burnin  = property ( fget=lambda self: self.__burnin, fset=set_burnin, doc="Burnin: Number of samples to be discarded at the beginning of each chain" )
    thin    = property ( fget=lambda self: self.__thin,   fset=set_thin,   doc="Thinning: Subsample chains to reduce autocorrelation" )
    pRpd    = property ( fget=getpRpd, doc="Determine posterior correlation of model predictions and data" )
    pRkd    = property ( fget=getpRkd, doc="Determine posterior correlation of model predictions and data" )
    pdeviance = property ( fget=getpdeviance , doc="Deviances of the posterior samples" )
    pthres  = property ( fget=getpthres, doc="posterior distribution of thresholds" )

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
        b.diagnostics()
    else:
        priors = ("Gauss(0,100)","Gamma(1,3)","Beta(2,100)")
        mcmc = BayesInference ( d, priors=priors )
        mcmc.thin = 10
        mcmc.burnin = 200
        mcmc.sample(start=(6,4,.3))
        mcmc.sample(start=(1,1,.1))
        mcmc.gof()
        print mcmc.getPI()

    p.show()

if __name__ == "__main__":
    main()
