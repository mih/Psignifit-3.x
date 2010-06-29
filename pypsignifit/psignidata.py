#!/usr/bin/env python
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

__docformat__ = "restructuredtext"

import sys,os
import operator
import numpy as N
import pylab as p
from scipy import stats,special,optimize
import _psipy as interface

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
    def __init__ ( self, plotting=None ):
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
        if plotting is None:
            self.__plotting        = {}
        else:
            self.__plotting = plotting
        defaults = {"label": "Psychometric function fit","color": "b", "linestyle": "-", "marker": "o", "linewidth": 1 }
        for k in defaults.keys():
            self.__plotting.setdefault ( k, defaults[k] )

    def evaluate ( self, x, prm=None ):
        """Evaluate the psychometric function model at positions given by x"""
        if prm==None:
            prm = self.estimate

        return N.array( interface.diagnostics ( x, prm, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ) )

    def getThres ( self, cut=0.5 ):
        """Get thresholds at cut"""
        if self.data == None:
            raise NotImplementedError
        return float(interface.diagnostics ( self.data, self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3])

    def __repr__ ( self ):
        return "< PsiInference object >"

    desc = property ( fget=lambda self: "sigmoid: %(sigmoid)s\ncore: %(core)s\nnAFC: %(nafc)d" % self.model,
            doc="A short description of the employed model")
    outl = property ( fget=lambda self: self.__outl, doc="A boolean array indicating whether or not a block was an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean array indicating whether or not a block was an influential observation" )

    @Property
    def label ():
        "Condition label used in plots"
        def fget ( self ):
            return self.__plotting["label"]
        def fset ( self, v ):
            self.__plotting["label"] = v
    @Property
    def color ():
        "Color used in plots"
        def fget ( self ):
            return self.__plotting["color"]
        def fset ( self, v ):
            self.__plotting["color"] = v
    @Property
    def linestyle ():
        "Linestyle used in plots"
        def fget ( self ):
            return self.__plotting["linestyle"]
        def fset ( self, v ):
            self.__plotting["linestyle"] = v
    @Property
    def linewidth ():
        "Linewidth used in plots"
        def fget ( self ):
            return self.__plotting["linewidth"]
        def fset ( self ):
            self.__plotting["linewidth"] = v
    @Property
    def marker ():
        "Data marker used in plots"
        def fget ( self ):
            return self.__plotting["marker"]
        def fset ( self, v ):
            self.__plotting["marker"] = v

##############################################################################################################################
class BootstrapInference ( PsiInference ):
    def __init__ ( self, data, sample=False, cuts=(.25,.5,.75), conf=(.025,.975), plotprm=None, **kwargs ):
        """Set up an object of bootstrapped data

        :Parameters:
            *data* :
                an array or a list of lists containing stimulus intensities in the
                first column, number of correct responses (nAFC) or number of YES-
                responses in the second column, and number of trials in the third
                column. Each row should correspond to one experimental block. In
                addition, the sequence of the rows is taken as the sequence of
                data aquisition. Alternatively, the relative frequencies of correct
                responses resp YES responses can be given.
            *sample* :
                if sample is True, bootstrap samples are drawn. If sample is an
                integer, it gives the number of samples that are drawn
            *sigmoid* :
                shape of the sigmoid function. Valid choices are
                    - 'logistic'   [Default]
                    - 'gauss'
                    - 'gumbel_l'
                    - 'gumbel_r'
                    - 'exp'
            *core* :
                term inside the sigmoid function. Valid choices are
                    - 'ab'         (x-a)/b        [Default]
                    - 'mw%g'       midpoint and width
                    - 'linear'     a+b*x
                    - 'log'        a+b*log(x)
                    - 'weibull'    2*s*m*(log(x)-log(m))/log(2) + log(log(2))    This will give you a weibull if combined with
                      the gumbel_l sigmoid and a reverse weibull if combined with the gumbel_r sigmoid.
                    - 'poly'       (x/a)**b   Will give you a weibull if combined with an exp sigmoid
            *priors* :
                a list of prior names. Valid choices are
                    - 'Uniform(%g,%g)'   Uniform distribution on an interval
                    - 'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                    - 'Beta(%g,%g)'      Beta distribution
                    - 'Gamma(%g,%g)'     Gamma distribution
                    - 'nGamma(%g,%g)'    Gamma distribution on the negative axis
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
            *parametric* :
                a boolean indicating wether or not parametric bootstrap should be used
            *plotprm* :
                a dictionary to take parameters for plotting data. Currently supported are the arguments
                'label', 'color', 'linestyle', 'linewidth' and 'marker'. These can all be set after creating
                an Inference instance, too. By using the respective properties.


        :Example:
        Estimate a psychometric function from some example data and derive bootstrap confidence
        intervals for a threshold

        >>> x = [0,2,4,6,8,10]
        >>> k = [26,30,38,46,50,49]
        >>> n = [50]*len(k)
        >>> B = BootstrapInference ( zip(x,k,n), priors=("flat","flat","Uniform(0,0.1)"), sample=True )
        >>> B.estimate
        array([ 3.80593409,  1.09308994,  0.00935698])
        >>> B.deviance
        2.5160989036891754
        >>> B.getThres()
        3.805934094097025
        >>> B.getCI(1)
        array([ 2.79484448,  4.73796576])
        """
        # Call the base constructor
        PsiInference.__init__(self,plotprm)
        self.__nsamples = 0

        start = kwargs.setdefault ( "start", None )
        kwargs.pop("start")

        # Store basic data
        self.data = N.array(data)
        if self.data[:,1].max() <= 1:
            # We have relative frequencies
            self.data[:,1] *= self.data[:,2]
            self.data[:,1] = N.floor ( self.data[:,1] )
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "ab"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }

        self.parametric = kwargs.setdefault ( "parametric", True )

        if self.model["core"][:2] == "mw":
            self.parnames = ["m","w"]
        elif self.model["core"] == "weibull":
            self.parnames = ["m","s"]
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
        self.estimate,self.fisher,self.thres,self.deviance = interface.mapestimate(self.data,cuts=self.cuts,start=start,**self.model)
        self.predicted,self.devianceresiduals,self.deviance,thres,self.Rpd,self.Rkd = interface.diagnostics(self.data,self.estimate, \
                nafc=self.model["nafc"],sigmoid=self.model["sigmoid"],core=self.model["core"])

        # The interface arrays are not numpy arrays
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
        self.__expanded  = False

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
        # print self.estimate
        self.__bdata,self.__bestimate,self.__bdeviance,self.__bthres,self.__th_bias,self.__th_acc,\
                self.__bRkd,self.__bRpd,self.__outl,self.__infl = interface.bootstrap(self.data,self.estimate,Nsamples,
                        cuts=self.cuts,**self.model)
        if not self.parametric:
            self.sample_nonparametric ( Nsamples )

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
        self.__infl      = N.array(self.__infl)

    def sample_nonparametric ( self, Nsamples=2000 ):
        """Draw nonparametric bootstrap samples

        :Parameters:
            *Nsamples* :
                number of bootstrapsamples to be drawn
        """
        self.__bdata,self.__bestimate,dev,self.__bthres,self.__th_bias,self.__th_acc,\
                Rkd,Rpd,outl,infl = interface.bootstrap(self.data,self.estimate,Nsamples,
                        cuts=self.cuts,parametric=False,**self.model)

    def getCI ( self, cut, conf=None ):
        """Determine the confidence interval of a cut

        :Parameters:
            *cut* :
                index(!) of the cut of interest
            *conf* :
                levels of confidence (default, levels taken from the object)
        """

        if conf is None:
            conf = self.conf
        elif isinstance ( conf, float ):
            conf = [conf]
        elif isinstance ( conf, int ):
            conf = [self.conf[conf]]

        if self.__expanded:
            ci = []
            for c in conf:
                k = self.__expandedConf.index(round(c,6))
                w = max(1-2*c,1-(1-c)*2)
                w = list(self.__expandedWidths).index(round(w,6))
                ci.append(self.__expandedCI[cut,w,k])
            return N.array(ci)

        bias = self.__th_bias[cut]
        acc  = self.__th_acc[cut]

        vals = []
        for pp in conf:
            vals.append(stats.norm.cdf( bias + ( stats.norm.ppf(pp) + bias ) / (1-acc*(stats.norm.ppf(pp) + bias )) ))

        return p.prctile ( self.__bthres[:,cut], 100*N.array(vals) )

    def __repr__ ( self ):
        return "< BootstrapInference object with %d blocks and %d samples >" % ( self.data.shape[0], self.nsamples )

    def sensitivity_analysis ( self, conf=0.95, Nsamples=2000, Npoints=8, verbose=True ):
        """Perform sensitivity_analysis to obtain expanded bootstrap intervals

        Sensitivity analysis is a strategy to expand bootstrap confidence intervals. Typically,
        the expanded bootstrap confidence intervals are more realistic than the unexpanded confidence
        intervals. The function fits a gaussian kernel density estimator to the joint distribution
        of the first two parameters (those that determine the shape of the psychometric function)
        and determines a number of points on the 68% contour of this estimated density. For each of
        these points bootstrap confidence intervals are estimated. The final bootstrap confidence
        intervals are then defined as the widest of these confidence intervals.
        After calling sensitivity_analysis() the getCI() method will give the expanded BCa confidence
        intervals.

        :Parameters:
            *conf* :
                desired confidence. Note that this is the "width" of the confidence interval, not the edges.
                This is necessary because sensitivity_analysis is used to expand the widths of the confidence
                intervals. It is ambiguous how percentiles of the bootstrap distribution would be modified
                by sensitivity_analysis.
            *Nsamples* :
                number of bootstrap samples per data point
            *Npoints* :
                number of points on the contour at which to perform a new bootstrap run

        :Output:
            *expandedCI*:
                an array of expanded confidence intervals
            *expansionPoints*:
                an array of coordinates of the points at which additional bootstrap samples were drawn
        """
        if self.__expanded:
            return N.array(self.__expandedCI),N.array(self._expansionPoints)

        if self.mcestimates is None:
            # We need an initial run
            self.sample ( Nsamples )
        al,bt = self.estimate[:2]
        prm0 = N.array([al,bt])

        if isinstance ( conf, float ):
            conf = [conf]

        # Now fit the confidence with a kde
        contour = self.mcdensity
        maxcont = contour.evaluate ( N.array( [al,bt] ) )

        # Determine unexpanded CI
        self.__expandedCI = []
        for l,cut in enumerate(self.cuts):
            self.__expandedConf = []
            self.__expandedCI.append([])
            for prob in conf:
                pprob = 1.-prob
                p1,p2 = 0.5*pprob,1-0.5*pprob
                self.__expandedCI[-1].append( self.getCI ( l, (p1,p2) )-self.thres[l] )
                self.__expandedConf += [round(p1,6),round(p2,6)]
        self.__expandedCI = N.array(self.__expandedCI)
        self.__expandedWidths = N.array(conf)

        # Expand
        self._expansionPoints = []
        for k in xrange(Npoints):
            # The next point in parameter space
            phi = float(2*N.pi*k) / Npoints
            searchaxis = N.array([N.cos(phi),N.sin(phi)])
            x0,x1 = 0.,max(self.mcestimates[:,0].max()-al,self.mcestimates[:,1].max()-bt)
            def f ( prm ):
                coeffs = prm0+searchaxis*prm
                return contour.evaluate(coeffs)-0.68*maxcont
            while f(x0)*f(x1)>=0:
                x0 =  x1
                x1 *= 2
            # We use bisections to find the point on the surface that gives 0.68 the density of the maximum
            try:
                a = optimize.bisect ( f, x0, x1 )
            except RuntimeError:
                # This is definitely a Hack: a has to be between x0 and x1
                a = 0.5*(x0+x1)
            self._expansionPoints.append(a*searchaxis+prm0)
            if verbose:
                sys.stderr.write("Bootstrapping point %d ... " % (k,))
                sys.stderr.flush()

            # Perform bootstrap on this next point
            fullprm = self.estimate.copy()
            fullprm[:2] = self._expansionPoints[-1]
            fullthres = interface.diagnostics(self.data,fullprm, nafc=self.model["nafc"],sigmoid=self.model["sigmoid"],core=self.model["core"],cuts=self.cuts)[3]
            bthres,th_bias,th_acc = interface.bootstrap(self.data,fullprm,Nsamples,cuts=self.cuts,**self.model)[3:6]
            thresholdCI = []
            for l,cut in enumerate(self.cuts):
                for pp,prob in enumerate(conf):
                    # Transform confidence to upper and lower limits
                    pprob = 1.-prob
                    p1,p2 = 0.5*pprob,1-0.5*pprob

                    # Determine BCa-confidence intervals
                    prob1 = stats.norm.cdf ( th_bias[l] + ( stats.norm.ppf(p1) + th_bias[l] ) / (1-th_acc[l]*(stats.norm.ppf(p1) + th_bias[l])) )
                    prob2 = stats.norm.cdf ( th_bias[l] + ( stats.norm.ppf(p2) + th_bias[l] ) / (1-th_acc[l]*(stats.norm.ppf(p2) + th_bias[l])) )
                    thresholdCI = p.prctile ( bthres[:,l], 100*N.array([prob1,prob2]) ) - fullthres[l] # We substract the threshold here to only obtain relative values
                    # If this confidence interval is larger than the original one, we expand the CI
                    if thresholdCI[0]<self.__expandedCI[l,pp,0]:
                        self.__expandedCI[l,pp,0] = thresholdCI[0]
                        if verbose:
                            sys.stderr.write("l- ")
                    if thresholdCI[1]>self.__expandedCI[l,pp,1]:
                        self.__expandedCI[l,pp,1] = thresholdCI[1]
                        if verbose:
                            sys.stderr.write("u+ ")


            if verbose:
                sys.stderr.write("\n")
                sys.stderr.flush()

        # Now we add the threshold back to the ci
        for l,cut in enumerate(self.cuts):
            self.__expandedCI[l,pp,:] += self.thres[l]

        # Store that we have expanded the CIs
        self.__expanded = True

        return N.array(self.__expandedCI),N.array(self._expansionPoints)


    outl = property ( fget=lambda self: self.__outl, doc="A boolean vector indicating whether a block should be considered an outlier" )
    infl = property ( fget=lambda self: self.__infl, doc="A boolean vector indicating whether a block should be considered an influential observation" )
    mcestimates = property ( fget=lambda self: self.__bestimate, doc="An array of bootstrap estimates of the fitted paramters" )
    mcdeviance = property ( fget=lambda self: self.__bdeviance, doc="A vector of bootstrapped deviances" )
    mcRpd = property ( fget=lambda self: self.__bRpd, doc="A vector of correlations between model prections and deviance residuals in all bootstrap samples" )
    mcRkd = property ( fget=lambda self: self.__bRkd, doc="A vector of correlations between block index and deviance residuals in all bootstrap samples" )
    mcthres = property ( fget=lambda self: self.__bthres, doc="Thresholds of the bootstrap replications" )
    mcdensity = property ( fget=lambda self: stats.kde.gaussian_kde ( self.mcestimates[N.logical_and(self.mcestimates[:,0]<10000,self.mcestimates[:,1]<10000),:2].T ),
            doc="A gaussian kernel density estimate of the joint density of the first two parameters of the model" )
    inference = property ( fget=lambda self: "CML-MC", doc="Type of inference performed by the object" )
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
    def __init__ ( self, data, sample=True, cuts=(.25,.5,.75), conf=(.025,.975), automatic=True, resample=False, plotprm=None, **kwargs ):
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
                    - 'exp'
            *core* :
                term inside the sigmoid function. Valid choices are
                    - 'ab'         (x-a)/b        [Default]
                    - 'mw%g'       midpoint and width
                    - 'linear'     a+b*x
                    - 'log'        a+b*log(x)
                    - 'weibull'    2*s*m*(log(x)-log(m))/log(2) + log(log(2))    This will give you a weibull if combined with
                      the gumbel_l sigmoid and a reverse weibull if combined with the gumbel_r sigmoid.
                    - 'poly'       (x/a)**b   Will give you a weibull if combined with an exp sigmoid
            *priors* :
                a list of prior names. Valid choices are
                    - 'Uniform(%g,%g)'   Uniform distribution on an interval
                    - 'Gauss(%g,%g)'     Gaussian distribution with mean and standard deviation
                    - 'Beta(%g,%g)'      Beta distribution
                    - 'Gamma(%g,%g)'     Gamma distribution
                    - 'nGamma(%g,%g)'    Gamma distribution on the negative axis
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
            *plotprm* :
                a dictionary to take parameters for plotting data. Currently supported are the arguments
                'label', 'color', 'linestyle', 'linewidth' and 'marker'. These can all be set after creating
                an Inference instance, too. By using the respective properties.

        :Example:
        Use MCMC to estimate a psychometric function from some example data and derive posterior
        intervals for a threshold

        >>> x = [0,2,4,6,8,10]
        >>> k = [26,30,38,46,50,49]
        >>> n = [50]*len(k)
        >>> mcmc = BayesInference ( zip(x,k,n), priors=("Gauss(0,5)","Gamma(1,4)","Beta(2,50)") )
        >>> mcmc.sample ( start=3*mcmc.estimate )
        >>> mcmc.sample ( start=0.1*mcmc.estimate )
        >>> mcmc.estimate
        array([ 3.64159245,  5.13138577,  0.02117899])
        >>> mcmc.deviance
        3.2953368439616186
        >>> mcmc.getThres()
        3.6522408270087667
        >>> mcmc.getCI()[1]
        array([ 2.65917603,  3.68535429,  4.56688308])
        """
        PsiInference.__init__(self,plotprm)

        # Store basic data
        self.data = N.array(data)
        if self.data[:,1].max() <= 1:
            # We have relative frequencies
            self.data[:,1] *= self.data[:,2]
            self.data[:,1] = N.floor ( self.data[:,1] )
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "mw0.1"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }
        self.retry = resample

        if self.model["core"][:2] == "mw":
            self.parnames = ["m","w"]
        elif self.model["core"] == "weibull":
            self.parnames = ["m","s"]
        else:
            self.parnames = ["a","b"]
        self.parnames.append("lambda")
        if self.model["nafc"]<2:
            self.parnames.append("guess")

        self.afac = kwargs.setdefault ( "afac", 0.4 )

        self.mapestimate,self.fisher,thres,self.mapdeviance = interface.mapestimate(self.data,start=None,**self.model)

        if cuts is None:
            self.cuts = (.25,.5,.75)
        else:
            self.cuts = cuts
        if isinstance(cuts,float):
            self.Ncuts = 1
        else:
            self.Ncuts = len(self.cuts)

        self.Rpd,self.Rkd = interface.diagnostics ( self.data, self.mapestimate, cuts=self.cuts, nafc=self.model["nafc"],
                sigmoid=self.model["sigmoid"], core=self.model["core"] )[4:]

        self.__meanestimate = None
        self.__meandeviance = None

        self.__mcmc_chains                         = []
        self.__mcmc_deviances                      = []
        self.__mcmc_posterior_predictives          = []
        self.__mcmc_posterior_predictive_deviances = []
        self.__mcmc_posterior_predictive_Rpd       = []
        self.__mcmc_posterior_predictive_Rkd       = []
        self.__mcmc_logposterior_ratios            = []

        self.__pRpd   = None
        self.__pRkd   = None
        self.__pthres = None

        self.conf = conf

        self.burnin = 0
        self.thin   = 1
        self.nsamples = None

        # We assume that parameter variation is proportional to the
        # estimated parameters
        self._steps = 0.1*self.mapestimate * 500/N.sum(self.data[:,2])
        # print self._steps

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
        chain,deviance,ppdata,ppdeviances,ppRpd,ppRkd,logpostratios = interface.mcmc ( self.data, start, Nsamples, stepwidths=self._steps, **self.model )
        self.__mcmc_chains.append(N.array(chain))
        self.__mcmc_deviances.append(N.array(deviance))
        self.__mcmc_posterior_predictives.append(N.array(ppdata))
        self.__mcmc_posterior_predictive_deviances.append(N.array(ppdeviances))
        self.__mcmc_posterior_predictive_Rpd.append (N.array(ppRpd))
        self.__mcmc_posterior_predictive_Rkd.append (N.array(ppRkd))
        self.__mcmc_logposterior_ratios.append (N.array(logpostratios) )

        self.__recomputeCorrelationsAndThresholds()

        # Resample if bad
        if self.retry:
            nprm = self.mapestimate.shape[0]
            resampled = 0
            while True:
                if resampled > 5:
                    raise SamplingError, "Resampling did not yield a converging chain after 5 tries"
                for k in xrange ( nprm ):
                    allok = True
                    ok,z,bad = self.geweke ( k )
                    if not bad is None:
                        # print "Resampling"
                        for b in bad:
                            self.resample ( b )
                        resampled += 1
                        break
                    else:
                        allok = allok and ok
                if allok==True:
                    break

    def resample ( self, chain, Nsamples=None, start=None ):
        """Replace a chain

        :Parameters:
            *chain* :
                index of the chain to be replaced
            *Nsamples* :
                number of posterior samples to be drawn. Default is to take the same number of samples
            *start* :
                starting value for the chain. Default is to take the same starting value as for the previous
                chain. If an integer is given as the starting value, this will be interpreted as the position
                of the old chain at which the new starting value can be found.
        """
        if isinstance (Nsamples,int):
            self.nsamples = Nsamples
        elif Nsamples is None:
            Nsamples = self.__mcmc_chains[chain].shape[0]
        else:
            Nsamples = self.nsamples

        if start is None:
            start = self.__mcmc_chains[chain][0,:]
        elif isinstance (Nsamples,int):
            start = self.__mcmc_chains[chain][start,:]

        mcchain,deviance,ppdata,ppdeviances,ppRpd,ppRkd,logpostratios = interface.mcmc ( self.data, start, Nsamples, stepwidths=self._steps, **self.model )
        self.__mcmc_chains[chain] = N.array(mcchain)
        self.__mcmc_deviances[chain] = N.array(deviance)
        self.__mcmc_posterior_predictives[chain] = N.array(ppdata)
        self.__mcmc_posterior_predictive_deviances[chain] = N.array(ppdeviances)
        self.__mcmc_posterior_predictive_Rpd[chain] = N.array(ppRpd)
        self.__mcmc_posterior_predictive_Rkd[chain] = N.array(ppRkd)
        self.__mcmc_logposterior_ratios[chain] = N.array(logpostratios)

        self.__recomputeCorrelationsAndThresholds()

    def bayesian_p ( self, quantity="deviance" ):
        """Bayesian p value associated with a given quantity

        The Bayesian p value of a model compares posterior predictives with the observed data.
        If the observed data are very unequal to the posterior predictives, this indicates that
        the model does not describe the data well. To compare observed data and simulated data
        (posterior predictives), it is common to derive a quantity of interest from the posterior
        predictives. The Bayesian p value is between 0 and 1 and values close to 0 and close to 1
        indicate bad model fit. This p value can be interpreted like a two sided test.

        :Parameters:
            *quantity* :
                This is the quantity do be derived. By default only deviance, Rpd and Rkd are available.
                If quantity is a function, this will be called on every data set and the respective
                p value will be calculated. The call on every data set takes two arguments:
                1. a nblocksX3 array of data and
                2. a parameter vector.
                This way any other transformation of the data can be realized.

        :Output:
            the bayesian p-value
        """
        if isinstance ( quantity, str ):
            if quantity.lower() == "deviance":
                return N.mean ( (self.ppdeviance-self.mcdeviance)>=0 )
            elif quantity.lower() == "rpd":
                return N.mean ( (self.ppRpd-self.mcRpd)>=0 )
            elif quantity.lower() == "rkd":
                return N.mean ( (self.ppRkd-self.mcRkd)>=0 )
            else:
                raise ValueError, "unsupported quantity for bayesian p value"
        elif operator.isCallable ( quantity ):
            d = self.data.copy()
            I = 0.
            for k in xrange ( self.Nsamples ):
                d[:,1] = self.posterior_predictive[k,:]
                I += double ( quantity ( d, self.mcestimates[k,:] ) - quantity ( self.data, self.mcestimates[k,:] >= 0 ) )
            return I/self.Nsamples

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

    def getmcdeviance ( self, chain=None, raw=False ):
        """Get samples from the posterior distribution of deviances

        :Parameters:
            *chain* :
                if chain is None, the samples are combined across all chains
                sampled so far. If chain is an integer, it is interpreted as
                the index of the chain to be returned
            *raw* :
                is true if deviances for all samples are to be returned (not
                respecting burnin and thinning). This only has an effect for
                single chains.

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
            if raw:
                return self.__mcmc_deviances[chain]
            else:
                return self.__mcmc_deviances[chain][self.burnin::self.thin]
        else:
            raise ValueError, "chain should be either None or an integer"

    def getppdata ( self, chain=None, raw=False ):
        """Get posterior predictive data

        Posterior predictive data are data samples from the joint posterior over
        data and parameters. These represent data that could be generated by the
        model. Comparison of posterior predictive data and the observed data forms
        the basis of bayesian model checking: If posterior predictive data differ
        systematically from the observed data, the fitted model does not capture
        all the structure in the data.

        :Parameters:
            *chain* :
                chain for which posterior predictive data should be returned
            *raw* :
                is true if all data (not respecting burnin and thinning) are to
                be returned (this only has an effect for single chains!)

        :Output:
            a array of nsamplesXncorrect predicted data
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_posterior_predictives:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance ( chain, int ):
            # Get a single chain
            if raw:
                # Get raw data
                return self.__mcmc_posterior_predictives[chain]
            else:
                return self.__mcmc_posterior_predictives[chain][self.burnin::self.thin]
        else:
            raise IndexError, "chain should be either None or an integer"

    def getppdeviance ( self, chain=None, raw=False ):
        """Get deviances associated with posterior predictive data

        Posterior predictive data are data samples from the joint posterior over data
        and parameters. Deviance of these samples is one possible transformation on
        which a comparison of these data with the observed data could be based.

        :Parameters:
            *chain* :
                chain index to be returned. If this is None (the default) all chains
                are combined.
            *raw* :
                if only a single chain is returned, it might be interesting to see
                the whole chain and ignore burnin and thinning. If raw==True, the
                chain is requested in this "raw" format.

        :Output:
            an array of nsamples deviances
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_posterior_predictive_deviances:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance ( chain, int ):
            if raw:
                return self.__mcmc_posterior_predictive_deviances[chain]
            else:
                return self.__mcmc_posterior_predictive_deviances[chain][self.burnin::self.thin]
        else:
            raise ValueError, "chain should be either None or an integer"

    def getppRpd ( self, chain=None, raw=False ):
        """Get correlations between psychometric function and deviance residuals associated with posterior predictive data

        Posterior predictive data are data samples from the joint posterior over data
        and parameters. Correlation between the psychometric function and the
        deviance residuals of these samples is one possible transformation on
        which a comparison of these data with the observed data could be based.

        :Parameters:
            *chain* :
                chain index to be returned. If this is None (the default) all chains
                are combined.
            *raw* :
                if only a single chain is returned, it might be interesting to see
                the whole chain and ignore burnin and thinning. If raw==True, the
                chain is requested in this "raw" format.

        :Output:
            an array of nsamples deviances
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_posterior_predictive_Rpd:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance ( chain, int ):
            if raw:
                return self.__mcmc_posterior_predictive_Rpd[chain]
            else:
                return self.__mcmc_posterior_predictive_Rpd[chain][self.burnin::self.thin]
        else:
            raise ValueError, "chain should be either None or an integer"

    def getppRkd ( self, chain=None, raw=False ):
        """Get correlations between block index and deviance residuals associated with posterior predictive data

        Posterior predictive data are data samples from the joint posterior over data
        and parameters. Correlation between the block index and the deviance residuals
        of these samples is one possible transformation on which a comparison of these
        data with the observed data could be based.

        :Parameters:
            *chain* :
                chain index to be returned. If this is None (the default) all chains
                are combined.
            *raw* :
                if only a single chain is returned, it might be interesting to see
                the whole chain and ignore burnin and thinning. If raw==True, the
                chain is requested in this "raw" format.

        :Output:
            an array of nsamples deviances
        """
        if chain==None:
            # Get all chains
            chains = []
            for chain in self.__mcmc_posterior_predictive_Rkd:
                chains.append ( chain[self.burnin::self.thin] )
            return N.concatenate ( chains, 0 )
        elif isinstance ( chain, int ):
            if raw:
                return self.__mcmc_posterior_predictive_Rkd[chain]
            else:
                return self.__mcmc_posterior_predictive_Rkd[chain][self.burnin::self.thin]
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
        indices = N.random.randint ( samples.shape[0], size=(Nsamples,) )
        # Scale deviance to 0,1
        deviances -= deviances[indices].min()
        deviances /= deviances[indices].max()
        deviances = N.clip(.4+4*deviances,0,1)
        for i in indices:
            # i = N.random.randint(samples.shape[0])
            psi = N.array(interface.diagnostics ( x, samples[i,:], sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
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
            return False,z,bad
        else:
            return True,z,None

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
    inference = property ( fget=lambda self: "MCMC", doc="Type of inference performed by the object" )
    mcthres = property ( fget=lambda self: self.__pthres, doc="posterior samples of the threshold" )
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
                    self.devianceresiduals,self.__meandeviance,self.thres,self.Rpd,self.Rkd = interface.diagnostics ( \
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
    posterior_predictive = property ( fget=getppdata, doc="Posterior predictive data associated with the MCMC samples" )
    ppdeviance = property ( fget=getppdeviance, doc="Deviances associated with the posterior predictive data" )
    ppRpd = property ( fget=getppRpd, doc="Correlations between psychometric function and deviance residuals associated with posterior predictive data" )
    ppRkd = property ( fget=getppRkd, doc="Correlations between block index and deviance residuals associated with posterior predictive data" )

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
    def nullevidence ():
        """model evidence for the corresponding null model

        This can be used for model selection: model evidence devided by null evidence gives the Bayes Factor
        for the comparison of the model agains the null model. This can be interpreted as "how much more
        probable is the given psychometric function than the null model for the present data. Also see the
        documentation for the evidence property.
        """
        def fget (self):
            # The null deviance can be directly calculated
            n = self.data[:,2].sum()
            k = self.data[:,1].sum()
            alpha,beta = 1.,1.    # flat prior for the null model
            return float(int(1./special.beta(k,n-k)))*special.beta(k+alpha,n-k+beta)/special.beta(alpha,beta)

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

    @Property
    def farstart ():
        """A proper starting value for the Rhat statistic

        This is a starting value for the mcmc process, that is relatively far away from the posterior density.
        In order to have a reasonably interpretable Rhat statistic. There should be multiple chains and these chains
        should have overdispersed starting values. farstart will always correspond to an overdispersed starting value.
        """
        def fget ( self ):
            k = N.random.randint(2)
            l = N.random.randint(2)
            x = self.mapestimate
            x[l] = p.prctile ( self.mcestimates[:,l], (2.5,97.5)[k] )
            # print x
            return x

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
        self._PsiInference__infl   = N.zeros(self.data.shape[0], 'd' )

        for k,theta in enumerate(samples):
            self.__pthres[k,:],self.__pRpd[k],self.__pRkd[k] = interface.diagnostics (\
                    self.data, theta, cuts=self.cuts, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )[3:]
        lpr = []
        for l in self.__mcmc_logposterior_ratios:
            lpr.append(l[self.burnin::self.thin,:])
        lpr = N.concatenate ( lpr, 0 )
        self._PsiInference__infl = -N.mean(lpr,0) + N.log(N.mean(N.exp(lpr),0))

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

        verbose = True

        # Determine size of initial test run
        if self.nsamples is None:
            NN = 0
            for q in self.conf:
                Nmin = pygibbsit.gibbsit ( q=q )["Nmin"]
                NN = max(NN,Nmin)
            self.nsamples = NN

        a = self.__roughvariance ()
        # asympvar = N.diag(fisherinv)
        # a = self.afac*N.sqrt(asympvar)
        # print a

        # chain,deviance,ppdata,ppdeviances,ppRpd,ppRkd,logpostratios = interface.mcmc ( self.data, self.mapestimate, NN, stepwidths=a, **self.model )
        # a = N.sqrt(N.diag(N.cov(chain.T)))
        # print a

        oldburnin = 0
        oldthin   = 1
        oldnsamples = NN
        for n in xrange ( noptimizations ):
            samples,deviances,ppdata,ppdeviances,ppRpd,ppRkd,logpostratios = interface.mcmc ( self.data, self.mapestimate, NN, stepwidths=a, **self.model )

            # Check all desired thresholds
            for q in self.conf:
                for k in xrange ( len(self.mapestimate) ):
                    try:
                        mcmcpars = pygibbsit.gibbsit ( samples[:,k], q=q )
                    except IndexError:
                        continue
                    self.burnin = max ( self.burnin, mcmcpars.burnin )
                    self.thin   = max ( self.thin,   mcmcpars.thin )
                    self.nsamples = max ( self.nsamples, mcmcpars.Nsamples )
            b = N.sqrt(N.diag(N.cov ( samples[self.burnin::self.thin].T )))
            if b.max() < 1e-10:
                a = self.__roughvariance ()
            else:
                a = b

            if verbose:
                print "Burnin:",self.burnin,"Thinning:",self.thin,"Nsamples:",self.nsamples
                print "Steps:",a
            if self.nsamples <= oldnsamples:
            # if oldburnin==self.burnin and oldthin==self.thin and oldnsamples==self.nsamples:
                break
            else:
                oldburnin,oldthin,oldnsamples = self.burnin,self.thin,self.nsamples
        self.mcmcpars = mcmcpars
        self._steps = a
        if verbose:
            print "Steps(final):",N.sqrt(N.diag(N.cov( samples[self.burnin::self.thin].T )))

    def __roughvariance ( self ):
        # Determine an initial variance estimate using the Fisher Information Matrix
        fisherI = N.matrix(self.fisher)
        try:
            fisherIinv = N.linalg.solve ( fisherI.T*fisherI+0.01*N.eye(fisherI.shape[0]), fisherI.T )
        except N.linalg.LinAlgError:
            # It seems as if the regularized fisher matrix can not be inverted
            # We directly get an estimate form bootstrap
            bsamples = interface.bootstrap ( self.data, self.estimate, 100, cuts=self.cuts, **self.model )[1]
            return bsamples.std(0)
        cond = abs(fisherI.A).sum(1).max() * abs(fisherIinv.A).sum(1).max()
        # print "Condition of Fisher Information Matrix:",cond
        # print fisherI

        if cond > 1e6:
            for k in xrange(20):
                localdata = self.data.copy()
                for l in xrange ( localdata.shape[0] ):
                    localdata[l,1] = N.random.binomial ( self.data[l,2], self.data[l,1].astype('d')/self.data[l,2] )
                fisherII = N.matrix(interface.mapestimate(localdata,start=None,**self.model)[1])
                fisherIIinv = N.linalg.solve ( fisherII.T*fisherII+0.01*N.eye(fisherII.shape[0]), fisherII.T )
                cond = abs(fisherII.A).sum(1).max() * abs(fisherIIinv.A).sum(1).max()
                # print "Condition of Fisher Information Matrix:",cond
                # print fisherI
                if cond < 1e6:
                    fisherI = fisherII
                    fisherIinv = fisherIIinv
                    break

        a = N.sqrt(N.diag(fisherIinv))
        # print "a =",a

        if abs(a).min() < 1e-10 or abs(a).max() > 1e10 or a[2] > 0.5:
            # It seems as if the Variance estimation via the Fisher Matrix failed
            bsamples = interface.bootstrap(self.data,self.estimate,100,cuts=self.cuts,**self.model)[1]
            a = bsamples.std(0)
            # print "a_boots =",a

        return a

if __name__ == "__main__":
    import doctest
    doctest.testmod()
