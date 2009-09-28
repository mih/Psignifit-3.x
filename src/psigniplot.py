#/usr/bin/env python

import pylab as p
import numpy as N
import pypsignifit
import re
from scipy import stats

__warnred = [.7,0,0]

def drawaxes ( ax, xtics, xfmt, ytics, yfmt, xname, yname ):
    """Draw x and y axes that look nicer than standard matplotlib

    This function deletes the pylab default axes and draws visually more
    pleasing axes. In addition it takes care that all graphics remain within
    the area covered by the axes.

    :Parameters:
        ax      the axes object in which the graphics should stay
        xtics   an array of x-ticks
        xfmt    format string for the x-ticks
        ytics   an array of y-ticks
        yfmt    format string for the y-ticks
        xname   label for the x-axis
        yname   label for the y-axis
    """

    # Data ranges
    yr = ytics.max()-ytics.min()
    xr = xtics.max()-xtics.min()

    # x axis
    yt = ytics.min()
    ax.plot ( [xtics.min(),xtics.max()], [yt-.05*yr]*2, 'k-' )
    for xt in xtics:
        ax.plot ( [xt]*2, [yt-.05*yr,yt-.03*yr], 'k-' )
        ax.text ( xt, yt-.06*yr, xfmt % (xt,), horizontalalignment="center", verticalalignment="top",fontsize=10 )
    ax.text ( xtics.mean(), yt-.12*yr, xname, horizontalalignment="center", verticalalignment="top", fontsize=16 )

    # y axis
    xt = xtics.min()
    ax.plot ( [xt-.05*xr]*2, [ytics.min(),ytics.max()] ,'k-' )
    for yt in ytics:
        ax.plot ( [xt-.05*xr,xt-.03*xr], [yt]*2, 'k-' )
        ax.text ( xt-.06*xr, yt, yfmt % (yt,), horizontalalignment="right", verticalalignment="center", fontsize=10 )
    ax.text ( xt-.18*xr, ytics.mean(), yname, horizontalalignment="right", verticalalignment="center", fontsize=16, rotation=90 )

    # Delete all the ugly parts of the axis
    p.setp(ax, frame_on=False,\
            xticks=(), xlim=(xtics.min()-.3*xr,xtics.max()+.1*xr), \
            yticks=(), ylim=(ytics.min()-.3*yr,ytics.max()+.1*yr) )

def plotRd ( InferenceObject, ax=None, regressor="p" ):
    """plot deviance residuals against a regressor

    Deviance residuals are used plotted agains either predicted performance or
    block index to check for systematic deviations of the data from the fitted
    function.

    :Parameters:
        InferenceObject a BootstrapInference or BayesInference object containing
                    the actual inference data
        ax          an axes object where the plot should go
        regressor   plot deviance residuals against model prediction (p) or
                    against block index (k)
    """
    if ax==None:
        ax = p.axes()

    # Plot the data points
    if regressor=="p":
        psi = InferenceObject.evaluate ( InferenceObject.data[:,0] )
    elif regressor=="k":
        psi = N.arange(len(InferenceObject.data[:,0]))
    else:
        raise ValueError,"regressor %s is unknown" % regressor
    psilims = N.array([psi.min(),psi.max()])
    devianceresiduals = InferenceObject.devianceresiduals
    ax.plot ( psi, devianceresiduals, "bo" )

    # Linear regression
    A = N.ones((len(psi),2),'d')
    A[:,1] = psi
    a,b = N.linalg.lstsq(A,devianceresiduals)[0]
    ax.plot(psilims,a+b*psilims,'b:')

    if regressor=="p":
        if InferenceObject.model["nafc"]==1:
            p.setp(ax,xlim=(0,1))
        else:
            p.setp(ax,xlim=(1./InferenceObject.model["nafc"],1))
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
        ax.text(psilims.mean(),ytics[-2],"Rpd=%.3f" % ( InferenceObject.Rpd, ) )
        xname = "model prediction"
    elif regressor=="k":
        ax.text(psilims.mean(),ytics[-2],"Rkd=%.3f" % ( InferenceObject.Rkd, ) )
        xname = "block index"

    drawaxes ( ax, xtics, "%g", ytics, "%g", xname, "deviance residuals" )

def plotHistogram ( simdata, observed, xname, shortname=None, ax=None, hideobserved=False ):
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
    drawaxes ( ax, xtics, "%g", ytics, "%d", xname, "number per bin" )

    # Write diagnostics
    yt = ytics.max()
    ax.text ( xtics.min(), yt+.1, "%s=%.3f, c(2.5%%)=%.3f, c(97.5%%)=%.3f" % (shortname,observed,p25,p975),\
            horizontalalignment="left",verticalalignment="bottom", fontsize=8 )

    if observed>p25 and observed<p975:
        return True
    else:
        return False

def plotPMF ( InferenceObject, xlabel_text="Stimulus intensity", ylabel_text=None,ax=None ):
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
    xmin = InferenceObject.data[:,0].min()
    xmax = InferenceObject.data[:,0].max()
    x = N.mgrid[xmin:xmax:100j]
    # psi = N.array(_psipy.diagnostics ( x, self.estimate, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
    psi = InferenceObject.evaluate ( x )
    ax.plot(x,psi,'b')

    # Plot the data
    xd = InferenceObject.data[:,0]
    pd = InferenceObject.data[:,1].astype("d")/InferenceObject.data[:,2]
    goodpoints = N.ones(pd.shape,bool)
    if not InferenceObject.outl==None:
        ax.plot(xd[InferenceObject.outl],pd[InferenceObject.outl],'r^')
        goodpoints = N.logical_and(goodpoints,N.logical_not(InferenceObject.outl))
    if not InferenceObject.infl==None:
        ax.plot(xd[InferenceObject.infl],pd[InferenceObject.infl],"rs")
        goodpoints = N.logical_and(goodpoints,N.logical_not(InferenceObject.infl))
    ax.plot(xd[goodpoints],pd[goodpoints],'bo')

    # Check axes limits
    if InferenceObject.model["nafc"]>1:
        ymin,ymax = 1./InferenceObject.model["nafc"]-.05,1.05
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
    if InferenceObject.model["nafc"]==1:
        for k,yt in enumerate(ytics):
            if yt<0 or yt>1:
                ytics.pop(k)
    else:
        for k,yt in enumerate(ytics):
            if yt<(1./InferenceObject.model["nafc"]) or yt>1:
                ytics.pop(k)
    ytics = N.array(ytics)

    drawaxes ( ax, xtics, "%g", ytics, "%g", xlabel_text, ylabel_text )

    # Write some model information
    if not InferenceObject.deviance is None:
        ax.text(0.5*(xmin+xmax),ymin+.05,"D=%g" % ( InferenceObject.deviance, ) )
    ax.text ( 0.5*(xmin+xmax),ymin+.1,InferenceObject.desc )

def plotThres ( InferenceObject, ax=None ):
    """Plot thresholds and confidence intervals

    :Parameters:
        InferenceObject either a BootstrapInference object or a BayesInference object
        ax              a pylab.axes object to be used for the plot.
    """
    if ax == None:
        ax = p.axes()

    for k,cut in enumerate(InferenceObject.cuts):
        c25,c975 = InferenceObject.getCI ( cut=k, conf=(.025,.975) )
        thres = InferenceObject.getThres ( cut )
        ylev = InferenceObject.evaluate ( [thres] )
        ax.plot ( [c25,thres,c975],[ylev]*3, 'b-|' )

def GoodnessOfFit ( InferenceObject, warn=True ):
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
        the mean of the posterior for BayesInference objects and the (constrained)
        maximum likelihood fit for BootstrapInference objects. To get an idea of the posterior
        resp. bootstrap distribution, credibility intervals are plotted at some positions (the
        location and width of the credibility intervals is given in the constructor). To make
        the posterior distribution for BayesInference objects really "plastic", a number of
        samples from the posterior distribution over psychometric functions are also drawn in
        light blue. The saturation of blue also codes the deviance of the respective function:
        the more saturated, the better the fit. For BootstrapInference objects, outliers and
        influential observations are marked as red triangles and red squares.
    2.  A histogram to approximate the posterior resp. bootstrap distribution of deviances.
        For BootstrapInference objects this histogram provides important information. It estimates
        the distribution of deviance that would be expected if the fitted model were perfectly
        valid. If the deviance of the fitted model is far in the tails of the deviance histogram,
        this typically indicates a bad fit. In that case, a warning is displayed if warn==True.
    3.  A plot of model predictions against deviance residuals. If there is no obvious
        interrelation between model prediction and deviance residuals, this indicates that the
        model describes the data reasonably well. To get an idea of the interrelation between
        model prediction and deviance residuals, the best fitting line is plotted as a dotted line.
    4.  A histogram of samples from the posterior resp. bootstrap distribution of correlations
        between model prediction and deviance residuals. The interpretation of this histogram
        differs for BootstrapInference and for BayesInference. For BayesInference the distibution
        should include 0. If the distribution is clearly shifted away from 0, this is strong
        evidence, that something is wrong with your model or your data. For BootstrapInference,
        The distribution shown corresponds to the distribution that would be expected if your
        fitted psychometric function would perfectly describe the data. Thus, if the maximum
        likelihood estimate (the vertical bold red line) is in the extremes of the distribution,
        this is strong evidence, that something is wrong with your model or your data.
    5,6 Similar to 3 and 4 but form correlations between block index and deviance residuals.
        Correlations between block index and deviance residuals indicate nonstationary
        data as should be found during e.g. perceptual learning.

    :Parameters:
        warn    if warn is set to True, red warning messages are displayed
                whenever the fit does not seem to describe the data well.
    """
    p.figure(figsize=(10,8))

    # First part: Data and fitted function, bottom deviance
    ax = p.axes([0,.5,.33,.5] )
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        InferenceObject.drawposteriorexamples ( ax=ax )
    plotThres ( InferenceObject, ax=ax )
    plotPMF   ( InferenceObject, ax=ax )
    good = plotHistogram ( InferenceObject.mcdeviance, InferenceObject.deviance, "posterior deviance", "D", p.axes ( [0,0,.33,.5] ) )
    if warn and not good:
        p.text ( p.gca().get.xticks().mean(), p.gca().get.yticks().mean(),
                "The fitted model is a bad\ndescription of the data!",
                fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    # Second part: Correlations between model prediction and residuals
    plotRd ( InferenceObject, p.axes([.33,.5,.33,.5]), "p" )
    good = plotHistogram ( InferenceObject.mcRpd, InferenceObject.Rpd, "posterior Rpd", "Rpd", p.axes([.33,0,.33,.5]) )
    if not good and warn==True:
        if InferenceObject.__repr__().split()[1] == "BootstrapInference":
            p.text ( 0, p.getp(p.gca(),'ylim').mean() , "Simulated Rpd differs from observed!\nModel deviates systematically from data", \
                    fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )
        elif InferenceObject.__repr__().split()[1] == "BayesInferenceObject":
            p.text ( 0, p.getp(p.gca(),'ylim').mean() , "Rpd is different from 0!\nModel deviates systematically from data", \
                    fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    # Third part: Correlations between model prediction and block index
    plotRd ( InferenceObject, p.axes([.66,.5,.33,.5]), "k" )
    good = plotHistogram ( InferenceObject.mcRkd, InferenceObject.Rkd, "posterior Rkd", "Rkd", p.axes([.66,0,.33,.5]) )
    if not good and warn==True:
        if isinstance ( InferenceObject, BootstrapInference ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean(), "Simulated Rkd differs from observed!\nData are nonstationary!",\
                    fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )
        elif isinstance ( InferenceObject, BayesInference ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean(), "Rkd is different from 0!\nData are nonstationary!",\
                    fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

def plotGeweke ( BayesInferenceObject, parameter=0, ax=None, warn=True ):
    """Geweke plot of moving average of samples

    :Parameters:
        BayesInferenceObject    a BayesInference object that contains all the
                    infromation about the sampling process
        parameter   index of the model parameter of interest
        ax          the pylab.axes object where the plot should go
        warn        should a warning message be displayed if non stationarity
                    of the samples is observed?
    """
    stationary,z = BayesInferenceObject.geweke ( parameter )

    if ax is None:
        ax = p.axes()

    for k in xrange ( z.shape[-1] ):
        p.plot(z[:,k],'o-')
    xtics = N.array(p.getp(ax,"xticks"))
    p.setp(ax,"xticks",xtics,"yticks",N.array((-3,-2,-1,0,1,2,3)))
    p.plot ( [xtics.min(),xtics.max()],[-2]*2,'k:')
    p.plot ( [xtics.min(),xtics.max()],[ 2]*2,'k:')
    drawaxes ( ax, xtics, "%g", N.array((-3,-2,-1,0,1,2,3)), "%g", "chain segment", "z-score" )

    if warn and not stationary:
        p.text(0.5*nsegments,0,"chains did not converge" , color=warnred, fontsize=16, rotation=45, verticalalignment="center", horizontalalignment="center" )

def plotChains ( BayesInferenceObject, parameter=0, ax=None, raw=False, warn=True ):
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
    for c in xrange(BayesInferenceObject.nchains):
        samples = BayesInferenceObject.getsamples ( c, raw=True )
        p.plot ( samples[:,parameter] )

    # Learn something about the axes
    xtics = N.array(ax.get_xticks())
    x0    = xtics.min()
    xr    = xtics.max()-xtics.min()
    ytics = ax.get_yticks()
    y0    = ytics.min()
    yr    = N.array(ytics.max()-ytics.min())

    p.text(x0+0.6*xr,y0+0.95*yr,"R^ = %.4f" % (BayesInferenceObject.Rhat ( parameter ) ) )

    if warn and BayesInferenceObject.Rhat(parameter)>1.1:
        p.text(x0+0.5*xr,y0+0.5*yr,"Chains do not seem to sample\nfrom the same distribution!",
                horizontalalignment="center",verticalalignment="center",fontsize=16,rotation=45,color=warnred)

    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", "sample #", BayesInferenceObject.parnames[parameter] )

def plotParameterDist ( InferenceObject, parameter=0, ax=None ):
    """Plot the distribution of parameters

    :Parameters:
        InferenceObject either a BootstrapInference object or a BayesInference object
                    containing the samples of the parameter distribtution
        parameter   index of the model parameter of interest
        ax          pylab.axes object where the plot should go
    """
    if ax is None:
        ax = p.axes()

    samples = InferenceObject.mcestimates[:,parameter]
    h,b,ptch = p.hist ( samples, bins=20, normed=True, histtype="step", lw=2 )

    if InferenceObject.__repr__().split()[1] == "BayesInference":
        priorstr = InferenceObject.model["priors"]
        if not priorstr is None:
            priorstr = priorstr[parameter]
            m = re.search (
                r"(\w+)\((-?\d*\.?\d*[eE]?-?\d*),(-?\d*\.?\d*[eE]?-?\d*)\)",
                priorstr )
            if not m is None:
                dist,prm1,prm2 = m.groups()
                prm1,prm2 = float(prm1),float(prm2)
                x = N.mgrid[b.min():b.max():100j]
                if dist.lower () == "gauss":
                    p.plot(x,stats.norm.pdf(x,prm1,prm2))
                elif dist.lower () == "beta":
                    p.plot(x,stats.beta.pdf(x,prm1,prm2))
                elif dist.lower () == "gamma":
                    p.plot(x,stats.gamma.pdf(x,prm2,scale=prm1))
                elif dist.lower () == "uniform":
                    p.plot(x,stats.uniform.pdf(x,prm1,prm2))

    # Store ticks
    xtics = ax.get_xticks()
    ytics = ax.get_yticks()

    # Highlight estimate and credibility intervals
    prm = InferenceObject.estimate[parameter]
    c25,c975 = p.prctile ( samples, (2.5,97.5) )
    ym = ax.get_yticks().max()
    p.plot ( [c25]*2,[0,ym],'b:', [c975]*2,[0,ym],'b:' )
    p.plot ( [prm]*2,[0,ym],'b' )
    p.text ( ax.get_xticks().mean(), ym, "%s^=%.3f, CI(95%%)=(%.3f,%.3f)" % ( InferenceObject.parnames[parameter],prm,c25,c975 ),
            fontsize=8, horizontalalignment="center",verticalalignment="bottom" )

    drawaxes ( ax, xtics, "%g", ytics, "%g", InferenceObject.parnames[parameter], "density estimate" )

def ParameterPlot ( InferenceObject ):
    """Show distributions and estimates for all parameters in the model

    :Parameters:
        InferenceObject a BootstrapInference or BayesInference object containing the
                desired data
    """
    nparams = len(InferenceObject.parnames)
    axw = 1./nparams
    fig = p.figure (figsize=(3*nparams,3))

    for k in xrange ( nparams ):
        ax = p.axes ( [axw*k,0,axw,1] )
        plotParameterDist ( InferenceObject, k, ax )

def ConvergenceMCMC ( BayesInferenceObject, parameter=0, warn=True ):
    """Diagram to check convergence of MCMC chains for a single parameter

    :Parameters:
        BayesInferenceObject    a BayesInference object containing all information about
                    the model and the posterior distribution
        parameter   model parameter of interest. So far, no model derived parameters such as
                    thresholds are supported
        warn        should warnings be displayed if the samples look suspicious?
    """
    fig = p.figure ( figsize=[9,3] )
    ax =  p.axes ( [0,0.,0.33,1] )
    plotChains ( BayesInferenceObject, parameter, ax, warn=warn )
    ax = p.axes ( [.33,0,.33,1] )
    plotGeweke ( BayesInferenceObject, parameter, ax, warn=warn )
    ax = p.axes ( [.66,0,.33,1] )
    plotParameterDist ( BayesInferenceObject, parameter, ax )

gof = GoodnessOfFit
