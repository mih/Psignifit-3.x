#/usr/bin/env python

import pylab as p
import numpy as N
from pypsignifit import BayesInference, BootstrapInference

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
#    def pmfanddata ( self, ax=None, xlabel_text="Stimulus intensity", ylabel_text=None ):
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
    """Plot thresholds and confidence intervals"""
    if ax == None:
        ax = p.axes()

    for k,cut in enumerate(InferenceObject.cuts):
        c25,c975 = InferenceObject.getCI ( cut=k, conf=(.025,.975) )
        thres = InferenceObject.getThres ( cut )
        ylev = InferenceObject.evaluate ( [thres] )
        # ylev  = _psipy.diagnostics ( [thres],   self.estimate, cuts=cut, nafc=self.model["nafc"], sigmoid=self.model["sigmoid"], core=self.model["core"] )
        ax.plot ( [c25,thres,c975],[ylev]*3, 'b-|' )

def GoodnessOfFit ( InferenceObject ):
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
    ax = p.axes([0,.5,.33,.5] )
    # if isinstance ( InferenceObject, BayesInference ):
    try:
        InferenceObject.drawposteriorexamples ( ax=ax )
    except:
        pass
    plotThres ( InferenceObject, ax=ax )
    plotPMF   ( InferenceObject, ax=ax )
    plotHistogram ( InferenceObject.mcdeviance, InferenceObject.deviance, "posterior deviance", "D", p.axes ( [0,0,.33,.5] ) )

    # Second part: Correlations between model prediction and residuals
    plotRd ( InferenceObject, p.axes([.33,.5,.33,.5]), "p" )
    good = plotHistogram ( InferenceObject.mcRpd, InferenceObject.Rpd, "posterior Rpd", "Rpd", p.axes([.33,0,.33,.5]) )
    if not good and warn==True:
        if isinstance ( InferenceObject, BootstrapInference ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean() , "Simulated Rpd differs from observed!\nModel deviates systematically from data", \
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )
        elif isinstance ( InferenceObject, BayesInferenceObject ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean() , "Rpd is different from 0!\nModel deviates systematically from data", \
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    # Third part: Correlations between model prediction and block index
    plotRd ( InferenceObject, p.axes([.66,.5,.33,.5]), "k" )
    good = plotHistogram ( InferenceObject.mcRkd, InferenceObject.Rkd, "posterior Rkd", "Rkd", p.axes([.66,0,.33,.5]) )
    if not good and warn==True:
        if isinstance ( InferenceObject, BootstrapInference ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean(), "Simulated Rkd differs from observed!\nData are nonstationary!",\
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )
        elif isinstance ( InferenceObject, BayesInference ):
            p.text ( 0, p.getp(p.gca(),'ylim').mean(), "Rkd is different from 0!\nData are nonstationary!",\
                    fontsize=16, color=warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

def plotGeweke ( BayesInferenceObject, ax=None ):
    raise NotImplementedError()

def plotChains ( BayesInferenceObject, ax=None ):
    raise NotImplementedError()

gof = GoodnessOfFit
