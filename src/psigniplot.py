#/usr/bin/env python

import pylab as p
import numpy as N

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

def plotPMF ( InferenceObject, ax=None ):
    raise NotImplementedError()

def plotGeweke ( BayesInferenceObject, ax=None ):
    raise NotImplementedError()

def plotChains ( BayesInferenceObject, ax=None ):
    raise NotImplementedError()

def GoodnessOfFit ( InferenceObject ):
    raise NotImplementedError()

gof = GoodnessOfFit
