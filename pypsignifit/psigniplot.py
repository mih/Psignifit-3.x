#/usr/bin/env python

import pylab as p
import numpy as N
import re
from scipy import stats
import psignidata

__all__ = ["GoodnessOfFit","ConvergenceMCMC","ParameterPlot","ThresholdPlot","plotSensitivity","plotInfluential"]
__warnred = [.7,0,0]

def drawaxes ( ax, xtics, xfmt, ytics, yfmt, xname, yname ):
    """Draw x and y axes that look nicer than standard matplotlib

    This function deletes the pylab default axes and draws visually more
    pleasing axes. In addition it takes care that all graphics remain within
    the area covered by the axes.

    :Parameters:
        *ax* :
            the axes object in which the graphics should stay
        *xtics* :
            an array of x-ticks
        *xfmt* :
            format string for the x-ticks
        *ytics* :
            an array of y-ticks
        *yfmt* :
            format string for the y-ticks
        *xname* :
            label for the x-axis
        *yname* :
            label for the y-axis
    """

    # Data ranges
    yr = ytics.max()-ytics.min()
    xr = xtics.max()-xtics.min()

    # x axis
    yt = ytics.min()
    ax.plot ( [xtics.min(),xtics.max()], [yt-.05*yr]*2, 'k-' )
    for xt in xtics:
        ax.plot ( [xt]*2, [yt-.05*yr,yt-.03*yr], 'k-' )
        ax.text ( xt, yt-.08*yr, xfmt % (xt,), horizontalalignment="center", verticalalignment="top",fontsize=10 )
    ax.text ( xtics.mean(), yt-.14*yr, xname, horizontalalignment="center", verticalalignment="top", fontsize=16 )

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
        *InferenceObject* :
            a BootstrapInference or BayesInference object containing
            the actual inference data
        *ax* :
            an axes object where the plot should go
        *regressor* :
            plot deviance residuals against model prediction (p) or
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

def plotppScatter ( simdata, observed, quantity, shortname=None, ax=None ):
    """plot a scatter diagram to compare observed and predicted data

    :Parameters:
        *simdata* :
            data simulated from the model (typically data from posterior predictives)
        *obseved* :
            observed data (transformed in the same way as the posterior predictives)
        *quantity* :
            name of the quantity that is checked
        *shortname* :
            abbreviation of the plotted quantity
        *ax* :
            pylab axes object where the plot should go.
    """
    if ax==None:
        ax = p.gca()

    ax.plot ( simdata, observed, '.' )
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    axmin = N.min ( list(xl)+list(yl) )
    axmax = N.max ( list(xl)+list(yl) )
    ax.plot ( [axmin,axmax],[axmin,axmax], 'k:' )
    ax.set_xlim ( axmin, axmax )
    ax.set_ylim ( axmin, axmax )

    xt = ax.get_xticks()
    yt = ax.get_yticks()

    drawaxes ( ax, xt, "%g", yt, "%g", "predicted "+quantity, "observed "+quantity )

    # Write diagnostics
    pval = N.mean( (simdata-observed)>=0 )
    ax.text ( xt.min(), yt.max()+.1, "Bayesian p (%s)=%.3f" % (shortname,pval),\
            horizontalalignment="left", verticalalignment="bottom", fontsize=8 )

    if pval<0.975 and pval>0.025:
        return True
    else:
        return False

def plotHistogram ( simdata, observed, xname, shortname=None, ax=None, hideobserved=False, reference="bootstrap" ):
    """plot a histogram and compare observed data to it

    :Parameters:
        *simdata* :
            an array of monte-carlo samples of the parameter of interest
        *observed* :
            observed value of the parameter of interest (for MCMC samples, it is often
            reasonable to use this as the value of 'no effect' or something)
        *xname* :
            name of the paramter of interest
        *shortname* :
            short name of the parameter of interest
        *ax* :
            axes object defining the area where the plot should go
        *hideobserved* :
            if this is True, the observed value is not plotted
        *reference* :
            reference of the data. Could be either a string 'bootstrap'/'bayes' or a number
            against which the histogram is tested

    :Output:
        returns a boolean value indicating whether or not the Null-Hypothesis that
            observed was drawn from the same distribution as simdata is true
    """
    if ax is None:
        ax = p.axes()

    if reference.lower()[:5]==  "boots":
        reference = observed
    elif reference.lower()[:5]=="bayes":
        reference = 0


    # Remove nan
    simdata = N.nan_to_num ( simdata )
    simdata = simdata[simdata!=0]

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

    if reference>p25 and reference<p975:
        return True
    else:
        return False

def plotPMF ( InferenceObject, xlabel_text="Stimulus intensity", ylabel_text=None,ax=None, showaxes=True, showdesc=True, **kwargs ):
    """Show the psychometric function and data in an axes system

    This function plots the best fitting psychometric function and with the
    corresponding data points. If data points are labelled influential, they
    are plotted as red squares, if data points are labelled as outliers, they
    are plotted as red triangles.
    The function uses its internal knowledge about the task (nAFC or Yes/No)
    to put the correct labels to the y-axis.

    :Parameters:
        *ax* :
            axes object in which the plot should go (default: current)
        *xlabel_text* :
            label for the x-axis
        *ylabel_text* :
            label for the y-axis, if this is None, the functions
            determines the correct label from its internal knowledge
            about the task
        *showaxes* :
            if this is set to False, no axes will be drawn
        *showdesc* :
            if this is set to False, no convergence description is drawn
        all other keyword arguments will be passed to the plot function.

    :Example:
    You can use this function to plot multiple psychometric functions. This
    is demonstrated by the example below:

    >>> d0 = [[0, 28, 50], [2, 33, 50], [4, 38, 50], [6, 45, 50], [8, 45, 50], [10, 49, 50]]
    >>> d1 = [[0, 22, 50], [2, 34, 50], [4, 31, 50], [6, 42, 50], [8, 42, 50], [10, 46, 50]]
    >>> d2 = [[0, 26, 50], [2, 31, 50], [4, 38, 50], [6, 47, 50], [8, 49, 50], [10, 49, 50]]
    >>> constraints = ("","","Uniform(0,.1)")
    >>> B0 = BootstrapInference ( d0, priors=constraints )
    >>> B1 = BootstrapInference ( d1, priors=constraints )
    >>> B2 = BootstrapInference ( d2, priors=constraints )
    >>> plotPMF ( B0, showaxes=False, showdesc=False, color='b', outliercolor='b', influentialcolor='b', linewidth=2 )
    >>> plotPMF ( B1, showaxes=False, showdesc=False, color='r', outliercolor='r', influentialcolor='r' )
    >>> plotPMF ( B2, showaxes=True,  showdesc=False, color='k', outliercolor='k', influentialcolor='k', linestyle='--' )

    Note that the last call to plotPMF sets showaxes to True and thus draws the axes.
    """
    if ax==None:
        ax = p.gca()

    # Plot the psychometric function
    xmin = InferenceObject.data[:,0].min()
    xmax = InferenceObject.data[:,0].max()
    x = N.mgrid[xmin:xmax:100j]
    # psi = N.array(_psipy.diagnostics ( x, self.estimate, sigmoid=self.model["sigmoid"], core=self.model["core"], nafc=self.model["nafc"] ))
    psi = InferenceObject.evaluate ( x )
    ax.plot(x,psi,
            color     = kwargs.setdefault ( 'color', 'b' ),
            linestyle = kwargs.setdefault ( 'linestyle', '-' ),
            linewidth = kwargs.setdefault ( 'linewidth', 1 ),
            label     = kwargs.setdefault ( 'label', 'psychometric function' )
            )

    # Plot the data
    xd = InferenceObject.data[:,0]
    pd = InferenceObject.data[:,1].astype("d")/InferenceObject.data[:,2]
    nd = InferenceObject.data[:,2]
    ax.scatter ( xd, pd, s=nd, c=kwargs.setdefault ( 'color', 'b' ),marker=kwargs.setdefault("markertype", "o") )

    # Check axes limits
    ymin,ymax = -.05,1.05
    if InferenceObject.model["nafc"]>1:
        ylabel_text = "P(correct)"
    else:
        if ylabel_text is None:
            ylabel_text = "P(Yes)"

    # Determine tics
    p.setp(ax,frame_on=False,ylim=(ymin,ymax))
    xtics = p.getp(ax,'xticks')
    ytics = list(p.getp(ax,'yticks'))
    # Clean up ytics
    if InferenceObject.model["nafc"]==1:
        for k,yt in enumerate(ytics):
            if yt<0 or yt>1:
                ytics.pop(k)
    else:
        for k,yt in enumerate(ytics):
            if yt<0 or yt>1:
                ytics.pop(k)
    ytics = N.array(ytics)

    if showaxes:
        drawaxes ( ax, xtics, "%g", ytics, "%g", xlabel_text, ylabel_text )

    # Write some model information
    if showdesc:
        if not InferenceObject.deviance is None:
            ax.text(0.5*(xmin+xmax),ymin+.05,"D=%g" % ( InferenceObject.deviance, ) )
        ax.text ( 0.5*(xmin+xmax),ymin+.1,InferenceObject.desc )

    return xtics.min(),xtics.max()

def plotThres ( InferenceObject, ax=None, color="b" ):
    """Plot thresholds and confidence intervals

    :Parameters:
        *InferenceObject* :
            either a BootstrapInference object or a BayesInference object
        *ax* :
            a pylab.axes object to be used for the plot.
        *color* :
            a pylab color marker
    """
    if ax == None:
        ax = p.gca()

    # Determine the range where the data live
    datarange = InferenceObject.data[:,0].min(),InferenceObject.data[:,0].max()
    dataw = datarange[1]-datarange[0]

    for k,cut in enumerate(InferenceObject.cuts):
        c25,c975 = InferenceObject.getCI ( cut=k, conf=(.025,.975) )
        thres = InferenceObject.getThres ( cut )
        ylev = InferenceObject.evaluate ( [thres] )
        if c25 < datarange[0]-dataw*0.2:
            bar = [ datarange[0]-dataw*0.2 ]
            markers = ["<"]
            c25out = True
            ax.text ( datarange[0]-dataw*0.2,ylev,"%g"%(c25,), horizontalalignment="center", fontsize=7 )
        else:
            bar = [ c25 ]
            c25out = False
            markers = ["|"]
        if thres > datarange[0] and thres < datarange[1]:
            bar.append(thres)
            markers.append("|")
        if c975 > datarange[1]+dataw*0.2:
            bar.append(datarange[1]+dataw*0.2)
            c975out = True
            ax.text ( datarange[1]+dataw*0.2,ylev,"%g"%(c975,), horizontalalignment="center", fontsize=7 )
            markers.append(">")
        else:
            bar.append(c975)
            c975out = False
            markers.append("|")
        ax.plot ( bar,[ylev]*len(bar), '-', color=color )
        for x,m in zip(bar,markers):
            ax.plot ( x, ylev, marker=m, color=color )

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
        *warn* :
            if warn is set to True, red warning messages are displayed
            whenever the fit does not seem to describe the data well.
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Goodness of fit diagnostics require monte carlo samples. Try to call the sample() method of your inference object."

    p.figure(figsize=(10,8))

    # First part: Data and fitted function, bottom deviance
    ax = p.axes([0,.5,.33,.5] )
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        InferenceObject.drawposteriorexamples ( ax=ax )
    plotThres ( InferenceObject, ax=ax )
    plotPMF   ( InferenceObject, ax=ax )
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        distname = "posterior"
        observed = -2*N.log(InferenceObject.nullevidence)
    else:
        distname = "bootstrap"
        observed = InferenceObject.deviance
    ax = p.axes ( [0,0,.33,.5] )
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        good = plotppScatter ( InferenceObject.ppdeviance, InferenceObject.mcdeviance, "deviance", "D", ax)
    elif InferenceObject.__repr__().split()[1] == "BootstrapInference":
        good = plotHistogram ( InferenceObject.mcdeviance, observed, distname+" deviance", "D", ax )
    if warn and not good:
        ax.text ( N.array(ax.get_xlim()).mean(), N.array(ax.get_ylim()).mean(),
                "The fitted model is a bad\ndescription of the data!",
                fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    # Second part: Correlations between model prediction and residuals
    ax = p.axes([.33,.5,.33,.5])
    plotRd ( InferenceObject, ax, "p" )
    ax = p.axes([.33,0,.33,.5])
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        good = plotppScatter ( InferenceObject.ppRpd, InferenceObject.mcRpd, "Rpd", "Rpd", ax )
    elif InferenceObject.__repr__().split()[1] == "BootstrapInference":
        good = plotHistogram ( InferenceObject.mcRpd, InferenceObject.Rpd, distname+" Rpd", "Rpd", ax, reference=InferenceObject.__repr__().split()[1] )
    if warn and not good:
        ax.text ( 0, N.mean(p.getp(ax,'ylim')) , "Simulated Rpd differs from observed!\nTry other sigmoid?", \
                fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

    # Third part: Correlations between model prediction and block index
    ax = p.axes([.66,.5,.33,.5])
    plotRd ( InferenceObject, ax, "k" )
    ax = p.axes([.66,0,.33,.5])
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        good = plotppScatter ( InferenceObject.ppRkd, InferenceObject.mcRkd, "Rkd", "Rkd", ax )
    elif InferenceObject.__repr__().split()[1] == "BootstrapInference":
        good = plotHistogram ( InferenceObject.mcRkd, InferenceObject.Rkd, distname+" Rkd", "Rkd", ax, reference=InferenceObject.__repr__().split()[1])
    if warn and not good:
        ax.text ( 0, N.mean(p.getp(ax,'ylim')), "Simulated Rkd differs from observed!\nData are nonstationary!",\
                fontsize=16, color=__warnred, horizontalalignment="center", verticalalignment="center", rotation=45 )

def plotGeweke ( BayesInferenceObject, parameter=0, ax=None, warn=True ):
    """Geweke plot of moving average of samples

    :Parameters:
        *BayesInferenceObject* :
            a BayesInference object that contains all the
            infromation about the sampling process
        *parameter* :
            index of the model parameter of interest
        *ax* :
            the pylab.axes object where the plot should go
        *warn* :
            should a warning message be displayed if non stationarity
            of the samples is observed?
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "Geweke MCMC convergence diagnostic requires monte carlo samples. Try to call the sample() method of your inference object."

    stationary,z,bad = BayesInferenceObject.geweke ( parameter )

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
        nsegments = z.shape[0]
        p.text(0.5*nsegments,0,"chains did not converge" , color=__warnred, fontsize=16, rotation=45, verticalalignment="center", horizontalalignment="center" )

def plotChains ( BayesInferenceObject, parameter=0, ax=None, raw=False, warn=True ):
    """Simply plot all chains for a single parameter

    :Parameters:
        *parameter* :
            index of the model parameter to plot
        *raw* :
            plot raw samples instead of thinned samples after burnin
        *ax* :
            axes in which to print
        *warn* :
            if True, warnings are written into the plot
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "Plotting MCMC chains requires monte carlo samples. Try to call the sample() method of your inference object."

    # Do we have an appropriate axis?
    if ax==None:
        ax = p.axes()

    # Plot the chains
    for c in xrange(BayesInferenceObject.nchains):
        samples = BayesInferenceObject.getsamples ( c, raw=raw )
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
                horizontalalignment="center",verticalalignment="center",fontsize=16,rotation=45,color=__warnred)

    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", "sample #", BayesInferenceObject.parnames[parameter] )

def plotParameterDist ( InferenceObject, parameter=0, ax=None ):
    """Plot the distribution of parameters

    :Parameters:
        *InferenceObject* :
            either a BootstrapInference object or a BayesInference object
            containing the samples of the parameter distribtution
        *parameter* :
            index of the model parameter of interest
        *ax* :
            pylab.axes object where the plot should go
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distribution of parameters requires monte carlo samples. Try to call the sample() method of your inference object."

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

def plotThresholdDist ( InferenceObject, cut=0, ax=None ):
    """Plot the distribution of thresholds

    :Parameters:
        *InferenceObjecxt* :
            a BootstrapInference or BayesInference object containing the desired
            data
        *cut* :
            index (!) of the desired cut
        *ax* :
            axes object to place the plot in.
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of thresholds requires monte carlo samples. Try to call the sample() method of your inference object."

    if ax is None:
        ax = p.axes()

    # Plot histogram
    mcthres = InferenceObject.mcthres[:,cut]
    h,b,ptch = p.hist ( mcthres, bins=20, normed=True, histtype="step", lw=2 )

    # Store ticks
    xtics = N.array(ax.get_xticks())
    ytics = N.array(ax.get_yticks())

    # Highlight estimate and credibility intervals
    thres = InferenceObject.getThres ( InferenceObject.cuts[cut] )
    c25,c975 = InferenceObject.getCI ( cut, (0.025,0.975) )
    ym = ytics.max()
    p.plot( [c25]*2,[0,ym],'b:', [c975]*2,[0,ym],'b:' )
    p.plot ( [thres]*2, [0,ym], 'b' )
    p.text ( xtics.mean(), ym, "F^{-1}(%.2f)=%.3f, CI(95%%)=(%.3f,%.3f)" % (InferenceObject.cuts[cut], thres, c25, c975 ),
            fontsize=8, horizontalalignment="center", verticalalignment="bottom" )

    drawaxes ( ax, xtics, "%g", ytics, "%g", "F^{-1}(%.2f)" % (InferenceObject.cuts[cut],), "density estimate" )

def ThresholdPlot ( InferenceObject ):
    """Show distributions and estimates for all thresholds

    :Parameters:
        *InferenceObject*
            a BootstrapInference or BayesInference object containing the
            desired data
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of thresholds requires monte carlo samples. Try to call the sample() method of your inference object."

    nthres = len(InferenceObject.cuts)
    axw = 1./nthres
    fig = p.figure ( figsize=(3*nthres,3) )

    for k in xrange ( nthres ):
        ax = p.axes ( [axw*k,0,axw,1] )
        plotThresholdDist ( InferenceObject, k, ax )

def ParameterPlot ( InferenceObject ):
    """Show distributions and estimates for all parameters in the model

    :Parameters:
        *InferenceObject* :
            a BootstrapInference or BayesInference object containing the
            desired data
    """

    if InferenceObject.mcestimates is None:
        raise ValueError, "Plotting distributions of parameters requires monte carlo samples. Try to call the sample() method of your inference object."

    nparams = len(InferenceObject.parnames)
    axw = 1./nparams
    fig = p.figure (figsize=(3*nparams,3))

    for k in xrange ( nparams ):
        ax = p.axes ( [axw*k,0,axw,1] )
        plotParameterDist ( InferenceObject, k, ax )

def ConvergenceMCMC ( BayesInferenceObject, parameter=0, warn=True ):
    """Diagram to check convergence of MCMC chains for a single parameter

    :Parameters:
        *BayesInferenceObject* :
            a BayesInference object containing all information about
            the model and the posterior distribution
        *parameter* :
            model parameter of interest. So far, no model derived parameters such as
            thresholds are supported
        *warn* :
            should warnings be displayed if the samples look suspicious?
    """

    if BayesInferenceObject.mcestimates is None:
        raise ValueError, "MCMC convergence diagnostics require monte carlo samples. Try to call the sample() method of your inference object."

    fig = p.figure ( figsize=[9,3] )
    ax =  p.axes ( [0,0.,0.33,1] )
    plotChains ( BayesInferenceObject, parameter, ax, warn=warn )
    ax = p.axes ( [.33,0,.33,1] )
    plotGeweke ( BayesInferenceObject, parameter, ax, warn=warn )
    ax = p.axes ( [.66,0,.33,1] )
    plotParameterDist ( BayesInferenceObject, parameter, ax )

def plotSensitivity ( BootstrapInferenceObject, ax=None ):
    """Visualize a sensitivity analysis to determine expanded bootstrap confidence intervals

    Sensitivity analysis is used for BootstrapInference objects to expand the confidence intervals
    in order to obtain more realistic coverage. This function calls the sensitivity_analysis() method
    of the BootstrapInferenceObject with default parameters. If other parameters are requested, the
    sensitivity_analysis() method should be called manually

    :Parameters:
        *BootstrapInferenceObject* :
            Inference object to be analyzed
        *ax* :
            pylab axes that should be used for plotting
    """
    if BootstrapInferenceObject.mcestimates is None:
        raise ValueError, "Sensitivity analysis requires monte carlo samples. Try to call the sample() method of your inference object."

    if ax==None:
        ax = p.axes()

    # Determine axes ranges
    prm1 = BootstrapInferenceObject.mcestimates[:,0]
    prm2 = BootstrapInferenceObject.mcestimates[:,1]
    ax.plot(prm1,prm2,'w.',markersize=1)
    xmin,xmax = ax.get_xlim()
    ymin,ymax = ax.get_ylim()
    ax.cla()

    # Plot the density estimate in the background
    x,y = N.mgrid[xmin:xmax:100j,ymin:ymax:100j]
    C = BootstrapInferenceObject.mcdensity(N.c_[N.ravel(x),N.ravel(y)].T)
    C.shape = 100,100
    ax.imshow( C.T,origin="lower",extent=(xmin,xmax,ymin,ymax), cmap=p.cm.gray_r )

    # Get the points and make sure, a sensitivity_analysis has indeed been run
    dummy,points = BootstrapInferenceObject.sensitivity_analysis(verbose=False)

    # plot the points
    ax.fill(points[:,0],points[:,1],fill=False,edgecolor="r",linewidth=2)
    ax.plot(prm1,prm2,"b.",markersize=2)
    ax.plot(points[:,0],points[:,1],'rd',markersize=5)
    ax.plot([BootstrapInferenceObject.estimate[0]],[BootstrapInferenceObject.estimate[1]],'ro',markersize=5)

    # plot marginal percentiles
    prm1lims = p.prctile ( BootstrapInferenceObject.mcestimates[:,0], (2.5,25,75,97.5) )
    prm2lims = p.prctile ( BootstrapInferenceObject.mcestimates[:,1], (2.5,25,75,97.5) )
    ax.plot( prm1lims, [ymin-0.05*(ymax-ymin)]*4, 'b-', [xmin-0.05*(xmax-xmin)]*4, prm2lims, 'b-' )
    ax.plot( prm1lims[1:3], [ymin-0.05*(ymax-ymin)]*2, 'b-', [xmin-0.05*(xmax-xmin)]*2, prm2lims[1:3], 'b-', linewidth=5 )

    # Draw axes
    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", BootstrapInferenceObject.parnames[0], BootstrapInferenceObject.parnames[1] )

def plotInfluential ( InferenceObject ):
    """Diagnostic plot for detecting influential observations

    Determining influential observations follows a different logic for bootstrap
    and for bayes inference. A block is labelled an influential observation if
    the fit for a dataset without that point is significantly different from the
    fit including that point. For BootstrapInference objects, this is quantified
    using a normed distance of the maximum likelihood fit including the block and
    withouth that block. This distance is normed in the following way: If the
    maximum likelihood fit for the reduced dataset remains inside the 95% confidence
    limits of the maximum likelihood fit for the full dataset, the influence
    value is below 1. Thus, influence values large than 1 indicate a problem with
    the data set. For BayesInference objects, the influence of a block is simply
    quantified as the Kullbach-Leibler divergence of the posterior for the full
    data set from the posterior for the reduced data set.

    :Parameters:
        *InferenceObject* :
            Data set for which the influential observations are to be plotted
    """
    maxinfl = N.argmax(InferenceObject.infl)
    ind = range ( InferenceObject.data.shape[0] )
    ind.pop(maxinfl)
    influencedDataset = psignidata.BootstrapInference( InferenceObject.data[ind,:], priors=InferenceObject.model["priors"],
            sample=False)

    p.figure ( figsize=(6,8) )
    ax = p.axes ( (0.0,.5,.9,.5) )
    if InferenceObject.__repr__().split()[1] == "BayesInference":
        InferenceObject.drawposteriorexamples ( ax=ax )
    plotPMF ( InferenceObject, ax=ax, showaxes=False, showdesc=False, color="b", linewidth=2 )
    ax.plot ( [InferenceObject.data[maxinfl,0]], [InferenceObject.data[maxinfl,1].astype("d")/InferenceObject.data[maxinfl,2]],
            'rx', markersize=20, markeredgewidth=5 )
    xl = plotPMF ( influencedDataset, ax=ax, showdesc=False, showaxes=True, color="r", markertype=([(0,0)],0), linewidth=2 )

    ax = p.axes ( (0.0, 0., .9, .5) )
    if InferenceObject.__repr__().split()[1] == "BootstrapInference":
        ax.plot ( [InferenceObject.data[:,0].min(),InferenceObject.data[:,0].max()], [1,1], 'k:' )
        yname = "Influence"
    else:
        yname = "D_KL( full || reduced )"
    ax.plot ( InferenceObject.data[:,0], InferenceObject.infl, 'bo' )
    ax.set_xlim(xl)
    drawaxes ( ax, ax.get_xticks(), "%g", ax.get_yticks(), "%g", "stimulus intensity", yname )


gof = GoodnessOfFit
