#/usr/bin/env python

import pylab as p

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
