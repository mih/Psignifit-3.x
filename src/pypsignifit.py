#!/usr/bin/env python

import numpy as N
import pylab as p
import _psipy

class PsiInference ( object ):
    def __init__ ( self ):
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
        self.outl              = None
        self.infl              = None

    def pmfanddata ( self, ax=None, xlabel_text="Stimulus intensity" ):
        """Show the psychometric function and data in an axes system"""
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
            ylabel_text = "P(correct)"
        else:
            ymin,ymax = -.05,1.05
            ylabel_text = "P(Yes)"
        lenx = xmax-xmin
        p.setp(ax,frame_on=False,ylim=(ymin,ymax))
        xtics = p.getp(ax,'xticks')
        ytics = p.getp(ax,'yticks')

        # xaxis
        ax.plot([xmin,xmax],[ymin+.01]*2,'k-',[xmin-.04*lenx]*2,[ymin+.05,1],'k-',linewidth=1)
        for xt in xtics:
            ax.plot([xt]*2,[ymin+.01,ymin+.02],'k-')
            ax.text(xt,ymin,"%2.1f" % (xt,), verticalalignment="top",horizontalalignment="center")
        ax.text(0.5*(xmin+xmax),ymin-.1,xlabel_text,fontsize=16,horizontalalignment="center")

        # yaxis
        for yt in ytics:
            if (yt>=ymin and yt <= 1):
                ax.plot([xmin-.04*lenx,xmin-.03*lenx],[yt]*2,'k-')
                ax.text(xmin-.05*lenx,yt,"%2.2f" % (yt,), verticalalignment="center",horizontalalignment="right")
        ax.text(xmin-.28*lenx,0.5*(ymin+ymax),ylabel_text,fontsize=16,verticalalignment="center",rotation=90)

        if not self.deviance is None:
            ax.text(0.5*(xmin+xmax),ymin+.05,"D=%g" % ( self.deviance, ) )
        ax.text ( 0.5*(xmin+xmax),ymin+.1,self.desc )

        p.setp(ax,xlim=(xmin-0.35*lenx,xmax+0.05*lenx),ylim=(ymin-.2,ymax),xticks=(),yticks=())

    def plotRd ( self, ax=None, regressor="p" ):
        """plot deviance residuals agains model prediction"""
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
            for k,xt in enumerate(xtics):
                if xtics[k]>1.:
                    xtics.pop(k)
        xtics = N.array(xtics)
        ytics = p.getp(ax,"yticks")

        if regressor=="p":
            ax.text(psilims.mean(),ytics[-2],"Rpd=%.3f" % ( self.Rpd, ) )
            xname = "model prediction"
        elif regressor=="k":
            ax.text(psilims.mean(),ytics[-2],"Rkd=%.3f" % ( self.Rkd, ) )
            xname = "block index"

        # x-axis
        yt = ytics.min()
        yr = ytics.max()-yt
        ax.plot([xtics[0],xtics[-1]],[yt-0.04*yr]*2,'k-')
        for xt in xtics:
            ax.plot([xt]*2,[yt-0.04*yr,yt-0.03*yr],'k-')
            ax.text(xt,yt-.05*yr,"%g" % (xt,), verticalalignment="top",horizontalalignment="center")
        ax.text(xtics.mean(),yt-.25*yr,xname,fontsize=16,horizontalalignment="center")

        # y-axis
        xt = xtics.min()
        xr = xtics.max()-xt
        ax.plot([xt-0.04*xr]*2,[ytics[0],ytics[-1]],'k-')
        for yt in ytics:
            ax.plot([xt-0.04*xr,xt-0.03*xr],[yt]*2,'k-')
            ax.text(xt-.05*xr,yt,"%.2f" % (yt,), verticalalignment="center",horizontalalignment="right")
        ax.text(xt-.3*xr,ytics.mean(),"deviance residuals",verticalalignment="center",fontsize=16,rotation=90)

        # Delete the uninteresting stuff
        p.setp(ax,frame_on=False,xticks=(),yticks=())
        p.setp(ax,xlim=(xtics.min()-0.4*xr,xtics.max()+0.04*xr),ylim=(ytics.min()-0.4*yr,ytics.max()+0.1*yr))

    desc = property ( fget=lambda self: "sigmoid: %(sigmoid)s\ncore: %(core)s\nnAFC: %(nafc)d" % self.model )

class BootstrapInference ( PsiInference ):
    def __init__ ( self, data, sample=False, **kwargs ):
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
        """
        self.data = data
        self.model = {
                "sigmoid": kwargs.setdefault("sigmoid","logistic"),
                "core":    kwargs.setdefault("core",   "ab"),
                "priors":  kwargs.setdefault("priors", None),
                "nafc":    kwargs.setdefault("nafc",    2)
                }
        self.estimate,self.deviance = _psipy.mapestimate(self.data,**self.model)
        self.predicted,self.devianceresiduals,self.deviance,self.Rpd,self.Rkd = _psipy.diagnostics(self.data,self.estimate)

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
                self.__bRkd,self.__bRpd,self.__outl,self.__infl = _psipy.bootstrap(self.data,self.estimate,Nsamples,**self.model)
        self.__outl = N.array(self.__outl,dtype=bool)
        self.__infl = N.array(self.__infl,dtype=bool)

    def plothistogram ( self, simdata, observed, xname, ax=None ):
        """plot a histogram and compare observed data to it"""
        if ax is None:
            ax = p.axes()

        if xname[0] == "R":
            ax.hist ( simdata, bins=N.arange(-1,1,.1) )
        else:
            ax.hist ( simdata, bins=20 )

        xtics = p.getp(ax,"xticks")
        ytics = p.getp(ax,"yticks")
        p25,p975 = p.prctile ( simdata, (2.5,97.5) )

        yr = ytics.max()-ytics.min()
        xr = xtics.max()-xtics.min()
        yy = [ytics.min(),ytics.max()+0.02*xr]
        ax.plot ( [observed]*2, yy, 'r', linewidth=2 )
        ax.plot ( [p25]*2, yy, 'r:', [p975]*2, yy, 'r:' )

        # y axes
        xt = xtics.min()
        ax.plot ( [xt-.02*xr]*2, [ytics.min(),ytics.max()], 'k-' )
        for yt in ytics:
            ax.plot ( [xt-.02*xr,xt-.01*xr], [yt]*2, 'k-' )
            ax.text ( xt-.03*xr,yt, "%d" % (yt,), verticalalignment="center", horizontalalignment="right" )
        ax.text ( xt-.25*xr, ytics.mean(), "number per bin", fontsize=16, verticalalignment="center", rotation=90 )

        # x axes
        yt = xtics.min()
        ax.plot ( [xtics.min(),xtics.max()], [yt-.04*yr]*2, 'k-' )
        for xt in xtics:
            ax.plot ( [xt]*2, [yt-.04*yr,yt-.03*yr], 'k-' )
            ax.text ( xt, yt-.05*yr, "%.1f" % (xt,), verticalalignment="top", horizontalalignment="center" )
        ax.text ( xtics.mean(), yt-.2*yr, xname, fontsize=16, horizontalalignment="center" )

        # Write diagnostics
        yt = ytics.max()
        ax.text ( xtics.min(), yt+.1, "%s=%.3f, c(2.5%%)=%.3f, c(97.5%%)=%.3f" % (xname,observed,p25,p975), horizontalalignment="left",verticalalignment="bottom", fontsize=8 )

        p.setp ( ax, frame_on=False, xticks=(), yticks=() )
        p.setp ( ax, xlim=(xtics.min()-.3*xr,xtics.max()+.01*xr), ylim=(ytics.min()-.3*yr,ytics.max()+.01*yr) )

    def diagnostics ( self ):
        """Make a diagnostic plot"""
        p.figure(figsize=(10,8))
        self.pmfanddata ( p.axes([0,.5,.33,.5] ) )
        self.plothistogram ( self.bdeviance, self.deviance, "deviance", p.axes([0,0,.33,.5]) )
        self.plotRd ( p.axes([.33,.5,.33,.5]), "p" )
        self.plothistogram ( self.bRpd, self.Rpd, "Rpd", p.axes([.33,0,.33,.5]) )
        self.plotRd ( p.axes([.66,.5,.33,.5]), "k" )
        self.plothistogram ( self.bRkd, self.Rkd, "Rkd", p.axes([.66,0,.33,.5]) )

    outl = property ( fget=lambda self: self.__outl )
    infl = property ( fget=lambda self: self.__infl )
    bdeviance = property ( fget=lambda self: self.__bdeviance )
    bRpd = property ( fget=lambda self: self.__bRpd )
    bRkd = property ( fget=lambda self: self.__bRkd )

def main ( ):
    x = [3]+[float(2*k) for k in xrange(6)]
    k = [25,34,32,40,48,50,48]
    n = [50]*7
    d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
    d = N.array(zip(x,k,n))
    priors = ("flat","flat","Uniform(0,0.1)")
    b = BootstrapInference ( d, sample=2000, priors=priors )

    # b.pmfanddata()
    # b.plotRd (regressor="k")
    # b.plothistogram ( b.bdeviance, b.deviance, "deviance" )
    b.diagnostics()

    p.show()

if __name__ == "__main__":
    main()
