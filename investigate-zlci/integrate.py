#!/usr/bin/env python

import numpy as np
from scipy.optimize import fmin,leastsq, brentq
import pypsignifit as pf
import pypsignifit.psignipriors as pfp
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import pypsignifit.psigobservers as pfo

def bounds ( mapest, pdata, ppmf, parameter="m" ):
    prm = mapest.copy()

    if mapest[2] < 0:
        prm[2] = 1e-5
    if len(mapest)>3:
        if mapest[3] < 0:
            prm[3] = 1e-5

    maxpost = -ppmf.neglpost ( prm, pdata )

    x = np.array ( [ pdata.getIntensity ( i ) for i in xrange ( pdata.getNblocks() ) ] )
    if parameter=="m":
        index = 0
        pmin = x.min ()
        pmax = x.max ()
    elif parameter=="w":
        index = 1
        x.sort ()
        pmin = 0.3*np.min ( np.diff(x) )
        pmax = 3*(x.max()-x.min())
    elif parameter=="lm":
        index = 2
        pmin = 0
        pmax = 1
    elif parameter=="gm":
        index = 3
        pmin = 0
        pmax = 1

    def error ( x, target=0.1*maxpost ):
        prm[index] = x
        lpost = -ppmf.neglpost ( prm, pdata )
        print "%s: %.3f %.3f %.3f %.3f %.3f" % ( parameter,lpost,maxpost,lpost-maxpost,target,x )
        return lpost - target

    lmax = pmax
    lmin = pmin
    print "Maxpost =",maxpost
    for r in [.1,.3,.5,.8,.9]:
        trg = maxpost+np.log(r)
        if error ( mapest[index], trg )*error ( pmax, trg ) < 0:
            lmax = brentq ( error, mapest[index], pmax, args=(trg,) )
            print parameter,r,"max"
            break

    for r in [.1,.3,.5,.8,.9]:
        trg = maxpost+np.log(r)
        if error ( mapest[index], trg ) * error ( pmin, trg ) < 0:
            lmin = brentq ( error, pmin, mapest[index], args=(trg,) )
            print parameter,r,"min"
            break

    if lmin==pmin:
        print "WARNING: did not optimize lower bound for parameter",parameter
    if lmax==pmax:
        print "WARNING: did not optimize upper bound for parameter",parameter


    return lmin,lmax

def integration_grid ( data, gridsize=9 ):
    data = np.array(data)
    mprior,mmin,mmax = pfp.default_mid ( data[:,0] )
    wprior,wmin,wmax = pfp.default_width ( data[:,0] )
    lprior,lmin,lmax = pfp.default_lapse ( )
    gprior,gmin,gmax = pfp.default_lapse ( )
    priors = (mprior,wprior,lprior,gprior)

    pdata,ppmf,pn = sfu.make_dataset_and_pmf ( data, 1, "logistic", "mw0.1", priors )
    mapest = pf.BootstrapInference ( data, priors, nafc=1 ).estimate
    mmin,mmax = bounds ( mapest, pdata, ppmf, "m" )
    wmin,wmax = bounds ( mapest, pdata, ppmf, "w" )
    lmin,lmax = bounds ( mapest, pdata, ppmf, "lm" )
    gmin,gmax = bounds ( mapest, pdata, ppmf, "gm" )

    print data[:,0]
    print "m",mmin,mmax
    print "w",wmin,wmax
    print "l",lmin,lmax
    print "g",gmin,gmax


    print gridsize*1j
    grid = np.reshape ( np.mgrid[mmin:mmax:1j*gridsize,wmin:wmax:1j*gridsize,lmin:lmax:1j*gridsize,gmin:gmax:1j*gridsize], (4,-1) )


    post = np.reshape ( np.array ( map ( lambda prm: ppmf.neglpost ( prm, pdata ), grid.T ) ), [gridsize]*pn ) # negative log posterior
    post = np.exp ( -post ) # posterior
    grid = np.reshape ( grid, [pn]+[gridsize]*pn )

    d = [ np.diff ( grid[i], axis=i )[0].max() for i in xrange ( pn ) ]

    fx = [ marginalize ( post, d, i ) for i in xrange ( pn ) ]
    x = []
    for i in xrange ( pn ):
        s =  [ i ] + [0]*pn
        s[i+1] = slice(0,grid.shape[i+1])
        x.append ( grid[s] )
    return x,fx,priors

def marginalize ( post, d, i ):
    p = reduce ( lambda x,y: x*y, d )
    fx = np.zeros ( post.shape[i] )
    for j in xrange ( post.shape[i] ):
        s = [slice(0,post.shape[0]),slice(0,post.shape[1]),slice(0,post.shape[2]),slice(0,post.shape[3])]
        s[i] = j
        fx[j] = post[s].sum()
    # return fx * p / d[i]
    return fx

def error_gauss ( prm, fx, x ):
    Z,mu,sg = prm
    sg = sg*sg
    return np.sum ( ( Z**2*np.exp ( -0.5*((x-mu)/sg)**2 ) - fx )**2 )
    # return Z*Z*np.exp ( -0.5*((x-mu)/sg**2 ) ) - fx

def error_gamma ( prm, fx, x ):
    Z,k,th = prm
    k = k*k
    th = th*th
    return np.sum ( ( Z**2*x**(k-1)*np.exp(-x/th) - fx )**2 )

def error_beta ( prm, fx, x ):
    Z,al,bt = prm
    return np.sum ( ( Z**2*x**(al-1)*(1-x)**(bt-1) - fx )**2 )

def fit_posterior ( fx, x ):
    post = []
    I = 10000
    N = 10000

    mu = x[0][np.argmax(fx[0])]
    fx[0] /= fx[0].max()
    mprm = fmin ( error_gauss, [1.,mu,1.5], args=(fx[0],x[0]), maxfun=N, maxiter=I )
    print mprm
    post.append ( "Gauss(%g,%g)" % ( mprm[1],mprm[2]**2 ) )

    fx[1] /= fx[1].max()
    wprm = fmin ( error_gamma, [1.,2,4], args=(fx[1],x[1]), maxfun=N, maxiter=I )
    post.append ( "Gamma(%g,%g)" % ( wprm[1]**2,wprm[2]**2 ) )

    fx[2] /= fx[2].max()
    lprm = fmin ( error_beta, [1.,2,20], args=(fx[2],x[2]), maxfun=N, maxiter=I )
    post.append ( "Beta(%g,%g)" % ( lprm[1],lprm[2] ) )

    if len(fx)>3:
        fx[3] /=  fx[3].max()
        gprm = fmin ( error_beta, [1.,2,20], args=(fx[3],x[3]), maxfun=N, maxiter=I )
        post.append ( "Beta(%g,%g)" % ( gprm[1],gprm[2] ) )

    return post

if __name__ == "__main__":
    import pylab as pl
    O = pfo.Observer ( 5,3,.05,.05, core="mw0.1", sigmoid="logistic", nafc=1 )
    data = O.DoAnExperiment ( [1,2,3,4,5,6,7,8,12], 30 )
    # data = [[1, 2, 30], [2, 2, 30], [3, 2, 30], [4, 5, 30], [5, 16, 30], [6, 22, 30], [7, 26, 30], [8, 27, 30], [12, 29, 30]]
    x,fx,priors = integration_grid ( data )
    print data

    post = fit_posterior(fx,x)
    print post
    f = [ sfu.get_prior ( p ) for p in post ]

    mapest = pf.BootstrapInference ( data, priors, core="mw0.1", nafc=1 ).estimate

    rng = [(0,10),(0,5),(0,.5),(0,.5)]

    for i,prm in enumerate ( ["m","w","lm","gm"] ):
        pl.subplot(221+i)
        xx = np.mgrid[rng[i][0]:rng[i][1]:1000j]
        g = np.array ( [f[i].pdf(x_) for x_ in xx] )
        pl.plot ( xx, g, '-' )
        mxind = np.argmax(fx[i])
        r = fx[i][mxind]/f[i].pdf(x[i][mxind])
        pl.plot ( x[i], fx[i]/r, 'o' )
        pl.plot ( [O.params[i]]*2,[0,f[i].pdf(O.params[i])], 'k' )
        pl.plot ( [mapest[i]]*2,[0,f[i].pdf(mapest[i])], 'r' )
        pl.title ( prm  )

    print mapest

    pl.show()
