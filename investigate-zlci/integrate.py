#!/usr/bin/env python

import numpy as np
from scipy.optimize import fmin,leastsq
import pypsignifit as pf
import pypsignifit.psignipriors as pfp
import swignifit.swignifit_raw as sfr
import swignifit.utility as sfu
import pypsignifit.psigobservers as pfo

def integration_grid ( data, gridsize=30 ):
    data = np.array(data)
    mprior,mmin,mmax = pfp.default_mid ( data[:,0] )
    wprior,wmin,wmax = pfp.default_width ( data[:,0] )
    lprior,lmin,lmax = pfp.default_lapse ( )
    gprior,gmin,gmax = pfp.default_lapse ( )
    priors = (mprior,wprior,lprior,gprior)

    print data[:,0]
    print "m",mmin,mmax
    print "w",wmin,wmax
    print "l",lmin,lmax
    print "g",gmin,gmax

    # mapest = pf.BootstrapInference ( data, priors ).estimate

    print gridsize*1j
    grid = np.reshape ( np.mgrid[mmin:mmax:1j*gridsize,wmin:wmax:1j*gridsize,lmin:lmax:1j*gridsize,gmin:gmax:1j*gridsize], (4,-1) )

    pdata,ppmf,pn = sfu.make_dataset_and_pmf ( data, 1, "logistic", "mw0.1", priors )

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

    mu = x[0][np.argmax(fx[0])]
    fx[0] /= fx[0].max()
    mprm = fmin ( error_gauss, [1.,mu,1.5], args=(fx[0],x[0]) )
    print mprm
    post.append ( "Gauss(%g,%g)" % ( mprm[1],mprm[2]**2 ) )

    fx[1] /= fx[1].max()
    wprm = fmin ( error_gamma, [1.,2,4], args=(fx[1],x[1]) )
    post.append ( "Gamma(%g,%g)" % ( wprm[1]**2,wprm[2]**2 ) )

    fx[2] /= fx[2].max()
    lprm = fmin ( error_beta, [1.,2,20], args=(fx[2],x[2]) )
    post.append ( "Beta(%g,%g)" % ( lprm[1],lprm[2] ) )

    if len(fx)>3:
        fx[3] /=  fx[3].max()
        gprm = fmin ( error_beta, [1.,2,20], args=(fx[3],x[3]) )
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

    showfits = True

    pl.subplot(221)
    pl.plot ( x[0], fx[0], 'o' )
    if showfits:
        xx = np.mgrid[x[0].min():x[0].max():100j]
        pl.plot ( xx, [f[0].pdf(x_) for x_ in xx], '-' )
    pl.title ( "m"  )

    pl.subplot(222)
    pl.plot ( x[1], fx[1], 'o')
    if showfits:
        xx = np.mgrid[x[1].min():x[1].max():100j]
        pl.plot ( xx, [f[1].pdf(x_) for x_ in xx], '-' )
    pl.title ( "w" )

    pl.subplot(223)
    pl.plot ( x[2], fx[2], 'o' )
    if showfits:
        xx = np.mgrid[x[2].min():x[2].max():100j]
        pl.plot ( xx, [f[2].pdf(x_) for x_ in xx], '-' )
    pl.title ( "lm" )

    pl.subplot(224)
    pl.plot ( x[3], fx[3], 'o')
    if showfits:
        xx = np.mgrid[x[3].min():x[3].max():100j]
        pl.plot ( xx, [f[3].pdf(x_) for x_ in xx], '-' )
    pl.title ( "gm" )
    pl.show()
