#!/usr/bin/env python

import pypsignifit
import pylab as p

def check_gof ( InferenceObject ):
    """Check whether the inference object fails in the goodness of fit test"""
    # if InferenceObject.inference == "CML-MC":
    if isinstance ( InferenceObject, pypsignifit.BootstrapInference ):
        # Model checking based on deviance
        return InferenceObject.deviance, p.prctile ( InferenceObject.mcdeviance, 97.5 )
    elif InferenceObject.inference == "MCMC":
        # Model checking based on Bayes Factor
        return InferenceObject.evidence, InferenceObject.nullevidence
    else:
        raise AttributeError, "Inference object has invalid attribute"

def check_trend ( InferenceObject ):
    """Check whether the data contain evidence for learning"""
    ci = p.prctile ( InferenceObject.mcRkd, (2.5, 97.5) )
    if InferenceObject.inference == "CML-MC":
        # A significant trend is if the observed correlation is *not* included
        if ci[0]>InferenceObject.Rkd or ci[1]<InferenceObject.Rkd:
            ok = False
        else:
            ok = True
    elif InferenceObject.inference == "MCMC":
        # A significant trend is if 0 correlation is *not* included
        if ci[0]>0 or ci[1]<0:
            ok = False
        else:
            ok = True
    else:
        raise AttributeError, "Inference object has invalid attribute"
    return ci[0],InferenceObject.Rkd,ci[1],ok

def check_descriptive ( InferenceObject ):
    """Check whether the data is well described by the fit"""
    ci = p.prctile ( InferenceObject.mcRpd, (2.5,97.5) )
    if InferenceObject.inference == "CML-MC":
        if ci[0]>InferenceObject.Rpd or ci[1]<InferenceObject.Rpd:
            ok = False
        else:
            ok = True
    elif InferenceObject.inference == "MCMC":
        if ci[0]>0 or ci[1]<0:
            ok = False
        else:
            ok = True
    else:
        raise AttributeError, "Inference object has invalid attribute"
    return ci[0],InferenceObject.Rpd,ci[1],ok

def check_influential ( InferenceObject ):
    """Check for influential observations"""
    if p.any(p.array(InferenceObject.infl)):
        return False
    else:
        return True

def check_outliers ( InferenceObject ):
    """Check for outliers"""
    if p.any(p.array(InferenceObject.outl)):
        return False
    else:
        return True

def check_ci ( InferenceObject ):
    """Determine confidence intervals for threshold"""
    if InferenceObject.inference == "CML-MC":
        ci0,ci1 = InferenceObject.getCI(1)
        InferenceObject.sensitivity_analysis (verbose=False)
        ce0,ce1 = InferenceObject.getCI(1)
        return InferenceObject.thres[1],ci0,ci1,ce0,ce1

def resultsline ( InferenceObject=None ):
    """Write a line with the results"""
    if InferenceObject is None:
        return "D D975 Rkd025 Rkd Rkd975 notrend Rpd025 Rpd Rpd975 nodescr infl outl thres ci025 ci975 ce025 ce975\n"
    elif InferenceObject.inference == "CML-MC":
        out = ""
        out += "%g\t%g"                % check_gof ( InferenceObject )
        out += "\t%g\t%g\t%g\t%d"      % check_trend ( InferenceObject )
        out += "\t%g\t%g\t%g\t%d"      % check_descriptive ( InferenceObject )
        out += "\t%d\t%d"              % ( check_influential ( InferenceObject ), check_outliers ( InferenceObject ) )
        out += "\t%g\t%g\t%g\t%g\t%g"  % check_ci ( InferenceObject )
        out += "\n"
        return out

if __name__ == "__main__":
    import pypsignifit
    import pypsignifit.psigobservers as observers
    d = observers.Observer ( 4,2,.02 ).DoAnExperiment ( [7,2,3,1,5,7,8] )
    B = pypsignifit.BootstrapInference ( d, priors=("","","Uniform(0,.1)"), sample=True )
    print resultsline()
    print resultsline(B)

    pypsignifit.GoodnessOfFit(B)
    p.show()
