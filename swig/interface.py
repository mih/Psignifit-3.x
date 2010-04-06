import numpy as np
import swignifit as sf
import operator as op

sig_dict = dict()

for subclass in sf.PsiSigmoid.__subclasses__():
    sig_dict[subclass.getDescriptor()] = subclass

class PsignifitException(Exception):
    pass

def get_sigmoid(descriptor):
    if not sig_dict.has_key(descriptor):
        raise PsignifitException("The sigmoid \'"+str(descriptor)+"\' you requested, is not available.")
    return sig_dict[descriptor]()

def get_cuts(cuts):
    if cuts is None:
        return sf.vector_double([0.5])
    elif op.isNumberType(cuts):
        return sf.vector_double([cuts])
    elif op.isSequenceType(cuts) and np.array([op.isNumberType(a) for a in cuts]).all():
        return sf.vector_double(cuts)
    else:
        raise PsignifitException("'cuts' must be either None, a number or a "+\
                "sequence of numbers.")

def available_sigmoids():
    print "The following sigmoids are available:"
    print sig_dict.keys()


def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True ):

    data = np.array(data).T
    x = sf.vector_double(data[0])
    k = sf.vector_int(data[1].astype(int))
    N = sf.vector_int(data[2].astype(int))
    data = sf.PsiData(x,N,k,nafc)
    sigmoid = get_sigmoid(sigmoid)
    core = sf.getcore(core, sigmoid.getcode(), data)
    pmf = sf.PsiPsychometric(nafc, core, sigmoid)
    Nparams = pmf.getNparams()
    cuts = get_cuts(cuts)
    # here we also need to set the priors
    # but again, there is no 'clean' way to do this at the moment
    # REMEMBER TO SOMEHOW SET PRIORS
    bs_list = sf.bootstrap(nsamples, data, pmf, cuts, start, True, parametric)
    jk_list = sf.jackknifedata(data, pmf)



x = [float(2*k) for k in xrange(6)]
k = [34,32,40,48,50,48]
n = [50]*6
d = [[xx,kk,nn] for xx,kk,nn in zip(x,k,n)]
priors = ('flat','flat','Uniform(0,0.1)')
bootstrap(d,nsamples=2000,priors=priors)
#samples,est,D,thres,bias,acc,Rkd,Rpd,out,influ = bootstrap(d,nsamples=2000,priors=priors)



