import numpy as np
import swignifit as sf

def bootstrap(data, start=None, nsamples=2000, nafc=2, sigmoid="logistic",
        core="ab", priors=None, cuts=None, parametric=True ):

    data = np.array(data).T
    x = sf.vector_double(data[0])
    k = sf.vector_int(data[1].astype(int))
    N = sf.vector_int(data[2].astype(int))
    data = sf.PsiData(x,N,k,nafc)
    sigmoid = sf.getsigmoid(sigmoid)
    core = sf.getcore(core, sigmoid.getcode(), data)
    pmf = sf.PsiPsychometric(nafc, core, sigmoid)
    Nparams = pmf.getNparams()
    # THIS IS JUST A PLACEHOLDER
    # REMOVE AS SOON AS THERE IS A GOOD WAY TO DETRMINE CUTS
    cuts = sf.vector_double([1.0, 0.5])
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



