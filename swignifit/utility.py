#/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################
import operator as op
import re
import numpy as np
import swignifit_raw as sfr

class PsignifitException(Exception):
    pass

def extract_subclasses(base, sub_func):
    to_visit = base.__subclasses__()
    subclasses = dict()
    for cl in to_visit:
        descriptor = eval("cl."+sub_func)
        if descriptor not in subclasses.keys():
            subclasses[descriptor] = cl
            to_visit.extend(cl.__subclasses__())
    return subclasses

def extract_subclasses_descriptor(base):
    return extract_subclasses(base, "getDescriptor()")

def extract_subclasses_reflection(base):
    return extract_subclasses(base, "__name__")

sig_dict = extract_subclasses_descriptor(sfr.PsiSigmoid)
core_dict = extract_subclasses_descriptor(sfr.PsiCore)
prior_dict = extract_subclasses_reflection(sfr.PsiPrior)
sampler_dict = extract_subclasses_reflection(sfr.PsiSampler)

def available_cores():
    print "The following cores are availabe:"
    print core_dict.keys()

def available_sigmoids():
    print "The following sigmoids are available:"
    print sig_dict.keys()

def available_priors():
    print "The following priors are available:"
    print prior_dict.keys()

def available_samplers():
    print "The following mcmc samplers are available:"
    print sampler_dict.keys()

def make_dataset(data, nafc):
    """ create a PsiData object from column based input """
    data = np.array(data).T
    x = sfr.vector_double(data[0])
    k = sfr.vector_int([int(i) for i in data[1].astype(int)])
    N = sfr.vector_int([int(i) for i in data[2].astype(int)])
    return sfr.PsiData(x,N,k,nafc)

def make_pmf(dataset, nafc, sigmoid, core, priors):
    sigmoid = get_sigmoid(sigmoid)
    core = get_core(core, dataset, sigmoid)
    pmf = sfr.PsiPsychometric(nafc, core, sigmoid)
    nparams = pmf.getNparams()
    set_priors(pmf,priors)
    return pmf, nparams

def make_dataset_and_pmf(data, nafc, sigmoid, core, priors):
    dataset = make_dataset(data, nafc)
    pmf, nparams = make_pmf(dataset, nafc, sigmoid, core, priors)
    return dataset, pmf, nparams

def get_sigmoid(descriptor):
    """ convert string represnetation of sigmoid to PsiSigmoid object """
    if not sig_dict.has_key(descriptor):
        raise PsignifitException("The sigmoid \'"+str(descriptor)+"\' you requested, is not available.")
    return sig_dict[descriptor]()

def get_core(descriptor, data, sigmoid):
    """ convert string representation of core to PsiCore object """
    descriptor, parameter = re.match('([a-z]+)([\d\.]*)', descriptor).groups()
    sigmoid_type = sigmoid.getcode()
    if descriptor not in core_dict.keys():
        raise PsignifitException("The core \'"\
                +str(descriptor)\
                +"\' you requested, is not available.")
    if len(parameter) > 0:
        return core_dict[descriptor](data, sigmoid_type, float(parameter))
    else:
        return core_dict[descriptor](data, sigmoid_type)

def get_prior(prior):
    """ convert string based representation of prior to PsiPrior object """
    try:
        prior = "sfr."+"Prior(".join(prior.split('('))
        return eval(prior)
    except Exception, e:
        return None

def set_priors(pmf, priors):
    if priors is not None:
        nparams = pmf.getNparams()
        if len(priors) != nparams:
            raise PsignifitException("You specified \'"+str(len(priors))+\
                    "\' priors, but there are \'"+str(nparams)+ "\' parameters.")
        for (i,p) in enumerate((get_prior(p) for p in priors)):
            if p is not None:
                pmf.setPrior(i, p)

def get_start(start, nparams):
    if len(start) != nparams:
            raise PsignifitException("You specified \'"+str(len(start))+\
                    "\' starting value(s), but there are \'"+str(nparams)+ "\' parameters.")
    else:
        return sfr.vector_double(start)

def get_params(params, nparams):
    if len(params) != nparams:
                raise PsignifitException("You specified \'"+str(len(start))+\
                        "\' parameters, but the model has \'"+str(nparams)+ "\' parameters.")
    else:
        return sfr.vector_double(params)

def get_cuts(cuts):
    if cuts is None:
        return sfr.vector_double([0.5])
    elif op.isSequenceType(cuts) and np.array([op.isNumberType(a) for a in cuts]).all():
        return sfr.vector_double(cuts)
    elif op.isNumberType(cuts):
        return sfr.vector_double([cuts])
    else:
        raise PsignifitException("'cuts' must be either None, a number or a "+\
                "sequence of numbers.")

