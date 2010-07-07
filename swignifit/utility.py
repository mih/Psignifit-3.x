#/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

"""Variety of utilities for working with swig generated code."""

__docformat__ = "restructuredtext"

import operator as op
import re
import numpy as np
import swignifit_raw as sfr

class PsignifitException(Exception):
    pass

def extract_subclasses(base, sub_func):
    """Recursively extract subclasses, given a swig base class.

    Parameters
    ----------
    base : swig class
        The base class from which to start.
    sub_func : string
        The function or attribute to use as name for subclass.

    Returns
    -------
    A dictionary mapping subclass names to constructors.

    """
    to_visit = base.__subclasses__()
    subclasses = dict()
    for cl in to_visit:
        descriptor = eval("cl."+sub_func)
        if descriptor not in subclasses.keys():
            subclasses[descriptor] = cl
            to_visit.extend(cl.__subclasses__())
    return subclasses

def extract_subclasses_descriptor(base):
    """Recursively extract subclasses, using the `getDescriptor()` method."""
    return extract_subclasses(base, "getDescriptor()")

def extract_subclasses_reflection(base):
    """Recursively extract subclasses, using the `__name__` attribute."""
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
    """Create a PsiData object from column based input.

    Parameters
    ----------
    data : sequence on length 3 sequences
        Psychometric data in colum based input,
        e.g.[[1, 1, 5], [2, 3, 5] [3, 5, 5]].
    nafc : int

    Returns
    -------
    PsiData object

    """
    data = np.array(data).T
    x = sfr.vector_double(data[0])
    k = sfr.vector_int([int(i) for i in data[1].astype(int)])
    N = sfr.vector_int([int(i) for i in data[2].astype(int)])
    return sfr.PsiData(x,N,k,nafc)

def make_pmf(dataset, nafc, sigmoid, core, priors):
    """Assemble PsiPsychometric object from model parameters.

    Parameters
    ----------
    nafc : int
    sigmoid : string
        Description of model sigmoid.
    core : string
        Description of model core.
    priors : sequence of strings
        The model priors.

    Returns
    -------
    (PsiPsychometric, int)
    PsiPsychometric object and number of free parameters in model.

    """
    sigmoid = get_sigmoid(sigmoid)
    core = get_core(core, dataset, sigmoid)
    pmf = sfr.PsiPsychometric(nafc, core, sigmoid)
    nparams = pmf.getNparams()
    set_priors(pmf,priors)
    return pmf, nparams

def make_dataset_and_pmf(data, nafc, sigmoid, core, priors):
    """Assemble PsiData and PsiPsychometric objects simultaneously.

    Parameters
    ----------
    see: make_dataset and make_pmf

    Returns
    -------
    (PsiData, PsiPsychometric, int)
    PsiData object, PsiPsychometric object and number of free parameters.

    """
    dataset = make_dataset(data, nafc)
    pmf, nparams = make_pmf(dataset, nafc, sigmoid, core, priors)
    return dataset, pmf, nparams

def get_sigmoid(descriptor):
    """Convert string representation of sigmoid to PsiSigmoid object."""
    if not sig_dict.has_key(descriptor):
        raise PsignifitException("The sigmoid \'"+str(descriptor)+"\' you requested, is not available.")
    return sig_dict[descriptor]()

def get_core(descriptor, data, sigmoid):
    """Convert string representation of core to PsiCore object."""
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
    """Convert string based representation of prior to PsiPrior object."""
    try:
        prior = "sfr."+"Prior(".join(prior.split('('))
        return eval(prior)
    except Exception, e:
        return None

def set_priors(pmf, priors):
    """Set the priors to be used in the model object.

    Parameters
    ----------
    pmf : PsiPsychometric object
        the model
    priors : list of strings of length of free parameters of `pmf`
        list of priors

    """
    if priors is not None:
        nparams = pmf.getNparams()
        if len(priors) != nparams:
            raise PsignifitException("You specified \'"+str(len(priors))+\
                    "\' priors, but there are \'"+str(nparams)+ "\' parameters.")
        for (i,p) in enumerate((get_prior(p) for p in priors)):
            if p is not None:
                pmf.setPrior(i, p)

def get_start(start, nparams):
    """Convert sequence of starting values to vector_double type."""
    if len(start) != nparams:
            raise PsignifitException("You specified \'"+str(len(start))+\
                    "\' starting value(s), but there are \'"+str(nparams)+ "\' parameters.")
    else:
        return sfr.vector_double(start)

def get_params(params, nparams):
    """Convert sequence of parameter values to vector_double type."""
    if len(params) != nparams:
                raise PsignifitException("You specified \'"+str(len(start))+\
                        "\' parameters, but the model has \'"+str(nparams)+ "\' parameters.")
    else:
        return sfr.vector_double(params)

def get_cuts(cuts):
    """ Convert `cuts` argument to vector_double type.

    Argument can be None, a number or a sequence of numbers. If None, there is
    only one cut at 0.5. If `cuts` is a number, function returns a vector_double
    with that number as a single element. If its a sequence, that sequence will
    be converted to vector_double type.

    Parameters
    ----------
    cuts : None, number or sequence of numbers

    Returns
    -------
    vector_double

    """
    if cuts is None:
        return sfr.vector_double([0.5])
    elif op.isSequenceType(cuts) and np.array([op.isNumberType(a) for a in cuts]).all():
        return sfr.vector_double(cuts)
    elif op.isNumberType(cuts):
        return sfr.vector_double([cuts])
    else:
        raise PsignifitException("'cuts' must be either None, a number or a "+\
                "sequence of numbers.")

