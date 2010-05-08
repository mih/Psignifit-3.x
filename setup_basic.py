#!/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" part of the build system for psignifit 3.x

Provides definitions commmon to all setup files.

"""

from distutils.core import setup, Extension
import numpy

psipp_sources = [
    "src/bootstrap.cc",
    "src/core.cc",
    "src/data.cc",
    "src/mclist.cc",
    "src/mcmc.cc",
    "src/optimizer.cc",
    "src/psychometric.cc",
    "src/rng.cc",
    "src/sigmoid.cc",
    "src/special.cc",
    "src/linalg.cc",
    "src/prior.cc"]

psipy_sources = ["psipy/psipy.cc"]
psipy = Extension ( "_psipy",
    sources = psipp_sources + psipy_sources,
    include_dirs=[numpy.get_include(), "src", "psipy"])

swignifit_sources = ["swignifit/swignifit_raw.cxx"]
swignifit = Extension('swignifit._swignifit_raw',
        sources = psipp_sources + swignifit_sources,
        include_dirs=["src"])

name = "pypsignifit"
version = "3.0beta"
author = "Ingo Fründ & Valentin Haenel"
description = "Statistical inference for psychometric functions"
packages = ["pypsignifit"]

