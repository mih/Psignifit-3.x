#!/usr/bin/env python
# encoding: utf-8
######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

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

setup ( name = "pypsignifit",
    version = "3.0beta",
    author = "Ingo Fründ & Valentin Haenel",
    description = "Statistical inference for psychometric functions",
    packages = ["pypsignifit"],
    ext_modules = [psipy] )

