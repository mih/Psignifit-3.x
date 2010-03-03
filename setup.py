#!/usr/bin/env python
# encoding: utf-8

from distutils.core import setup, Extension
import numpy

psippmodule = Extension ( "_psipy",
        sources = [
            "src/bootstrap.cc",
            "src/core.cc",
            "src/data.cc",
            "src/mclist.cc",
            "src/mcmc.cc",
            "src/optimizer.cc",
            "src/psipy.cc",
            "src/psychometric.cc",
            "src/rng.cc",
            "src/sigmoid.cc",
            "src/special.cc",
            "src/linalg.cc",
            "src/prior.cc"],
            include_dirs=[numpy.get_include()]
        )

setup ( name = "pypsignifit",
        version = "3.0",
        author = "Ingo Fründ & Valentin Haenel",
        description = "Statistical inference for psychometric functions",
        packages = ["pypsignifit"],
        ext_modules = [psippmodule] )

