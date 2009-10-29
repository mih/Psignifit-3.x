#!/usr/bin/env python
# encoding: utf-8

from distutils.core import setup, Extension

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
            "src/linalg.cc"]
        )

setup ( name = "pypsignifit",
        version = "2.0",
        author = "Ingo Fründ & Valentin Haenel",
        description = "Statistical inference for psychometric functions",
        packages = ["pypsignifit"],
        ext_modules = [psippmodule] )

