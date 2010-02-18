
#!/usr/bin/env python

"""
setup.py file for SWIG example
"""

from distutils.core import setup, Extension


swignifit = Extension('_swignifit',
        sources = [
            "swignifit_wrap.cxx",
            "../src/bootstrap.cc",
            "../src/core.cc",
            "../src/data.cc",
            "../src/mclist.cc",
            "../src/mcmc.cc",
            "../src/optimizer.cc",
            "../src/psychometric.cc",
            "../src/rng.cc",
            "../src/sigmoid.cc",
            "../src/special.cc",
            "../src/linalg.cc",
            "../src/prior.cc"],
            include_dirs=["../src"]
        )

setup (name = 'swignifit',
       version = '0.1',
       author      = "Valentin Haenel",
       description = """experimental swig wrapping for psignifit""",
       ext_modules = [swignifit],
       py_modules = ["swignifit"],
       )

