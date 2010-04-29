#!/usr/bin/env python
# encoding: utf-8
# vi: set ft=python sts=4 ts=4 sw=4 et:

######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

""" setup: part of the build system for psignifit 3.x

This file will build the pypsignifit module, and compile ONLY the psipy
extension.

"""

from setup_basic import *

setup ( name = name,
    version = version,
    author = author,
    description = description,
    packages = packages,
    ext_modules = [psipy] )

