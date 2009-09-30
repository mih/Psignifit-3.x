=============
psignifit API
=============

Data Types
==========
.. automodule:: psignidata

.. autoclass:: BootstrapInference
   :members:

.. autoclass:: BayesInference
   :members:

Diagnostic Plots
================

Default functions
-----------------

The following functions will be imported by default:

.. automodule:: psigniplot
   :members:

Subfunctions
------------

To get access to all plot functions in isolation, they can also be imported separately. Here is the documentation

.. autofunction:: psigniplot.drawaxes
.. autofunction:: psigniplot.plotRd
.. autofunction:: psigniplot.plotHistogram
.. autofunction:: psigniplot.plotPMF
.. autofunction:: psigniplot.plotThres
.. autofunction:: psigniplot.plotGeweke
.. autofunction:: psigniplot.plotChains
.. autofunction:: psigniplot.plotParameterDist

Errors
======

.. automodule:: psignierrors
    :members:
