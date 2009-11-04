=============
psignifit API
=============

Data Types
==========
.. automodule:: pypsignifit.psignidata

.. autoclass:: BootstrapInference
   :members:

.. autoclass:: BayesInference
   :members:

Diagnostic Plots
================

Default functions
-----------------

The following functions will be imported by default:

.. automodule:: pypsignifit.psigniplot
   :members:

Subfunctions
------------

To get access to all plot functions in isolation, they can also be imported separately. Here is the documentation

.. autofunction:: pypsignifit.psigniplot.drawaxes
.. autofunction:: pypsignifit.psigniplot.plotRd
.. autofunction:: pypsignifit.psigniplot.plotHistogram
.. autofunction:: pypsignifit.psigniplot.plotPMF
.. autofunction:: pypsignifit.psigniplot.plotThres
.. autofunction:: pypsignifit.psigniplot.plotGeweke
.. autofunction:: pypsignifit.psigniplot.plotChains
.. autofunction:: pypsignifit.psigniplot.plotParameterDist

Simulated Observers
===================

psignifit allows to simulate a number of observers to access stability of psychometric functions.

.. automodule:: pypsignifit.psigobservers

.. autoclass:: Observer
   :members:

.. autoclass:: LinearSystemLearner
   :members:

.. autoclass:: CriterionSettingObserver
   :members:

.. autoclass:: BetaBinomialObserver
   :members:

Errors
======

.. automodule:: pypsignifit.psignierrors
    :members:
