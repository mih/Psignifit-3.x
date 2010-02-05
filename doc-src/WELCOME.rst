========================
Welcome to psignifit 3.0
========================

psignifit is a toolbox that allows the user to fit psychometric functions and to test
hypotheses about psychometric data. Compared to the "classical" version of psignifit,
the new package comes with a number of additional features for fitting psychometric functions.
These include the detection of outliers, and the identification of influential
observations, and a full bayesian treatment of psychometric functions including bayesian
model selection and goodness of fit assessment. A new philosophy of defining the
shape of the psychometric function allows for considerably more flexibility in specifiing
psychometric function models.

The primary interface for psignifit 3.0 is now python, instead of matlab. In contrast to
matlab, python is a complete programming language that supports virtually all features you
might wish for. python allows you to perform numerical computations as flexible and fast as
in matlab. In addition, python provides a number of additional features like object
oriented programming and simple creation of graphical user interfaces. python is easy to
learn: even users with no prior programming experience often master python within two weeks.
However, even though python is similar to matlab, it is not the same. Therefore, we also
provide a matlab version of psignifit 3.0. However, we do not guarantee support for this
toolbox in future releases and we encourage users to use the python version as this will be up
to date. python also offers a further advantage to matlab: python is free [1]_.

TODO: matlab toolbox

We also noted that a growing number of statistical toolboxes are designed for the statistics
environment R and there might be some users that are interested in using a R version of psignifit.
Similar to the matlab interface, we provide a basic R version of psignifit 3.0. Again, we do not
guarantee support for this R library.

A very rudimentary R library can be found in the folder 'rpsignifit'. Be warned that this is work
in progress. The major aim of current development is to have a proper version of psignifit running
at all. This version will be in python as mentioned above. However, after the python version is
running properly, there will also be an R version.

We would be glad to find developers that are interested in supporting these non-python interfaces
to the psignifit engine.


How to install
==============

If you want to install psignifit on your computer make sure that you have all the dependencies installed.
Currently this documentation only deals with the installation of the python version of psignifit.
Installation instructions for matlab and R are going to follow as soon as these toolboxes are
ready for use.

The C++ core of psignifit does the real work. It is completely coded in C++ and does not require any
additional libraries installed.

The python version of psignifit requires the following packages installed:

    * `numpy <http://numpy.scipy.org/>`_ / `scipy <http://www.scipy.org/>`_
    * `matplotlib <http://matplotlib.sourceforge.net/>`_

Download the current version of psignifit from

`<http://sourceforge.net/projects/psignifit/>`_

To also build the documentation yourself, you will need `sphinx <http://sphinx.pocoo.org/>`_.

Linux
-----

On the command line, navigate to the root directory of the psignifit distribution. By default,
the installation process will install the psignifit documentation into the root directory into
a folder called doc-html . To change this behavior, you might want to modify the Makefile (this
should be self-explaining). Now, you can simply type::

    make install

as root and everythin will be installed to the right place. If you don't want the documentation
be installed, you can also say::

    make python-install

which will internally call::

    python setup.py install

So, if you want a special flavor of the python installation and are familiar with using python
setup-scripts, you can also use special options for the installation.

Mac OSX
-------

A simple::

    python setup.py install

should install the python toolbox for you. However, keep in mind that you need the abovementioned
dependencies.

.. [1] That means both, free as in "free beer" and free as in "free speech".
