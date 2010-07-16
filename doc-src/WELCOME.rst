========================
Welcome to Psignifit 3.0
========================

Psignifit is a toolbox that allows the user to fit psychometric functions and to test
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
However, even though python is similar to matlab, it is not the same. Therefore, we also plan
to provide a matlab version of psignifit 3.0 with the first official release. However, we will
not guarantee support for this toolbox in future releases and we encourage users to use the
python version as this will be up to date. python also offers a further advantage to matlab:
python is free [1]_.

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

Getting in touch
================

To contact the authors and current maintainers please use:

    psignifit-users@lists.sourceforge.net


This list can be used to ask questions about usage and installation, report
bugs, and request new features. If you use psignifit we recommend you subscribe
to this list.

How to install
==============

If you want to install psignifit on your computer make sure that you have all the dependencies installed.
Currently this documentation only deals with the installation of the python version of psignifit.
Installation instructions for matlab and R are going to follow as soon as these toolboxes are
ready for use.

The C++ core of psignifit does the real work. It is completely coded in C++ and does not require any
additional libraries installed.

Download the current version of psignifit from:

    `<http://sourceforge.net/projects/psignifit/>`_

For additional information about the structure of the code, the build system and
version control, see: :doc:`CONTRIBUTING`

Dependencies
------------

This section lists all dependencies. The version numbers are the versions we
used during development.

When using Debian
.................

If you are using `Debian <http://www.debian.org/>`_, the following packages will
satisfy the dependencies listed below:

* ``make``
* ``gcc``
* ``python``
* ``python-dev``
* ``python-numpy (provides python-numpy-dev)``
* ``swig``
* ``python-scipy``
* ``python-matplotlib``
* ``python-sphinx``
* ``doxygen``
* ``python-nose``

Compile-Time
............
* `Make <http://www.gnu.org/software/make/>`_ (3.81-8)
    for building the software
* `Gcc <http://gcc.gnu.org/>`_ (4:4.4.3-1)/
* `Python <python http://www.python.org/>`_ (2.5.5-6)/
* `Python/C API <http://docs.python.org/c-api/>`_ (2.5.5-2)
    for compiling the ``psipy`` and ``swignifit`` interface for python
* `Numpy/C API <http://docs.scipy.org/doc/numpy/reference/c-api.html>`_ (1:1.3.0-3)
    for compiling the ``psipy`` interface for python
* `Simplified Wrapper and Interface Generator (SWIG) <http://www.swig.org/>`_ (1.3.40-2)
    for compiling the ``swignifit`` interface for python

Run-Time
........
* `Python <python http://www.python.org/>`_ (2.5.5-6)/
* `Numpy <http://numpy.scipy.org/>`_  (1:1.3.0-3)/
* `Scipy <http://www.scipy.org/>`_ (0.7.1-1)/
* `Matplotlib <http://matplotlib.sourceforge.net/>`_ (0.99.1.2-3)
    for the python version

Documentation
.............

`sphinx <http://sphinx.pocoo.org/>`_ (0.6.5-1)
    to generate the html documentation
`doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ (1.6.3-1)
   to generate the C++ documentation
`epydoc <http://epydoc.sourceforge.net/>`_ (3.0.1-5)
   to generate the Python API documentation

Testing
.......

`nosetest <http://somethingaboutorange.com/mrl/projects/nose/0.11.2/>`_
(0.11.1-1)
  for running some of the unit tests


Linux
-----

On the command line, navigate to the root directory of the psignifit distribution. By default,
the installation process will install the psignifit documentation into the root directory into
a folder called doc-html . To change this behavior, you might want to modify the Makefile (this
should be self-explaining). Now, you can simply type::

    make install

as root and everything will be installed to the right place. If you don't want the documentation
be installed, you can also say::

    make python-install

which will internally call::

    python setup.py install

So, if you want a special flavor of the python installation and are familiar with using python
setup-scripts, you can also use special options for the installation.

In some cases, you may want to install psignifit locally in your users home
directory. For details about this, see `Install into users home directory`_.

Mac OSX
-------

A simple::

    python setup.py install

should install the python toolbox for you. However, keep in mind that you need the abovementioned
dependencies.

In some cases, you may want to install psignifit locally in your users home
directory. For details about this, see `Install into users home directory`_.

Install into users home directory
---------------------------------

In some cases, you do not have root/admin rights on the computer you are working
on. This would prevent you from installing psignifit in the system path as
described above. As a workaround, the setup routine allows installation into a
users home-directory using the command::

    python setup.py install --home=$HOME

where you should replace ``$HOME`` with the name of your own home-directory.
This command will install psignifit into ``$HOME/lib/python/psignifit``.
To use psignifit from this path, you will also have to set the ``$PYTHONPATH``
variable. Either you invoke the python interpreter from the commandline by
calling::

    PYTHONPATH=$HOME/lib/python python

or you set the ``$PYTHONPATH`` variable in your ``.bashrc`` (or equivalent) file
by adding the line::

    export PYTHONPATH=$HOME/lib/python

The last way to set the ``$PYTHONPATH`` variable is to set it directly from the
python interpreter using the ``os`` module.

Execute without Installation
----------------------------

If you wish to build and execute pypsignifit in place, you must add the results
of the build process to the ``$PYTHONPATH``.

Build with::

    python setup.py build

Execute with::

    PYTHONPATH=build/lib.macosx-10.3-i386-2.5 python -c "import pypsignifit"

But remember to replace ``lib.macosx-10.3-i386-2.5`` with whatever is appropriate to
your operating system. You will find this in the ``build`` directory.

.. [1] That means both, free as in "free beer" and free as in "free speech".
