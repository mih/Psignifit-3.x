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

.. raw :: html

    <div class="admonition-see-also admonition seealso">
    <p class="first admonition-title">Important</p>
    <p class="last">
    Psignifit3 is beta software. This means that it is still under heavy development. Occasionally
    it might not work as expected. This might have two reasons:
    either, the documentation was ambiguous, or you have discovered a programming error (bug).
    In any case, you would help us a lot if you write an email to our mailing list (see below)
    or <a href="mailto:igordertigor@users.sourceforge.net,valentin-haenel@users.sourceforge.net?subject=[psignifit]">personally to us</a>.
    Typically we can solve your problem within hours.
    </p>
    </div>

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
If you are looking for installation instructions for matlab, see `Installing mpsignifit and the command line interface`_.
Installation instructions for R are going to follow as soon as these toolboxes are
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

When using EPD
..............

If you are using Windows or MacOS, the easiest way to use psignifit is to use the
`Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.
In that case, you might want to download one of the 'swigged' archives.

Extract the archive to a folder. Navigate to this folder from the command line and type::

    python setup.py install

This should give you the python version of psignifit (and only the python version!).

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

`nosetest <http://somethingaboutorange.com/mrl/projects/nose/0.11.2/>`_ (0.11.1-1)
  for running some of the unit tests


Linux and Mac OSX
-----------------

On the command line, navigate to the root directory of the psignifit distribution. By default,
the installation process will install the psignifit documentation into the root directory into
a folder called doc-html . To change this behavior, you might want to modify the Makefile (this
should be self-explaining). Now, you can simply type::

    make install

as root and everything will be installed to the right place.

To generate the documentation use::

    make doc

If you want to try psignifit without installing it into your system, you might
consider reading the section `Execute without Installation`_ below.

If you want a special flavor of the python installation and are familiar with using python
setup-scripts, you can also use special options for the installation, by
executing the ``setup.py`` script explicitly. Note however, that in this case
you will first have to generate the swig interface. An example can be found in
the section `Install into users home directory`_.

Install into users home directory
---------------------------------

In some cases, you do not have root/admin rights on the computer you are working
on. This would prevent you from installing psignifit in the system path as
described above. As a workaround, the setup routine allows installation into a
users home-directory. In this case you must first generate the ``swig``
interface::

    make swig

After this you may install psignifit locally by typing::

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

If you wish to build and execute pypsignifit in place, simply type::

    make

This will build everything into the current working directory, and allow you to
import psignifit as long as you remain in the current working directory.

Testing your installation
-------------------------

To run a number of tests on your installation, you can call::

    make test

This will call a rather large suite of tests for psignifit.

Installing mpsignifit and the command line interface
----------------------------------------------------

mpsignifit is a matlab version for psignifit. As mentioned above, psignifit is developed as a
python tool, thus in most cases the python version will be more up to date and have more features.
There were technical reasons to switch the development of psignifit from matlab to python.
To overcome these technical problems, we had to separate the workhorse functions of psignifit
completely from the matlab environment. Psignifit now comes with a very rudimentary command line
interface. The matlab version of psignifit, mpsignifit will then internally call commands from
the command line and integrate the results in matlab. This means that in order to use psignifit
from within matlab, you have to install both, the command line interface as well as mpsignifit.
This section describes how to do so.

Be aware however that the matlab version of psignifit provides significantly less features than
the python version.

Installing the command line interface on Mac OSX or Linux
.........................................................

Download psignifit from `sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and
extract the compressed file to a folder in your home directory. Navigate into the folder.
You have two installation options. By default, the command line interface will be installed to a
folder called ``bin`` in your home directory. You can change this behavior by editing the
``Makefile``. At the beginning of the ``Makefile``, you find a line::

    CLI_INSTALL=$(HOME)/bin

replace this by e.g. ``/usr/bin/`` for system wide installation.

Once you have the Makefile in your desired shape type::

    make cli-install

If the installation directory is not on your system search path, you may have to add it.
To do so, add::

    export PATH=$PATH:$HOME/bin

to your ``.bashrc`` (if you use bash). If you use zsh, the same line should be in your
``.zshrc.local`` file.

Now, you should be able to call::

    psignifit-mcmc -h
    psignifit-diagnostics -h
    psignifit-bootstrap -h
    psignifit-mapestimate -h

And see some usage messages after each call.

Installing the command line interface on Windows
................................................

Download the file ``psignifit_cli_3_beta_installer.exe`` form
`sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and run it.
Follow the instructions on the screen. At the end of the installation, you will be asked whether
you want to add psignifit-cli to your environment path. You should leave this box checked. You
will not be able to use psignifit from within matlab if you uncheck this box!

Installing the matlab files
...........................

If you have not yet obtained a copy of the psignifit sources, get one now (see above).
The file will most probably be a file ending either with ``.tar.gz`` or with ``.zip``.
Unpack the file and navigate to the unpacked folder. Within that folder there is (amoung
other things) one folder called ``mpsignifit``. Copy this folder to a save place (e.g. the
``toolbox`` folder in your matlab installation directory).
Now you have to make matlab aware that the new files are there. To do so, start matlab and
type::

    addpath path\to\mpsignfit\files

where you replace ``path\to\mpsignifit\files`` with the path where you copied the ``mpsignifit``
folder. You might now want to call::

    savepath

to avoid having to call the above command everytime you start matlab.

You can check that everything went fine by calling::

    test_psignifit



.. [1] That means both, free as in "free beer" and free as in "free speech".
