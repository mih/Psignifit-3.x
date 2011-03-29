========================
Welcome to Psignifit 3.0
========================

Psignifit is a toolbox that allows the user to fit psychometric functions and to test
hypotheses about psychometric data. Compared to the "classical" version of psignifit,
the new package comes with a number of **additional features**:
    * full bayesian treatment of psychometric functions including bayesian model selection and goodness of fit assessment
    * identification of influential observations
    * detection of outliers (potentially to be excluded)
    * new philosophy of defining the shape of the psychometric function allowing for considerably more flexibility in specifiing psychometric function models.

The primary interface for psignifit 3.0 is now `python <http://www.python.org/>`_, instead of matlab. In contrast to
matlab, python is a complete programming language that supports virtually all features you
might wish for. Python allows you to perform numerical computations as flexible and fast as
in matlab. In addition, python provides features like object oriented programming and simple creation of graphical user interfaces. Python is easy to learn: even users with no prior programming experience can master python within weeks.
Finally python is free to use because of its OSI-approved open source license.

However, even though python is similar to *matlab*, it is not the same. Therefore, we also plan
to provide a matlab version of psignifit 3.0 with the first official release. However, we will
not guarantee support for this toolbox in future releases and we encourage users to use the
python version as this will be up to date.

We also noted that a growing number of statistical toolboxes are designed for the statistics
environment *R* and there might be users that are interested in using a R version of psignifit.
Similar to the matlab interface, we provide a basic R version of psignifit 3.0. Again, we do not
guarantee support for this R library. A very rudimentary R library can be found in the folder 'rpsignifit'. Be warned that this is work in progress.


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


****************
Getting in touch
****************

To contact the authors and current maintainers please use:

    psignifit-users@lists.sourceforge.net


This list can be used to ask questions about usage and installation, report
bugs, and request new features. If you use psignifit we recommend you subscribe
to this list.


**************
How to install
**************

Download the current version of psignifit from:

    `<http://sourceforge.net/projects/psignifit/>`_

In the following, installation instructions are provided for different operating systems.  You can find a detailed listing of `Required packages`_ at the end of this document, but the dependencies are also considered in the following sections:

* :ref:`operating_system_linux`
* :ref:`operating_system_mac`
* :ref:`operating_system_windows`

Currently this documentation only describes the installation of the python version of psignifit.
If you are looking for installation instructions for matlab, see `Installing mpsignifit`_. Installation instructions for R are going to follow as soon as these toolboxes are
ready for use.

The C++ core of psignifit does the real work. It is completely coded in C++ and does not require any
additional libraries installed.

For additional information about the structure of the code, the build system and
version control, see: :doc:`CONTRIBUTING`


.. _operating_system_linux:

Linux
=====

If you are using `Debian <http://www.debian.org/>`_, the following packages need to be installed:

* ``make``
* ``gcc``
* ``python``
* ``python-dev``
* ``python-numpy (provides python-numpy-dev)``
* ``python-scipy``
* ``python-matplotlib``
* ``python-sphinx``
* ``doxygen``
* ``python-nose``
* ``swig``

In order to check whether or not you have the packages already installed, type::

    sudo aptitude search make gcc python ... swig

Packages that are installed on your machine are listed with a leading <i>

In order to install missing packages, type::

    sudo aptitude install make gcc python ... swig



System-wide installation
------------------------
On the command line, navigate to the root directory of the psignifit distribution. Now, you can simply type::

    make install

as root and everything will be installed to the right place. By default, the psignifit documentation will be installed into the root directory into a folder called doc-html. To change this behavior, you might want to modify the Makefile (this should be self-explaining).

To generate the documentation use::

    make doc


If you want to try psignifit without installing it into your system, you are referred to the section `Execute without Installation`_ below.

If you want a special flavor of the python installation and are familiar with using python
setup-scripts, you can also use special options for the installation, by
executing the ``setup.py`` script explicitly. Note however, that in this case
you will first have to generate the swig interface. An example can be found in
the section `Install into users home directory`_.


Install into users home directory
---------------------------------
If you do not have root/admin rights on your computer the setup routine allows installation into your home-directory. In this case you must first generate the ``swig`` interface::

    make swig

After this you may install psignifit locally by typing::

    python setup.py install --home=$HOME

where ``$HOME`` should be replace by the name of your own home-directory.
This command will install psignifit into ``$HOME/lib/python/psignifit``.
To use psignifit from this path, you will also have to set the ``$PYTHONPATH``
variable. Either you invoke the python interpreter from the commandline by
calling::

    PYTHONPATH=$HOME/lib/python python

or you set the ``$PYTHONPATH`` variable in your ``.bashrc`` (or equivalent) file
by adding the line::

    export PYTHONPATH=$HOME/lib/python

Yet another option is to set the ``$PYTHONPATH`` variable directly from the
python interpreter using the ``os`` module.


Installing the command line interface
-------------------------------------

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



.. _operating_system_mac:

Mac OSX
=======

The easiest way to use psignifit is to use the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.
In that case, you might want to download one of the 'swigged' archives from:

`Psignifit3 Downloads <http://sourceforge.net/.projects/psignifit/files/>`_

(The 'swigged' archives are the ones with the string 'swigged' in the filename)

Extract the archive to a folder. Navigate to this folder from the command line and install::

    unzip psignifit3.0_beta_swigged_24-03-2011.zip
    cd psignifit3.0_beta_swigged_24-03-2011
    python setup.py install


Installing the command line interface
-------------------------------------

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



.. _operating_system_windows:

Windows
=======

The easiest way to use psignifit is to use the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.
In that case, you might want to download one of the 'swigged' archives from:

`Psignifit3 Downloads <http://sourceforge.net/.projects/psignifit/files/>`_

(The 'swigged' archives are the ones with the string 'swigged' in the filename)

Extract the archive to a folder. Navigate to this folder from the command line and install::

    unzip psignifit3.0_beta_swigged_24-03-2011.zip
    cd psignifit3.0_beta_swigged_24-03-2011
    python setup.py install


Installing the command line interface
-------------------------------------

Download the file ``psignifit_cli_3_beta_installer.exe`` from
`sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and run it.
Follow the instructions on the screen. At the end of the installation, you will be asked whether
you want to add psignifit-cli to your environment path. You should leave this box checked. You
will not be able to use psignifit from within matlab if you uncheck this box!


Testing your installation
=========================

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.



Execute without installation
============================

If you wish to build and execute pypsignifit in place, simply type::

    make

This will build everything into the current working directory, and allow you to
import psignifit as long as you remain in the current working directory.



Installing mpsignifit
=====================

In order to use psignifit from within matlab (mpsignifit), you have to install both, the command line interface (of your respective operating system) as well as mpsignifit.

Installing the mpsignifit files
-------------------------------

If you have not yet obtained a copy of the psignifit sources, get one from `sourceforge <http://sourceforge.net/projects/psignifit/files/>`_.
The file will most probably be a file ending either with ``.tar.gz`` or with ``.zip``.
Unpack the file and navigate to the unpacked folder. Within that folder there is (among
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



Required packages
-----------------

This section lists all dependencies. The version numbers are the versions we
used during development.

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

* `sphinx <http://sphinx.pocoo.org/>`_ (0.6.5-1)
    to generate the html documentation
* `doxygen <http://www.stack.nl/~dimitri/doxygen/>`_ (1.6.3-1)
   to generate the C++ documentation
* `epydoc <http://epydoc.sourceforge.net/>`_ (3.0.1-5)
   to generate the Python API documentation

Testing
.......

* `nosetest <http://somethingaboutorange.com/mrl/projects/nose/0.11.2/>`_ (0.11.1-1)
  for running some of the unit tests

