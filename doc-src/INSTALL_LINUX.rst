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


Testing your installation
=========================

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

