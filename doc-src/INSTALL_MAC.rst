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


Testing your installation
=========================

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

