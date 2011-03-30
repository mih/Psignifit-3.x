Windows
=======

This is the installation instruction for the python version, not the matlab version.

The easiest way to use psignifit is to use the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.
In that case, you might want to download one of the 'swigged' archives from:

`Psignifit3 Downloads <http://sourceforge.net/.projects/psignifit/files/>`_

(The 'swigged' archives are the ones with the string 'swigged' in the filename)

Extract the archive to a folder. Navigate to this folder from the command line and install::

    unzip psignifit3.0_beta_swigged_24-03-2011.zip
    cd psignifit3.0_beta_swigged_24-03-2011
    python setup.py install


Installing the command line interface (optional)
-------------------------------------------------

Download the file ``psignifit_cli_3_beta_installer.exe`` from
`sourceforge <http://sourceforge.net/projects/psignifit/files/>`_ and run it.
Follow the instructions on the screen. At the end of the installation, you will be asked whether
you want to add psignifit-cli to your environment path. You should leave this box checked. You
will not be able to use psignifit from within matlab if you uncheck this box!


Testing your installation
-------------------------

To check whether your installation has been successful and pypsignifit is working properly, you can call::

    make test

This will call the standard test suite for psignifit.

