Psignifit for matlab users
==========================

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
folder. If you now call::

    savepath

you avoid having to call the above command everytime you start matlab.

You can check that everything went fine by calling::

    test_psignifit


