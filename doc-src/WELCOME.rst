========================
Welcome to psignifit 2.0
========================

psignifit is a toolbox to fit psychometric functions an test hypotheses on psychometric data.
Compared to the old "classical" version of psignifit, this package a number of  new features
for fitting psychometric functions. These include detection of outliers and influential
observations and a full bayesian treatment of psychometric functions that includes bayesian
model selection and goodness of fit assessment. Furthermore, a new philosophy of defining the
shape of the psychometric function allows for considerably more flexibility in specifiing
psychometric function models.

The primary interface for psignifit 2.0 is no longer matlab but is python. In contrast to
matlab, python is a complete programming language that supports virtually all features you
might wish for. python allows you to perform numerical computations as flexible and fast as
in matlab. In addition, python allows for a number of additional features like object
oriented programming and simple creation of graphical user interfaces. python is easy to
learn: even users with no prior programming experience often master python within two weeks,
python is similar to matlab but it is not the same. Therefore, we also provide a matlab
version of psignifit 2.0. However, we do not guarantee support for this toolbox in future
releases and we encourage users to use the python version as this will usually be more up
to date. python also offers a further advantage to matlab: python is free (and that means
both, free as in "free beer" and free as in "free speech").

TODO: matlab toolbox

We also noted that a growing number of statistical toolboxes are designed for the statistics
environment R and there might be some users that are interested in using a R version of psignifit.
Similar to the matlab interface, we provide a basic R version of psignifit 2.0. Again, we do not
guarantee support for this R library.

A very rudimentary R library is found in the folder rpsignifit. this is definitely work in progress.
The current aim of development is mainly to get a proper version of psignifit running
at all. This version will be in python as mentioned above. However, After the python version is
running properly, there will also be an R version.

We would be glad to find developers that are interested in supporting these non-python interfaces
to the psignifit engine.


How to install
==============

Linux
-----

On the command line, navigate to the root directory of the psignifit distribution. By default,
the installation process will install the psignifit documentation into the root directory into
a folder called doc-html . To change this behavior, you might want to modify the Makefile (this
should be self-explaining). Now, you can simply type

make install

as root and everythin will be installed to the right place. If you don't want the documentation
be installed, you can also say

make python-install

which will internally call

python setup.py install

So, if you want a special flavor of the python installation and are familiar with using python
setup-scripts, you can also use special options for the installation.

Mac OSX
-------

TODO: Valentin, you did that, didn't you?

