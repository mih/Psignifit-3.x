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

TODO: R library

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


Contributing to psignifit 2.0
=============================

psignifit 2.0 is free software. You are free to modify the software under the terms of the license
that is distributed with psignifit 2.0. We welcome developers that want to contribute to psignifit 2.0.
A development snapshot of psignifit can be obtained like this

TODO: command to clone git repo

Commits
-------

To make it easier to keep track of the development of psignifit, we use a number of marks for commits.
Use the following marks for commits:

[NF]    new feature
[BF]    bug fix
[RF]    refactoring
[FO]    code formatting (adding spaces etc.)
[UT]    unit tests
[DOC]   documentation

Extending
---------

In principle every part of the library can be replaced. This is generally done by deriving from the fundamental base classes.
An exception is adding a new sigmoid:

Adding a new sigmoid
....................

Adding a new sigmoid requires two steps:
1. Write a new class that inherits from PsiSigmoid
2. If you want your new sigmoid to work with mwCore objects, you have to add a label for that, too.
   The constructor for the mwCore class looks roughly like this:

    mwCore::mwCore ( int sigmoid, double al )
            : sigmtype(sigmoid), alpha(al), zshift(0) {
        switch (sigmoid) {
        case 1:
            ...
            break;
        /////////////// here ////////////////
        default:
            throw NotImplementedError();
        }
    }

    The mwCore class scales the second parameter such that w can be directly interpreted as the
    width of the region in which the sigmoid still rises significantly. What to "rise significantly"
    means is parameterized by the parameter alpha of the mwCore. The default alpha==0.1 indicates
    that w is the width of the range over which the sigmoid rises from 0.1 to 0.9. Thus, the scaling
    of the second parameter obviously depends on the sigmoid. At the position marked by
        /////////////// here ////////////////
    in the above code example, you should add a new case that defines all the scaling parameters
    depending on your sigmoid. zalpha scales w to the correct range, zshift is an additional
    shift to ensure the the sigmoid has an output value of 0.5 at an input value of 0.
