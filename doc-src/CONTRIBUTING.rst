=============================
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
    The mwCore class scales the second parameter such that w can be directly interpreted as the
    width of the region in which the sigmoid still rises significantly. What to "rise significantly"
    means is parameterized by the parameter alpha of the mwCore. The default alpha==0.1 indicates
    that w is the width of the range over which the sigmoid rises from 0.1 to 0.9. Thus, the scaling
    of the second parameter obviously depends on the sigmoid.
    The constructor for the mwCore class looks roughly like this::

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

    At the position marked by::

        /////////////// here ////////////////

    in the above code example, you should add a new case that defines all the scaling parameters
    depending on your sigmoid. zalpha scales w to the correct range, zshift is an additional
    shift to ensure the the sigmoid has an output value of 0.5 at an input value of 0.
