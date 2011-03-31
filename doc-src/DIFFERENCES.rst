=================================
Changes between psignifit 2 and 3
=================================

This page contains some of the changes in conversion between psignifit 2 and 3. As of now, this list is nots complete, while we are in the process of updating this list, please let us know if you find any differences that we have not added so far.

Order of Parameters
-------------------

In psignifit 2 the different parameters were used in the following oder:

a b gamma lambda

Psignifit 3 changes the order of the parameters slightly when you are using the classical notation 

a b lambda gamma

This is because you do not have to use gamma in all situations. When you are analysing a 2AFC task you will only have to specify 3 priors. If you are analysing a yes-no task, for example, you will have to specify your 4th prior for gamma as well.


SPecifying Psychometric function
Old: give function (one of the options)
New: now, simoid and core specify what maps onto what
OLD: always ab
NEW: still an option, now mw is also a option, and lots of different other combinations.

Priors specified differently
NEW: struct, if you want to change one thing you have to make all explicit either completely default or all specified
OLD: each parameter separately, so you can keep most of them default and only change one


