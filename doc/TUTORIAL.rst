============================
A quick start to pypsignifit
============================

This document presents two example analyzes of psychometric function data using pypsignifit.
The first example explains how to fit a psychometric function using constrained maximum
likelihood as described in the papers by Wichmann and Hill (2001a,b). The second example
deals with a bayesian approach to fitting a psychometric function. Parts of the ideas for
this can be found in the paper by Kuss et al (2005), however most of this example is new
at the time of this writing.

To get you starting with pypsignifit, open a python interpreter and type the following:

>>> from pypsignifit import *
>>> dir()
['BayesInference', 'BootstrapInference', 'ConvergenceMCMC', 'GInitiallyUnowned', 'GoodnessOfFit', 'ParameterPlot', '__builtins__', '__doc__', '__name__', 'show']

As you see, there is a number of functions and data types imported in the current workspace.
To view documentation about one of these functions, you can use the online python help by typing
help ( function name ). For instance,

>>> help ( BayesInference )

will show the documentation of the BayesInference object class.

We will now create an example data set for which we want to estimate a psychometric function.
We assume that the data are from a 2AFC task

>>> nafc = 2
>>> stimulus_intensities = [0.0,2.0,4.0,6.0,8.0,10.0]
>>> number_of_correct = [34,32,40,48,50,48]
>>> number_of_trials  = [50]*len(stimulus_intensities)
>>> data = zip(stimulus_intensities,number_of_correct,number_of_trials)

The last line creates a list of tuples of the form (stimulus intensity, number of correct responses,
number of trials). Each tuple summarizes data from a single experimental block. We will assume that
the data have been acquired in the same sequence in which they are entered, i.e. in the sequence
of ascending stimulus intensity.

Example 1: Constrained Maximum Likelihood and Bootstrap Inference
=================================================================

Constrained maximum likelihood provides a way to estimate parameters from a psychometric function
using maximum likelihood estimation while imposing constraints on some of the parameters.
Typically, the psychometric function is parameterized as

Psi(x) = guess + (1-guess-lapse) * F ( x | a,b ),

where guess is the guessing rate, lapse is the lapsing rate, F is a sigmoid function and a and
b are parameters governing the shape of the sigmoid function. In this example, we will try to fit
the above data with a logistic function, using the same parameterization as in the original
psignifit software (Hill, 2001):

F ( x | a,b ) = 1. / ( 1 + exp ( - (x-a)/b ) )

This is the default setting.

For a 2AFC task, the guessing rate is fixed at 0.5. Thus, our model has three free parameters:
a, b, and lapse. We want to keep a and b unconstrained and restrict lapse to values between
0 and 0.1:

>>> constraints = ( 'unconstrained', 'unconstrained', 'Uniform(0,0.1)' )

Now we can fit the psychometric function

>>> B = BootstrapInference ( data, priors=constraints )
>>> print B
< BootstrapInference object with 6 blocks and 0 samples >
>>> B.estimate
array([ 2.7517686 ,  1.45723724,  0.01555636])
>>> B.deviance
8.071331367479198

Thus, the a is approximately 2.75, b is approximately 1.46, and lambda is approximately 0.016.
How well do these parameters describe the data? The deviance is approximately 8.07. Is this a
high or a low value? To know this, we have to draw a number of bootstrap samples:

>>> B.sample()
>>> print B
< BootstrapInference object with 6 blocks and 2000 samples >

We see that B has changed: instead of 0 samples, we now have 2000 parametric bootstrap samples
in the object. We can use these samples to assess the goodness of fit

>>> GoodnessOfFit(B)

in an interactive session, this should open a window that looks like the following. (In some
cases, you may have to type show() before you see the window).

.. image:: BootstrapGoodnessOfFit.png

The panel in the upper left displays the fitted psychometric function and the data points.
In addition, some information about the fitted model is displayed and confidence intervals for
thresholds at three levels are shown. The panel in the lower left displays a histogram of
the deviances that were to be expected if the fitted model was perfectly correct. In addition,
there are 95% confidence limits (dotted lines) and the observed deviance. If the observed
deviance is outside th 95% confidence limits, this is an indication of a bad fit. The plot
in the middle on top plot deviance residuals against the predicted correct response rate of
the model. This plot helps to detect systematic deviations between model and data. Trends in
this graph indicate systematic differences between model and data. The dotted line is the
best linear fit that relates deviance residuals to the predicted correct response rate.
The correlation between model prediction and deviance residuals is shown in the plot. The
plot in the middle at the bottom shows a histogram of these correlations under the assumption
that the fitted model is perfectly correct. Again the dotted lines denote 95% intervals
of the correlations and the solid line marks the observed correlation between model prediction
and deviance residuals. If the obsered correlation is outside the 95% interval, this indicates
a systematic deviation of the model from the data. The two plots on the right follow the same
logic as the plots in the middle. The difference is that in this case deviance residuals are
plotted agains block index, i.e. the sequence in which the data were acquired. If the observer
was still learning the task, this should be visible in this plot.

We can also get a graphical representation of the fitted parameters:

>>> ParameterPlot(B)

in an interactive session, this should again open a window showing estimated densities of
the model parameters as shown below. (Again, you might have to type show() to see the window).

.. image:: BootstrapParameters.png

Each of these plots shows the estimated density of one of the model parameters. In addition,
the estimated parameter is marked by a solid vertical line and the 95% confidence interval is
marked by dotted vertical lines. The confidence interval limits and the estimates are written
on top of the graph.

Reparameterizing the model
--------------------------

pypsignifit reformulates the function F ( x | a,b ) by means of two separate functions f: \R->\R
and g: \R^3->\R. We can think of f as the nonlinear part of the psychometric function, while
g is in most cases linear in x. Often g can be changed without seriously altering the possible
model shapes. In pypsignifit f is called the 'sigmoid' and g is called the 'core'. Using different
combinations of sigmoid and core allows a high flexibility of model fitting. For instance
Kuss, et al (2005) used a parameterization in terms of the 'midpoint' m of the sigmoid and the
'width' w. Here width is defined as the distance F^(-1) ( 1-alpha ) - F^(-1) ( alpha ). To
perform BootstrapInference for this model we can proceed as follows

>>> Bmw = BootstrapInference ( data, sample=2000, priors=constraints, core="mw0.1" )

# TODO: This is not finished yet and it does not work! Deviances are far too high

Example 2: Bayesian inference
=============================

References
==========
Hill, NJ (2001): Testing Hypotheses About Psychometric Functions. PhD Thesis, Oxford.
Kuss, M, Jäkel, F, Wichmann, FA (2005): Bayesian inference for psychometric functions. J Vis, 5,
    478-492.
Wichmann, FA, Hill, NJ (2001a): The psychometric function: I. Fitting, sampling, and goodness of fit.
    Perc Psychophys, 63(8), 1293-1313.
Wichmann, FA, Hill, NJ (2001b): The psychometric function: II. Bootstrap-based confidence intervals
    and sampling. Perc Psychophys, 63(8), 1314-1329.
