==================================================
Specification of Models for Psychometric functions
==================================================

Typically, we want to use psignifit to fit data from a psychophysical experiment. These data
can come from different experimental designs and can display different shapes. These points are
not directly dependent on the actual fitting process.

Specifiing the experimental design
==================================

If nothing else is requested, psignifit will assume that the data come from a 2 alternatives
forced choice experiment. That means, that on each trial two stimuli were presented and the
observer knows that one and only one of these stimuli was the target stimulus. Although this
design has considerable theoretical advantages [Green_and_Swets]_ there might be practical
reasons to collect data in another way. One modification might be in showing more than two
stimuli. In this case, you should create your Inference object with the keyword 'nafc' set
to the number of stimulus alternatives that were presented. Setting the keyword 'nafc' to
a value of 2 or larger results in a fixed guessing rate of 1/nafc. These designs would all
inherit most of the theoretical advantages from the 2 alternatives forced choice experiment
and they would all be called forced choice experiments (with 4 alternatives, we would call
it a four alternatives forces choice experiment for example).

Another modification of the standard acquisition procedure could be to present only one
stimulus per trial. The observers might then have to indicate whether the target stimulus
was presented or not (typically called yes-no task). In some discrimination experiments,
the observers have to indicate whether a target stimulus was to the left or the right of
a mark or whether it was presented for a longer or shorter time than a reference stimulus.
In all these experiments we will record which response an observer chose and we will then
fit the number of "stimulus present", "stimulus left", "stimulus longer" responses (or
whatever is suitable in the present context). We will summarize these designs as "yes-no
designs" although the term yes-no is typically restricted to detection like tasks. The
crucial difference between yes-no designs and forced choice designs for fitting
psychometric functions is that yes-no designs allow for arbitrarily "guessing" rates.
For instance, in a detection task, the observer might be very conservative and virtually
never report the presence of a target for low stimulus intensities. Or the observer might
always respond "stimulus left" if the stimulus is presented sufficiently for to the left
of a mark. In all these situations, the lower asymptote of the psychometric function will
be a free parameter. As in all these situations only one stimulus is presented, we can
make the lower asymptote of the psychometric function a free parameter by setting the
keyword 'nafc' to 1.

A note on yes-no experiments on detection
-----------------------------------------

In a detection experiment, we typically have two types of trials. In one case a target
stimulus is presented, in the other case, no target stimulus is presented. In terms of signal
detection theory, these two cases are called "signal+noise" and "noise only". The observer
responds to these two cases with either "yes, a signal was present" or "no, only noise
was presented". This results in four different outcomes on each trial: hits (the observer
correctly reports the presence of a signal), false alarms (the observer reports a signal
although "noise only" was presented), misses (the observer reports no signal
although a signal was presented), and correct rejections (the observer correctly reports
the absence of a signal). These different experimental outcomes are discussed in large detail
in standard books on signal detection theory [Green_and_Swets]_. What matters with respect
to fitting psychometric functions is that depending on the strategy, both, hit rate as well
as the correct response rate might change. So which of these two should be fitted with a
psychometric function? There are two general objectives that an optimal observer could
follow in a yes-no task.

1. maximize the hit rate while keeping a fixed false alarm rate. In this case, we would
   like to fit the hit rate with a psychometric function (the false alarm rate is
   constant anyhow). Thus, if we want to fit the hit rate with a psychometric function,
   we should check that the observers maintained a more of less fixed false alarm
   rate. A future release of psignifit will contain more formal tools for this check.
2. maximize the number of correct responses. In this case, we would like to fit
   the correct response rate with a psychometric function. To check that an observer
   really uses this strategy, we should check that the false alarm rate decreases with
   the hit rate.

Signal detection theory also offers a number of criterion free discriminability parameters,
like area under the ROC curve and the famous d' index. However, these indices can not generally
be assumed to have binomial variance (or anything similar to that). Therefore, psignifit
does not attempt to fit such data.

Specifiing the shape of the psychometric function
=================================================


References
==========
.. [Green_and_Swets] Green, DM and Swets, JA (1966): Signal Detection Theory and Psychophysics. New York: Wiley.
