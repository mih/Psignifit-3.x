\name{PsigniSetup}
\alias{PsigniSetup}
\title{Data setup for fitting psychometric functions}
\description{
    In order to fit psychometric functions, there are a number of that need to be fixed.
    This function returns a list that contains everything that needs to be known to fit a psychometric
    function.
}
\usage{
    PsigniSetup ( x, k, n, priors=list("","","Uniform(0,.1)"), sigmoid="logistic", core="mw0.1", number.of.alternatives=2, cuts=0.5 )
}
\arguments{
\item{x} Stimulus intensities at which stimuli were presented. This should be a sequence in the order in which the stimuli were presented
    during the experiment.
\item{k} Number of correct responses (in nAFC trials) or number of YES responses (in Yes/No trials). These numbers should correspond to the
    entries in x.
\item{n} Number of repetitions of the respective stimuli.
\item{priors} a list of priors for the respective parameters. Psychometric functions for nAFC trials are governed by three parameters, the first two
    determine the shape of the sigmoid function and the last determines the upper asympote. Psychometric functions for Yes/No trials are described by
    four parameters. These are the same as for the nAFC case and an additional parameter for the lower asympote.
    Valid priors are "Gauss(mu,sigma)", "Gamma(
}
