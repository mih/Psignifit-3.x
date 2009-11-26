\name{bootstrap.goodness.of.fit}
\alias{bootstrap.gof}
\title{Goodness of fit diagnostic plot for bootstrap inference}
\description{This function creates a goodness of fit plot for bootstrap inference. A plot with 6 panels in two rows is generated. The top row
    shows information about the constrained maximum likelihood point estimate, the bottom row shows information about the distribution of
    bootstrap samples. Going through the panels by column we have the following information:
    1. The fitted psychometric function with the fitted data points. The size of the data points codes for the number of trials presented.
    In addition to the psychometric function, the bootstrap BCa confidence intervals are shown.
    2. A histogram of deviances expected for the fitted psychometric function. The deviance of the observed data is also marked. If this deviance is in the
    extremes of the expected deviances, this suggests that the fitted psychometric function is not a good description of the data.
    3. A scatter plot of deviance residuals against model prediction. If the fitted model describes all the structure in the data, there should be no
    systematic trend in these data. To aid intuition about trends, the best fitting line is plotted with the data points.
    4. A histogram of the correlations of deviance residuals and model prediction as expected if the fitted psychometric function was correct. The
    correlation for the observed data should not be in the extremes of this histogram. A 95\% confidence interval for the case where the fitted function
    is correct is marked by dotted lines and the correlation for the observed data is marked by a solid line.
    5. A scatter plot of deviance residuals against block index. This is an indicator of learning. If there is no learning, there should be no trend in this
    plot. The rest is analog to 3.
    6. This is analog to 4 but for correlations of deviance residuals and block index.
}
\usage{bootstrap.goodness.of.fit ( inference )}
\arguments{
\item{inference} A psiginference object as returned by PsigBootstrap
}
\references{
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293-1313.
}


\name{bayes.goodness.of.fit}
\alias{bayes.gof}
\title{Goodness of fit diagnostic plot for bayes inference}
\description{This function creates a goodness of fit plot for bayesian inference based on MCMC. A plot with 6 panels in two rows is generated. The top
    row shows information about the mean estimate, the bottom row shows information about selected properties of the posterior distribution.
    Going through the panels by column we have the following information:
    1. The fitted psychometric function (mean of the posterior) with the fitted data points. The size of the data points codes for the number of
    trials presented. In addition to the psychometric function two measures of spread of the posterior in function space are shown. (i) There are
    horizontal lines marking the 95\% posterior intervals of the respective performance levels. (ii) In light gray 20 posterior samples are displayed
    as psychometric functions.
    2. Posterior predictive plot of the deviance of the posterior. For each posterior parameter sample, a data set is generated. These data sets
    are called posterior predictive data sets. The plot shows deviances for the observed data set (for all posterior parameter samples) plotted
    against deviances for the posterior predictive data sets. If the observed data set is likely to be generated from the fitted model, the points
    in this plot should scatter around the diagonal. The diagonal is drawn as a dotted line in this plot.
    3. A scatter plot of deviance residuals against model prediction (for the mean of the posetrior). If the fitted model describes all the
    structure in the data, there should be no systematic trend in these data. To aid intuition about trends, the best fitting line is plotted with
    the data points.
    4. Posterior predictive plot of the correlation between deviance residuals and model prediction. This is similar to plot (2.) but for correlation
    betwenn deviance residuals and model prediction instead of deviance.
    5. A scatter plot of deviance residuals against block index. This is an indicator of learning. If there is no learning, there should be no trend in this
    plot. The rest is analog to 3.
    6. This is analog to 4 but for correlations of deviance residuals and block index.
}
\usage{bayes.goodness.of.fit ( inference )}
\arguments{
\item{inference} A psiginference object as returned by PsigBayes
}
\references{
    Wichmann, FA & Hill, NJ ( 2001 ): The psychometric function: I. Fitting, sampling, and goodness of fit. Perception & Psychophysics, 63(8), 1293-1313.
}


