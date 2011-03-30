
==========================
Reparameterizing the model
==========================

pypsignifit reformulates the function :math:`F ( x | a,b )` by means of two separate functions :math:`f: \mathbb{R}\to\mathbb{R}`
and :math:`g: \mathbb{R}^3\to\mathbb{R}`. We can think of :math:`f` as the nonlinear part of the psychometric function, while
:math:`g` is in most cases linear in x. Often g can be changed without seriously altering the possible
model shapes. In pypsignifit :math:`f` is called the 'sigmoid' and :math:`g` is called the 'core'. Using different
combinations of sigmoid and core allows a high flexibility of model fitting. For instance
Kuss, et al (2005) used a parameterization in terms of the 'midpoint' :math:`m` of the sigmoid and the
'width' :math:`w`. Here width is defined as the distance :math:`F^{-1} ( 1-\alpha ) - F^{-1} ( \alpha )`. To
perform BootstrapInference for this model we can proceed as follows

>>> Bmw = BootstrapInference ( data, sample=2000, priors=constraints, core="mw0.1", nafc=nafc )
>>> Bmw.estimate
array([ 2.75176858,  6.40375494,  0.01555636])
>>> Bmw.deviance
8.0713313674704921
>>> Bmw.getThres()
2.7517685843037913
>>> Bmw.cuts
(0.25, 0.5, 0.75)
>>> Bmw.getCI(1)
array([ 1.4842732 ,  4.06407509])

Note that this model has the same deviance as the model fitted above. Also the obtained thresholds are the same.
However, as the parameterization is different, the actual fitted parameter values are different.
More details on sigmoids and cores and how they can be used to specify models can be found in the section
about _`Specification of Models for Psychometric functions`
