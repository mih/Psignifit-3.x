x <- c(0,2,4,6,8,10)
k <- c(24,22,33,50,47,50)
n <- rep(50,length(k))
priors <- vector("character",3)
priors[1] <- "Gauss(0,10)"
priors[2] <- "Gamma(1,4)"
priors[3] <- "Beta(2,40)"

dyn.load("rpsignifit.so")
map <- .C( "mapestimate",
    stimulus.intensities=as.double(x),
    number.of.correct=as.integer(k),
    number.of.trials=as.integer(n),
    number.of.blocks=as.integer(length(k)),
    sigmoid=as.character("logistic"),
    core=as.character("mw0.1"),
    number.of.alternatives=as.integer(2),
    priors=as.character(priors),
    number.of.parameters=as.integer(3),
    estimate=as.double(vector("numeric",3)),
    deviance=as.double(0),
    Rpd=as.double(0),
    Rkd=as.double(0)
    )

boots <- .C ( "performbootstrap",
    stimulus.intensities=as.double(x),
    number.of.correct=as.integer(k),
    number.of.trials=as.integer(n),
    number.of.blocks=as.integer(length(k)),
    sigmoid=as.character("logistic"),
    core=as.character("mw0.1"),
    number.of.alternatives=as.integer(2),
    priors=as.character(priors),
    generating=as.double(map$estimate),
    number.of.parameters=as.integer(3),
    number.of.samples=as.integer(1000),
    cuts=as.double(0.5),
    number.of.cuts=as.integer(1),
    bootstrap.data=as.integer(vector("numeric",length(k)*1000)),
    bootstrap.estimates=as.double(vector("numeric",3*1000)),
    bootstrap.deviances=as.double(vector("numeric",1000)),
    bootstrap.Rpd=as.double(vector("numeric",1000)),
    bootstrap.Rkd=as.double(vector("numeric",1000)),
    bootstrap.thresholds=as.double(vector("numeric",1*1000)),
    influential=as.double(vector("numeric",length(k))),
    acceleration=as.double(vector("numeric",1)),
    bias=as.double(vector("numeric",1))
    )
bootstrap.data       <- matrix(boots$bootstrap.data,boots$number.of.samples,boots$number.of.blocks,TRUE)
bootstrap.estimate   <- matrix(boots$bootstrap.estimates,boots$number.of.samples,boots$number.of.parameters,TRUE)
bootstrap.thresholds <- matrix(boots$bootstrap.thresholds,boots$number.of.samples,boots$number.of.cuts,TRUE)

proposal <- vector("numeric",3)
proposal[1] = .4
proposal[2] = .4
proposal[3] = .01

mcmc <- .C ( "performmcmc",
    stimulus.intensities=as.double(x),
    number.of.correct=as.integer(k),
    number.of.trials=as.integer(n),
    number.of.blocks=as.integer(length(k)),
    sigmoid=as.character("logistic"),
    core=as.character("mw0.1"),
    number.of.alternatives=as.integer(2),
    priors=as.character(priors),
    proposal=as.double(proposal),
    start=as.double(map$estimate),
    number.of.parameters=as.integer(3),
    number.of.samples=as.integer(2000),
    cuts=as.double(0.5),
    number.of.cuts=as.integer(1),
    posterior.samples=as.double(vector("numeric",3*2000)),
    posterior.deviances=as.double(vector("numeric",2000)),
    predicted.data=as.integer(vector("numeric",length(k)*2000)),
    predicted.Rpd=as.double(vector("numeric",2000)),
    predicted.Rkd=as.double(vector("numeric",2000)),
    predicted.deviances=as.double(vector("numeric",2000)),
    logratios=as.double(vector("numeric",length(k)*2000))
    )

diag <- .C ( "getdiagnostics",
    stimulus.intensities=as.double(x),
    number.of.correct=as.integer(k),
    number.of.trials=as.integer(n),
    number.of.blocks=as.integer(length(k)),
    sigmoid=as.character("logistic"),
    core=as.character("mw0.1"),
    number.of.alternatives=as.integer(2),
    priors=as.character(priors),
    number.of.paramters=as.integer(3),
    cuts=as.double(0.5),
    number.of.cuts=as.integer(1),
    parameters=as.double(map$estimate),
    deviance=as.double(0),
    Rpd=as.double(0),
    Rkd=as.double(0),
    thres=as.double(vector("numeric",1)),
    deviance.residuals=as.double(vector("numeric",length(k)))
    )

xx <- seq(0,10,length.out=100)
Fx <- .C ( "pmfevaluate",
    stimulus.intensities=as.double(xx),
    number.of.intensities=as.integer(length(xx)),
    parameters=as.double(map$estimate),
    number.of.parameters=as.integer(3),
    sigmoid=as.character("logistic"),
    core=as.character("mw0.1"),
    number.of.alternatives=as.integer(2),
    f.x=as.double(vector("numeric",length(xx)))
    )

