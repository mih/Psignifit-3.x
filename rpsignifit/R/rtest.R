x <- c(0,2,4,6,8,10)
k <- c(24,22,33,50,47,50)
n <- rep(50,length(k))
priors <- vector("character",3)
priors[1] <- "Gauss(0,10)"
priors[2] <- "Gamma(1,4)"
priors[3] <- "Beta(2,40)"

dyn.load("rpsignifit.so")
out <- .C( "mapestimate",
	stimulus.intensities=as.double(x),
	number.of.correct=as.integer(k),
	number.of.trials=as.integer(n),
	number.of.blocks=as.integer(length(k)),
	sigmoid=as.character("logistic"),
	core=as.character("mw0.1"),
	number.of.alternatives=as.integer(2),
	estimate=as.double(vector("numeric",3)),
	number.of.parameters=as.integer(3),
	deviance=as.double(0),
	priors=as.character(priors)
	)

print (out)
