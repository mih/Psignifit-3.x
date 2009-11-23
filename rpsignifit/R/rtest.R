dyn.load("rpsignifit.so")

###################### PsigniSetup ###################################
# Basic data information and fitting directions
#
PsigniSetup <- function ( x, k, n,
    priors=list("","","Uniform(0,.1)"),
    sigmoid="logistic",
    core="mw0.1",
    number.of.alternatives=2 ) {

    data <- list ( stimulus.intensities=as.double(x),
        number.of.correct=as.integer(k),
        number.of.trials=as.integer(n),
        number.of.blocks=as.integer(length(k)),
        sigmoid=as.character(sigmoid),
        core=as.character(core),
        number.of.alternatives=as.integer(number.of.alternatives),
        priors=priors,
        cuts=0.5)
    attr(data,"class") <- "psignisetup"

    return (data)
}

print.psignisetup <- function ( data ) {
    nprm <- if (data$number.of.alternatives<2) 4 else 3
    parnames <- c("alpha","beta","lambda","gamma")
    if ( substr(data$core,1,2)=="mw" ) {
        parnames[1] = "m    "
        parnames[2] = "w    "
    }
    cat ( paste ("Data from ", data$number.of.alternatives, "-AFC experiment with ", data$number.of.blocks," blocks\n", sep="") )
    cat ( paste ("Intended fit is with",data$sigmoid,"sigmoid and",data$core,"core\n") )
    cat ( "Priors:\n" )
    for ( i in 1:nprm ) {
        cat ( paste ( "  ", parnames[i], "\t", priors[i], "\n" ) )
    }

    print ( data.frame ( stimulus.intensities=data$stimulus.intensities, number.of.correct=data$number.of.correct, number.of.trials=data$number.of.trials ) )
}

############################################################
#            psiginference methods                         #
############################################################



############################################################
#             psiginference types                          #
############################################################

MAPestimation <- function ( psignidata ) {
    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3
    map <- .C ( "mapestimate",
        stimulus.intensities=as.double(psignidata$stimulus.intensities),
        number.of.correct=as.integer(psignidata$number.of.correct),
        number.of.trials=as.integer(psignidata$number.of.trials),
        number.of.blocks=as.integer(psignidata$number.of.blocks),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(number.of.alternatives),
        priors=as.character(psignidata$priors),
        number.of.parameters=as.integer(nprm),
        estimate=as.double( vector("numeric", nprm) ),
        deviance=as.double(0),
        Rpd=as.double(0),
        Rkd=as.double(0)
        )
    map$data.samples <- NULL
    map$parameter.samples <- NULL
    map$deviance.samples <- NULL
    map$Rpd.samples <- NULL
    map$Rkd.samples <- NULL
    map$threshold.samples <- NULL
    map$influential <- NULL
    map$acceleration <- NULL
    map$bias <- NULL
    attr(map,"class") <- "psiginference"
    attr(map,"inference") <- "point"
    return (map)
}

PsigBootstrap <- function ( psignidata, number.of.samples=2000, generating=-999 ) {
    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3
    boots <- .C ( "performbootstrap",
        stimulus.intensities=as.double(psignidata$stimulus.intensities),
        number.of.correct=as.integer(psignidata$number.of.correct),
        number.of.trials=as.integer(psignidata$number.of.trials),
        number.of.blocks=as.integer(psignidata$number.of.blocks),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(psignidata$number.of.alternatives),
        priors=as.character(psignidata$priors),
        generating=as.double(generating),
        number.of.parameters=as.integer(nprm),
        number.of.samples=as.integer(number.of.samples),
        cuts=as.double(psignidata$cuts),
        number.of.cuts=as.integer(length(psignidata$cuts)),
        data.samples=as.integer(vector("numeric",psignidata$number.of.blocks*number.of.samples)),
        parameter.samples=as.double(vector("numeric",nprm*number.of.samples)),
        deviance.samples=as.double(vector("numeric",number.of.samples)),
        Rpd.samples=as.double(vector("numeric",number.of.samples)),
        Rkd.samples=as.double(vector("numeric",number.of.samples)),
        threshold.samples=as.double(vector("numeric",length(psignidata$cuts)*number.of.samples)),
        influential=as.double(vector("numeric",psignidata$number.of.blocks)),
        acceleration=as.double(vector("numeric",length(psignidata$cuts))),
        bias=as.double(vector("numeric",1),length(psignidata$cuts))
        )
    boots$data.samples <- matrix(boots$data.samples, boots$number.of.samples,boots$number.of.blocks, TRUE)
    boots$parameter.samples <- matrix(boots$parameter.samples, boots$number.of.samples, boots$number.of.parameters, TRUE)
    boots$threshold.samples <- matrix(boots$threshold.samples, boots$number.of.samples, boots$number.of.cuts, TRUE)

    map <- MAPestimation ( psignidata )
    boots$estimate <- map$estimate
    boots$deviance <- map$deviance
    boots$Rpd <- map$Rpd
    boots$Rkd <- map$Rkd
    boots$logratios <- NULL
    boots$proposal <- NULL
    boots$deviance.predictions <- NULL
    attr(boots,"class") <- "psiginference"
    attr(boots,"inference") <- "bootstrap"
    return (boots)
}

PsigBayes <- function ( psignidata, number.of.samples=2000, start=NULL, proposal=NULL ) {
    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3
    if (proposal==NULL) {
        # Use Raftery/Lewis instead? Would have to be implemented...
        proposal <- if (nprm==4) c(.4,.4,.01,.01) else c(.4,.4,.01)
    } else {
        if (length(proposal<nprm) {
            cat ("Error in PsigBayes: Wrong length of proposal argument\n")
            return (NULL)
        }
    }

    if (start==NULL) {
        start <- MAPestimation$estimate
    }

    mcmc <- .C ( "performmcmc",
        stimulus.intensities=as.double(psignidata$stimulus.intensities),
        number.of.correct=as.integer(psignidata$number.of.correct),
        number.of.trials=as.integer(psignidata$number.of.trials),
        number.of.blocks=as.integer(psignidata$number.of.blocks),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(psignidata$number.of.alternatives),
        priors=as.character(psignidata$priors),
        proposal=as.double(proposal),
        generating=as.double(start),
        number.of.parameters=as.integer(nprm),
        number.of.samples=as.integer(number.of.samples),
        cuts=as.double(psignidata$cuts),
        number.of.cuts=as.integer(length(psignidata$cuts)),
        parameter.samples=as.double(vector("numeric",nprm*number.of.samples)),
        deviance.samples=as.double(vector("numeric",number.of.samples)),
        data.samples=as.integer(vector("numeric",number.of.blocks*number.of.samples),
        Rpd.samples=as.double(vector("numeric",number.of.samples)),
        Rkd.samples=as.double(vector("numeric",number.of.samples)),
        deviance.predictions=as.double(vector("numeric",number.of.samples)),
        logratios=as.double(vector("numeric",number.of.blocks*number.of.samples))
        )
    mcmc$data.samples <- matrix(mcmc$data.samples, mcmc$number.of.samples, mcmc$number.of.blocks, TRUE)
    mcmc$parameter.samples <- matrix(mcmc$parameter.samples, mcmc$number.of.samples, mcmc$number.of.blocks, TRUE)

# TODO: MCMC mit thresholds?
    mcmc$threshold.samples <- NULL

    mcmc$logratios <- matrix(mcmc$logratios, mcmc$number.of.samples, mcmc$number.of.blocks, TRUE)

    return (mcmc)
}

PsigDiagnostics <- function ( parameters, psignidata ) {
    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    diag <- .C ( "getdiagnostics",
        stimulus.intensities=as.double(psignidata$stimulus.intensities),
        number.of.correct=as.integer(psignidata$number.of.correct),
        number.of.trials=as.integer(psignidata$number.of.trials),
        number.of.blocks=as.integer(psignidata$number.of.blocks),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(psignidata$number.of.alternatives),
        priors=as.character(psignidata$priors),
        number.of.parameters=as.integer(nprm),
        cuts=as.double(psignidata$cuts),
        number.of.cuts=as.integer(length(psignidata$cuts)),
        parameters=as.double(parameters),
        deviance=as.double(0),
        Rpd=as.double(0),
        Rkd=as.double(0),
        thres=as.double(vector("numeric",length(psignidata$cuts))),
        deviance.residuals=as.double(vector("numeric",psignidata$number.of.blocks))
        )
    attr(diag,"class") <- "psigdiagnostics"

    return (diag)
}

PsigEvaluate <- function ( parameters, psignidata, x=NULL ) {
    if (is.null(x)) x <- seq(min(psignidata$stimulus.intensities), max(psignidata$stimulus.intensities), length.out=100)
    nprm <- if (psignidata$number.of.alternatives<2) 4 else 3

    Fx <- .C ( "pmfevaluate",
        stimulus.intensities=as.double(x),
        number.of.intensities=as.integer(length(x)),
        parameters=as.double(parameters),
        number.of.parameters=as.integer(nprm),
        sigmoid=as.character(psignidata$sigmoid),
        core=as.character(psignidata$core),
        number.of.alternatives=as.integer(psignidata$number.of.alternatives),
        f.x=as.double(vector("numeric",length(x)))
        )

    return (Fx$f.x)
}

