## QIC
##
##' Generic function calculating the Quasi-likelihood under Independence
##' Criterion for one or several fitted model objects.
##'
##' @title QIC: Quasi-likelihood under Independence Criterion
##' @param mod A \code{coxph} or \code{clogit} model.
##' @param ... Optionally more fitted model objects.
##' @return If \code{details = FALSE}, simply returns the QIC. If
##' \code{details = TRUE}, returns a data frame presenting the QIC,
##' the quasi-likelihood, the number of observations, the number of
##' events, the number of paramaters and the trace. If several models
##' are provided, two additional columns present the delta QIC and the
##' model weights.
##' @author Mathieu Basille \email{basille@@ase-research.org} and
##' Thierry Duchesne
##' @export
QIC <- function(mod, ...)
{
    UseMethod("QIC")
}

##' @rdname QIC
##' @param details Logical, whether to provide detailed outputs
##' (turned automatically to \code{TRUE} when several models are
##' fitted).
##' @export
QIC.coxph <- function(mod, ..., details = FALSE)
{
    ## If several models are provided
    if (!missing(...)) {
        ## Check if all models are CoxPH or a CLogit
        if (any(!sapply(list(mod, ...), inherits, "coxph")))
            stop("Object of class 'coxph' or 'clogit' expected.")
        ## Check if robust variances were estimated for all models
        if (any(!sapply(list(mod, ...), function(x) exists("naive.var",
            x))))
            stop("QIC can be computed only if robust variances are estimated.")
        ## Compute the trace term for all models
        trace <- sapply(list(mod, ...), function(x) sum(diag(solve(x$naive.var) %*%
            x$var)))
        ## Extract the quasi-likelihood for all models
        quasi <- sapply(list(mod, ...), function(x) x$loglik[2])
        ## Compute the QIC
        QIC <- -2 * quasi + 2 * trace
        ## Extract the number of observations for all models
        n <- sapply(list(mod, ...), function(x) x$n)
        ## Extract the number of events for all models
        nevent <- sapply(list(mod, ...), function(x) x$nevent)
        ## Extract the number of parameters for all models
        K <- sapply(list(mod, ...), function(x) length(x$coefficients))
        ## Check the number of observations and events
        if (any(n != n[1L]) | any(nevent != nevent[1L]))
            warning("Models are not all fitted to the same number of observations or events.")
        ## Prepare a detailed output
        val <- data.frame(QIC = QIC, QuasiLL = quasi, n = n,
            nevent = nevent, K = K, Trace = trace, deltaQIC = QIC -
                min(QIC), weight = modWeights(mod, ..., criterion = "QIC",
                names = FALSE))
        ## Assign the names of the models as row names
        row.names(val) <- as.character(match.call()[-1L])
        ## Return the data frame
        return(val)
    }
    else {
        ## Check if model is CoxPH or a CLogit
        if (!inherits(mod, "coxph"))
            stop("Object of class 'coxph' or 'clogit' expected.")
        ## Check if robust variances were estimated
        if (!exists("naive.var", mod))
            stop("QIC can be computed only if robust variances are estimated.")
        ## Compute the trace term
        trace <- sum(diag(solve(mod$naive.var) %*% mod$var))
        ## Extract the quasi-likelihood
        quasi <- mod$loglik[2]
        ## If 'details', return a detailed output
        if (details) {
            val <- data.frame(QIC = -2 * quasi + 2 * trace, QuasiLL = quasi,
                n = mod$n, nevent = mod$nevent, K = length(mod$coefficients),
                Trace = trace)
            ## Assign the name of the model as row names
            row.names(val) <- as.character(match.call()$mod)
            ## Return the data frame
            return(val)
        }
        ## Otherwise, just return the QIC
        else return(-2 * quasi + 2 * trace)
    }
}
