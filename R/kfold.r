## kfold
##
##' Cross-validation for regression models.
##'
##' Note: needs complete names in the coxph/clogit call, and not
##' 'cos(var)' or 'I(var*10)', except for 'strata()' and
##' 'cluster()'.
##'
##' Also needs complete case in the model (e.g. it fails if
##' there is at least one NA in the observed steps for any variable of
##' the model). Returns an error if some stratas have no case.
##' @title kfold
##' @param mod A fitted model for which there exists a \code{kfold}
##' method (currently \code{coxph} and \code{clogit} models).
##' @param k The number of equal size subsamples of the partition.
##' @param nrepet The number of repetitions.
##' @param jitter Logical, whether to add some random noise to the
##' predictions (useful when the model is fitted on categorical
##' variables, which can produces error in the ranking process).
##' @param reproducible Logical, whether to use a fixed seed for each
##' repetition.
##' @param details Logical, whether to return details of each
##' repetition (useful for debugging).
##' @return A data frame with the correlations (\code{cor}) and the
##' type of value (\code{type}).
##' @author Mathieu Basille \email{basille@@ase-research.org}, with
##' the help of Terry Therneau and Guillaume Bastille-Rousseau
##' @export
kfold <- function(mod, k = 5, nrepet = 100, jitter = FALSE,
    reproducible = TRUE, details = FALSE)
{
    UseMethod("kfold")
}


##' @rdname kfold
kfold.coxph <- function(mod, k = 5, nrepet = 100, jitter = FALSE,
    reproducible = TRUE, details = FALSE)
{
    ## Check the class of the model (should be "coxph" or "clogit")
    if (!inherits(mod, "coxph"))
        stop("Model of class 'coxph' expected.")
    ## Load survival
    require(survival)
    ## Try to retrieve the data
    dt <- try(model.frame(mod), silent = TRUE)
    ## If it failed, stop and give a solution
    if (class(dt) == "try-error")
        stop("'model.frame' was unable to retrieve the data.",
          "Use 'model = TRUE' in the 'coxph' or 'clogit' call.")
    ## The first column is named 'srv' instead of 'Surv(faketime,
    ## case)'
    names(dt)[1] <- "srv"
    ## Which column is the strata?
    nstr <- attr(terms(mod), "specials")$strata
    ## Ugly regexp to extract and apply the strata variable name
    names(dt)[nstr] <- namestr <- sub("strata\\((.*)\\)", "\\1",
        names(dt)[nstr])
    ## If there is a cluster term...
    if (!is.null(attr(terms(mod), "specials")$cluster)) {
        ## Which column is the cluster?
        nclu <- attr(terms(mod), "specials")$cluster
        ## Ugly regexp to extract and apply the cluster variable name
        names(dt)[nclu] <- sub("cluster\\((.*)\\)", "\\1",
            names(dt)[nclu])
    }
    ## Is it really a problem?
    ## ncase <- table(tapply(dt$srv[, 2], dt[, nstr], function(x) sum(x == 1)))
    ## if (any(names(ncase) == "0"))
    ##     stop(paste("Some stratas had no case.",
    ##       "It is likely that NAs were present in the variables for some cases."))
    ## Prepare the 'kfold', 'rd' and 'warn' objects
    kfold <- rd <- warn <- numeric(length = nrepet)
    ## 'dbg' object for debugging when 'details = TRUE'
    if (details)
        dbg <- list()
    ## The core of the kfold, each repetition
    for (i in 1:nrepet) {
        ## Create a 'set' column, which defaults to "train"
        dt$sets <- "train"
        ## Allows for reproducibility
        if (reproducible)
            set.seed(i)
        ## Sample the "test" data set
        dt$sets[dt[, namestr] %in% sample(unique(dt[, namestr]),
            length(unique(dt[, namestr]))/k)] <- "test"
        ## Update the regression using the training data
        reg <- update(mod, srv ~ ., data = subset(dt, sets ==
            "train"), model = TRUE)
        ## Extract the "test" data set
        dtest <- droplevels(subset(dt, sets == "test"))
        ## And compute the predictions associated to this data set
        ## using the training regression
        dtest$predall <- exp(predict(reg, type = "lp", newdata = dtest,
            reference = "sample"))
        ## In case of equality among predictions (e.g. categorical
        ## variable), add some noise to the predictions
        if (jitter) {
            ## Allows for reproducibility
            if (reproducible)
                set.seed(i)
            dtest$predall <- jitter(dtest$predall)
        }
        ## The function to compute the rank within a strata
        samplepred <- function(df) {
            ## Number of controls
            nrand <- sum(df$srv[, 2] == 0)
            ## Rank of the case (among case + controls)
            obs <- rank(df$predall)[df$srv[, 2] == 1]
            ## Rank of a random control (among controls only!)
            if (reproducible)
                set.seed(i)
            rand <- sample(rank(df$predall[df$srv[, 2] == 0]),
                1)
            return(data.frame(obs = obs, rand = rand, nrand = nrand))
        }
        ## Compute the ranks for each strata and bind them together
        ranks <- do.call(rbind, by(dtest, dtest[, namestr], samplepred))
        ## Is there the same number of controls per strata?
        nrand <- unique(ranks$nrand)
        ## If no, use the greatest number of controls (and keep track
        ## of it)
        if (length(nrand) != 1) {
            nrand <- max(nrand)
            warn[i] <- 1
        }
        ## Compute the Spearman correlation on the ranks for the cases
        kfold[i] <- cor(1:(nrand+1), table(factor(ranks$obs,
            levels = 1:(nrand+1))), method = "spearman")
        ## Same for the random controls
        rd[i] <- cor(1:(nrand), table(factor(ranks$rand, levels = 1:(nrand))),
            method = "spearman")
        ## Store the ranks for debugging if 'details'
        if (details)
            dbg[[i]] <- ranks
    }
    ## Create a data frame with the correlations and the type of value
    res <- data.frame(cor = c(kfold, rd), type = rep(c("obs",
        "rand"), each = nrepet))
    ## Return the debugging information as an attribute if 'details'
    if (details)
        attr(res, "details") <- dbg
    ## If the number of controls is not the same for each strata, send
    ## a warning
    if (sum(warn) > 0)
        warning(paste("The number of controls was not equal among stratas for",
          sum(warn), "repetitions. Correlations might be biased.",
          "Use 'details = TRUE' to get more details."))
    ## Return the result data frame
    return(res)
}
