## rdSteps
##
##' Draw random steps from a ltraj object using empirical distributions of
##' step lengths and turning angles.
##'
##' Note that 1) only complete steps are kept (i.e. steps characterized by
##' start and end points, and turning angle (i.e. relative angle to the
##' previous step); and 2) the information stored in infolocs is transfered
##' to all random steps within a strata (i.e. it assumes it is the same for
##' all steps within a strata).
##' @title Draw random steps
##' @param x A ltraj object.
##' @param nrs The number of random steps to draw for each observed step
##' (default = 10).
##' @param distMax Only draw step lengths and turning angles using steps
##' shorter than this threshold. Default is \code{Inf}, i.e. all steps are
##' kept.
##' @param only.others Logical, draws step lengths and turning angles from
##' all other individuals, excluding the current one.
##' @param simult Logical, whether to draw step lengths and turning angles
##' simultaneously, i.e. both measurements come from a single observed
##' step, instead of being drawn independently.
##' @param reproducible Logical,
##' @return A data frame, with new columns \code{case} (1 for observed
##' steps and 0 for random steps) and \code{strata} (an common integer for
##' paired observed and random steps).
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## Load the data
##' data(puechcirc)
##'
##' ## Simple example to check the distributions of step lengths and turning
##' ## angles
##' bla <- rdSteps(puechcirc)
##' boxplot(bla$rel.angle ~ bla$case)
##' boxplot(bla$dist ~ bla$case)
##'
##' ## Reproducibility
##' head(rdSteps(puechcirc, reproducible = TRUE))
##' head(rdSteps(puechcirc, reproducible = TRUE))
rdSteps <- function(x, nrs = 10, distMax = Inf, others = FALSE,
    simult = FALSE, reproducible = FALSE)
{
    ## Check if ltraj
    if (!inherits(x, "ltraj"))
        stop("x should be an object of class ltraj")
    ## Convert the ltraj to a data frame, and only keep complete steps
    ## (shortcut for na.omit(complete.steps = TRUE))
    xdf <- subset(ld(x), !(is.na(x) | is.na(dx) | is.na(rel.angle)))
    ## Stop if 'case' already exist
    if (!is.null(xdf$case))
        stop("The variable 'case' already exists in the ltraj.")
    ## case = 1 for observed steps
    xdf$case <- 1
    ## Stop if 'strata' already exist
    if (!is.null(xdf$strata))
        stop("The variable 'strata' already exists in the ltraj.")
    ## strata = 1:nrow(xdf) for observed steps
    xdf$strata <- 1:nrow(xdf)
    ## Prepare the empirical distributions of step lengths, turning angles,
    ## and the associated IDs (note that we only keep steps <= distMax)
    empdis <- subset(xdf, dist <= distMax, select = c("dist",
        "rel.angle", "id"))
    ## Split the data frame according to the ID
    xl <- split(xdf, xdf$id)
    ## Function to create a strata, i.e. random steps associated to a
    ## single observed step
    rdStep <- function(st, dis, seed) {
        ## Reproducible?
        if (reproducible)
            set.seed(seed)
        ## Step lengths and turning angles sampled simultaneously
        if (simult) {
            sp <- sample(1:nrow(dis), nrs, replace = TRUE)
            slrd <- dis$dist[sp]
            rhord <- dis$rel.angle[sp]
        }
        ## Not simultaneously
        else {
            ## Sample step lengths
            slrd <- sample(dis$dist, nrs, replace = TRUE)
            ## Sample turning angles
            rhord <- sample(dis$rel.angle, nrs, replace = TRUE)
        }
        ## Compute absolute angles
        alphard <- st$abs.angle - st$rel.angle + rhord
        ## Prepare the random data frame
        rd <- st[rep(1, nrs), ]
        ## Compute dx as cos(abs.angle)*step length
        rd$dx <- cos(alphard) * slrd
        ## Compute dy as sin(abs.angle)*step length
        rd$dy <- sin(alphard) * slrd
        ## Store step lengths, absolute angles and turning angles
        rd$dist <- slrd
        rd$abs.angle <- alphard
        rd$rel.angle <- rhord
        ## case = 0 for random steps
        rd$case <- 0
        ## Row names of the form 1.1, 1.2, etc.
        row.names(rd) <- paste(row.names(st), 1:nrs, sep = ".")
        ## Bind the step with the associated random steps
        return(rbind(st, rd))
    }
    ## Function that works over each id
    rdId <- function(lti) {
        ## Check id, and remove it from the distributions
        if (others)
            empdis <- subset(empdis, id != lti$id[1])
        ## Call rdStep on each line, bind the results together in a data
        ## frame
        return(do.call(rbind, lapply(1:nrow(lti), function(i)
            rdStep(lti[i, ], dis = empdis, seed = i))))
    }
    ## Call rdId over the ltraj, bind the results together in a data frame
    return(do.call(rbind, lapply(xl, rdId)))
}
