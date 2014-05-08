## makeCluster
##
##' Create independant clusters from a sequence of (observed + random)
##' steps.
##'
##' There must be one and only one case per strata. In addition, the
##' sequences must be orderd in ascending order (there is no check on
##' this), and the case and strata should correspond to that order.
##'
##' It is only necessary to provide \code{nclust} to get at least that
##' many clusters (if possible). The argument \code{nloc} can be set
##' instead, if one wants exactly a number of successive steps in a
##' cluster.
##' @title Create independent clusters
##' @param seq A numeric, indicating the continuous sequences of
##' steps. See Details.
##' @param strata A numeric, indicating the strata of each observed
##' and random steps.
##' @param case A numeric, with 1 for observed steps and 0 for random
##' steps.
##' @param ndrop The number of steps to drop between each
##' cluster. They will be NAs in the resulting vector.
##' @param nclust The ideal number of clusters to obtain. See Details.
##' @param nstep The number of steps consisting in a single
##' cluster. See Details.
##' @return A vector of the same length as \code{seq}, \code{strata}
##' and \code{case}, with the cluster number for each step, and three
##' attributes giving the number of clusters (\code{nclust}), the
##' number of successive steps kept per cluster (\code{nstep}) and the
##' number of successive steps dropped between each cluster
##' (\code{ndrop}).
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## case (1 observed + 9 random)
##' case <- rep(rep(1:0, c(1, 9)), 100)
##' ## 100 stratas of 10 steps (1 observed + 9 random)
##' strata <- rep(1:100, each = 10)
##' ## 5 sequences of 20*10 steps
##' seq <- rep(1:5, each = 200)
##' head(data.frame(case, strata, seq), 22)
##' ## Make at least 15 clusters by dropping 2 steps between each of them
##' makeCluster(seq, strata, case, ndrop = 2, nclust = 15)
##' ## Same result with 'nstep = 7' direcly
##' makeCluster(seq, strata, case, ndrop = 2, nstep = 7)
makeCluster <- function(seq, strata, case, ndrop, nclust = NULL,
    nstep = NULL)
{
    ## Check that there is exactly 1 case for each strata, else stop
    tabcase <- table(case, strata)
    if (any(tabcase["1", ] != 1))
        stop("Some stratas with no or more than one case.")
    ## One and only one of 'nclust' or 'nstep' must be provided
    if ((is.null(nclust) & is.null(nstep)) | (!is.null(nclust) &
        !is.null(nstep)))
        stop("One, and only one, of 'nclust' or 'nstep' must be provided.")
    ## Sub-function to create the clusters
    cluster <- function(seq, strata, case, ndrop, nstep) {
        ## Create the pattern [(0*nstep),(NA*ndrop)]
        pat <- rep(c(0, NA), times = c(nstep, ndrop))
        ## Initiate a counter
        compt <- 0
        ## Initiate the 'cluster' vector
        clust <- numeric()
        ## For each sequence, do:
        for (i in unique(seq)) {
            ## How many steps (observed+random) per strata
            ## corresponding to this sequence?
            tabi <- table(strata[seq == i])
            ## How many clusters can fit in this sequence? ('ceiling'
            ## because it suffices that a cluster starts)
            maxi <- ceiling(sum(seq[case == 1] == i)/(nstep +
                ndrop))
            ## If more than one cluster, repeat the pattern, and set
            ## the corresponding cluster number to each kept step; cut
            ## it to the number of cases
            if (maxi > 1)
                clusti <- (rep(pat, times = maxi) + rep(compt +
                  1:maxi, each = (nstep + ndrop)))[1:sum(seq[case ==
                  1] == i)]
            ## If only one cluster, add 'maxi' (1) to the counter, and
            ## repeat it as many times as there are cases
            else clusti <- rep(compt + maxi, sum(seq[case ==
                1] == i))
            ## Repeat the cluster/case given the number of steps and
            ## add it to the already existing vector of clusters
            clust <- c(clust, rep(clusti, times = tabi))
            ## Increment the counter using the max number of clusters
            compt <- compt + maxi
        }
        ## Return the vector of clusters
        return(clust)
    }
    ## If 'clust' is not provided, use 'nstep' directly:
    if (is.null(nclust))
        clust <- cluster(seq = seq, strata = strata, case = case,
            ndrop = ndrop, nstep = nstep)
    ## Else iterate 'nstep' in 'cluster' until the number of cluster
    ## is lower than
    else {
        ## Automatic approximate 'nstep'
        nstep <- floor(sum(case)/(nclust)) - ndrop
        clust <- cluster(seq = seq, strata = strata, case = case,
            ndrop = ndrop, nstep = nstep)
        ## If the number of clusters is equal or greater than the
        ## required one, iterate through increasing 'nstep'
        if (max(clust, na.rm = TRUE) >= nclust) {
            for (i in (nstep + 1):sum(case)) {
                clusttmp <- cluster(seq = seq, strata = strata,
                  case = case, ndrop = ndrop, nstep = i)
                ## If the number of clusters is no more equal or
                ## greater than the required one, breaks
                if (max(clusttmp, na.rm = TRUE) < nclust) {
                  break
                }
                ## Else store 'clusttmp' and 'nstep' (and go to the
                ## next iteration)
                else {
                  clust <- clusttmp
                  nstep <- i
                }
            }
        }
        ## Else if the number of clusters is smaller than the required
        ## one, iterate through decreasing 'nstep'
        else {
            for (i in (nstep - 1):0) {
                clusttmp <- cluster(seq = seq, strata = strata,
                    case = case, ndrop = ndrop, nstep = i)
                ## If the number of clusters is no more smaller than
                ## the required one, breaks
                if (max(clusttmp, na.rm = TRUE) >= nclust) {
                  break
                }
                ## Else store 'clusttmp' and 'nstep' (and go to the
                ## next iteration)
                else {
                  clust <- clusttmp
                  nstep <- i
                }
            }
        }
        if (max(clust, na.rm = TRUE) < nclust)
            stop("Impossible to return that many clusters. Try with a smaller 'nclust'.")
    }
    ## Set 'nclust', 'nstep' and 'ndrop' attributes
    attr(clust, "nclust") <- max(clust, na.rm = TRUE)
    attr(clust, "nstep") <- nstep
    attr(clust, "ndrop") <- ndrop
    return(clust)
}
