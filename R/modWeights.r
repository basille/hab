## modWeights
##
##' Compute model weights using Information Theory criteria.
##'
##' @title Model weights using Information Theory
##' @param ... All fitted model objects.
##' @param criterion The function to be used as the Information Theory
##' criterion. Only tested with \code{QIC} (default).
##' @param names Logical, whether to assign names to the result.
##' @return A numeric vector with the corresponding weights.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
modWeights <- function(..., criterion = QIC, names = TRUE)
{
    ## Compute the criterion of choice on each model
    CritMod <- unlist(lapply(list(...), criterion))
    ## Compute the difference to the minimum criterion
    dCritMod <- CritMod - min(CritMod)
    ## Compute the weights
    weights <- exp(-0.5 * dCritMod)/sum(exp(-0.5 * dCritMod))
    ## If 'names', assign the names of the models
    if (names)
        names(weights) <- as.character(match.call()[-1L])
    ## Return the weights
    return(weights)
}
