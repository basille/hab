## lincircor
##
##' Compute Linear-circular correlation, and test it.
##'
##' @title Linear-circular correlation
##' @param xl A numeric, providing the linear variable.
##' @param theta A numeric, providing the circular variable, in
##' radians (\code{rad = deg * pi / 180}).
##' @return \code{lincircor} returns a numeric providing the
##' linear-circular correlation. \code{lincircor.test} returns a list
##' with class \code{"htest"} containing the following components:
##' \itemize{
##' \item \code{statistic}: the value of the test statistic
##' \item \code{parameter}: the degrees of freedom of the test
##' statistic, which follows a F distribution if \code{x}, and
##' \code{theta} are independent and \code{x} is normally distributed
##' \item \code{p.value}: the p-value of the test
##' \item \code{estimate}: the estimated measure of linear-circular
##' correlation
##' \item \code{null.value}: the value of the correlation measure
##' under the null hypothesis, always \code{0}
##' \item \code{alternative}: a character string describing the
##' alternative hypothesis
##' \item \code{method}: a character string indicating how the
##' correlation was measured
##' \item \code{data.name}: a character string giving the names of the
##' data}
##' @references Mardia, K. V. & Jupp, P. E. (2000) Directional
##' statistics. Wiley, 429 pp.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' ## From Mardia & Jupp (2000), p.246
##' xl <- c(28, 85.2, 80.5, 4.7, 45.9, 12.7, 72.5, 56.6, 31.5, 112,
##'     20, 72.5, 16, 45.9, 32.6, 56.6, 52.6, 91.8, 55.2)
##' theta <- c(327, 91, 88, 305, 344, 270, 67, 21, 281, 8, 204, 86,
##'     333, 18, 57, 6, 11, 27, 84)
##' ## 'xl' is not in radians (but the test is computed anyway)
##' lincircor(xl, theta)
##' ## Conversion to radians
##' lincircor(xl, theta*pi/180)
##' ## Test
##' lincircor.test(xl, theta*pi/180)
lincircor <- function(xl, theta) {
    ## Check the lengths of 'xl' and 'theta'
    if (length(xl) != length(theta))
        stop("x and theta should be vector of same length.")
    ## Quick check of 'theta' in radians (but calculations are done
    ## anyway)
    if (max(theta, na.rm = TRUE) > 2 * pi | min(theta, na.rm = TRUE) <
        -2 * pi)
        warning("theta should be in radians. Calculations done anyway.")
    ## Remove NAs
    nna <- is.na(xl) | is.na(theta)
    if (sum(nna) > 0)
        warning(paste(sum(nna), "NA(s) removed from the data."))
    xl <- xl[!nna]
    theta <- theta[!nna]
    ## Correlation arguments
    rxc <- cor(xl, cos(theta))
    rxs <- cor(xl, sin(theta))
    rcs <- cor(cos(theta), sin(theta))
    ## Return the linear-circular correlation
    return((rxc^2 + rxs^2 - 2 * rxc * rxs * rcs)/(1 - rcs^2))
}


## lincircor.test
##
##' @export
##' @rdname lincircor
lincircor.test <- function(xl, theta) {
    ## Computes the linear-circular correlation
    lincir <- lincircor(xl, theta)
    names(lincir) <- "cor"
    ## Df of the test
    df2 <- sum(!is.na(xl) & !is.na(theta)) - 3
    parm <- c(2, df2)
    names(parm) <- c("df1", "df2")
    ## Compute the test statistic
    stat <- df2 * lincir/(1 - lincir)
    names(stat) <- "stat"
    ## P-value of the test
    pv <- 1 - pf(stat, parm[1], parm[2])
    names(pv) <- "p.value"
    ## Null value of the test
    nv <- 0
    names(nv) <- "cor"
    ## Data names
    Call <- match.call()
    nm <- paste(as.character(Call[2]), as.character(Call[3]),
        sep = " and ")
    ## Make a 'htest' test (package stats) for pretty printing
    test <- list(statistic = stat, parameter = parm, p.value = pv,
        estimate = lincir, null.value = nv, method = "Linear-circular correlation",
        data.name = nm, alternative = "greater")
    class(test) <- "htest"
    return(test)
}
