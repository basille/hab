## makeSeq
##
##' Define sequences of a trajectory based on the chronology of steps,
##' allowing for gaps.
##'
##' @title Define sequences of a trajectory
##' @param x A ltraj object.
##' @param gap The maximum time interval between two successives steps
##' before starting a new sequence (note that the time interval needs
##' to be strictly greater than \code{gap} for a new sequence to
##' start).
##' @param units The unit of \code{gap} (default is \code{hour}).
##' @param na.omit Logical, whether to remove missing locations to
##' form the sequence (default).
##' @param name A character string indicating the column name of the
##' sequence to be stored in infolocs.
##' @return For \code{makeSeq}, a ltraj with a new variable
##' \code{name} in infolocs, with the sequence number for each
##' location, in the form of an integer series for each burst.
##' @author Mathieu Basille \email{basille@@ase-research.org}
##' @export
##' @examples
##' data(puechcirc)
##'
##' puechcirc <- makeSeq(puechcirc, gap = 15, units = "mins", name = "seq15")
##' infolocs(puechcirc, which = "seq15", simplify = TRUE)
##' summarySeq(puechcirc, "seq15")
##'
##' puechcirc <- makeSeq(puechcirc, gap = 25, units = "mins", name = "seq25")
##' summarySeq(puechcirc, "seq25")
makeSeq <- function(x, gap = 1, units = c("hours", "mins", "secs",
    "days", "weeks"), na.omit = TRUE, name = "seq")
{
    ## Check that x is a ltraj object
    if (!inherits(x, "ltraj"))
        stop("`x` should be of class `ltraj`")
    ## Match `units`
    units <- match.arg(units)
    ## Function for one burst
    lseq <- function(df, gap, units, na.omit) {
        ## Get date
        ldate <- df[["date"]]
        ## Remove missing locs
        if (na.omit) {
            nas <- !is.na(df[["x"]])
            ldate <- subset(ldate, nas)
        }
        ## Compute the time intervals between locs, in unit `units`
        ddiff <- difftime(ldate[-1], ldate[-length(ldate)],
            units = units)
        ## Check that the `date`` is correctly ordered
        if (any(ddiff < 0))
            stop("The `date` must be chronologically ordered.")
        ## The sequence starts with 1, and is the cumulative sum of the
        ## intervals greater than `gap`
        seq <- cumsum(c(1, ddiff > gap))
        ## Add NAs for missing locs
        if (na.omit) {
            tmp <- rep(NA, nrow(df))
            tmp[which(nas)] <- seq
            seq <- tmp
        }
        ## Intervals between each sequence are given by `ddiff` greater
        ## than `gap`, and is stored in the attribute `inter`
        attr(seq, "inter") <- ddiff[ddiff > gap]
        ## The gap is stored in the attribute `gap`
        attr(seq, "gap") <- gap
        ## The unit for gap is stored in the attribute `gap-unit`
        attr(seq, "gap-unit") <- units
        ## Put everything in a data frame
        seq <- data.frame(seq, row.names = row.names(df))
        names(seq) <- name
        return(seq)
    }
    ## If there is no infolocs, just compute the seq and store it as
    ## the new infolocs
    if (is.null(infolocs(x)))
        infolocs(x) <- lapply(x, lseq, gap = gap, units = units,
            na.omit = na.omit)
    ## Otherwise, bind it to the existing infolocs using mapply
    else infolocs(x) <- mapply(cbind, infolocs(x), lapply(x,
        lseq, gap = gap, units = units, na.omit = na.omit),
        SIMPLIFY = FALSE)
    return(x)
}


## summarySeq
##
##' @rdname makeSeq
##' @param x A ltraj object.
##' @param name A character string indicating the column name of the
##' sequence computed by \code{makeSeq}.
##' @return For \code{summarySeq}, an object of class \code{summarySeq}.
##' @export
summarySeq <- function(x, name = "seq")
{
    ## Retrieve the sequence
    seq <- infolocs(x, name, simplify = TRUE)
    ## `summ` is a list giving the number of sequences (max seq
    ## #), the length of each sequence and the intervals between
    ## sequences, constrained by `id`
    summ <- lapply(seq, function(i) list(nr = max(i, na.rm = TRUE),
        length = unclass(table(i)), inter = attr(i, "inter")))
    ## Names of `summ` given by `id`
    names(summ) <- burst(x)
    ## `units` is stored in an attribute
    attr(summ, "units") <- attr(seq[[1]], "gap-unit")
    ## `summ` is of class `summarySeq`
    class(summ) <- "summarySeq"
    return(summ)
}


## print.summarySeq
##
##' @rdname makeSeq
##' @param x An object of class \code{summarySeq}.
##' @export
print.summarySeq <- function(x, ...)
{
    ## Check that `x` is of class `x`
    if (!inherits(x, "summarySeq"))
        stop("x should be of class `summarySeq`")
    ## Print the total number of sequences
    cat("Total number of sequences:", sum(unlist(lapply(x,
        function(i) i$nr))), "\n")
    ## Print the number of sequences per individual (mean, sd)
    cat("\nNumber of sequences per individual (mean = ",
        round(mean(unlist(lapply(x, function(i) i$nr)), na.rm = TRUE), 2),
        ", sd = ", round(sd(unlist(lapply(x, function(i) i$nr)),
          na.rm = TRUE), 2), "):\n", sep = "")
    print(unlist(lapply(x, function(i) i$nr)))
    ## Print a summary of sequence lengths (all indidividuals combined)
    cat("\nLength of sequences (all individuals, steps):\n")
    ## Print a summary of sequence intervals (all indidividuals
    ## combined)
    print(summary(unlist(lapply(x, function(i) i$length))))
    cat("\nTime intervals between sequences (all individuals, ",
        attr(x, "units"), "):\n", sep = "")
    print(summary(unlist(lapply(x, function(i) i$inter))))
}
